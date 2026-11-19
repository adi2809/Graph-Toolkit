//! Goldberg–Tarjan cost scaling with signed integer costs and finite capacities.

use std::collections::VecDeque;
use std::error::Error;
use std::fmt;

use super::{Residual, max_flow, residual_at, residual_neighbors};
use crate::{EdgeId, Graph, GraphBuilder, GraphError, GraphKind, NodeId};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ArcCost {
    pub capacity: u64,
    pub unit_cost: i64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CostFlow {
    pub total_cost: i128,
    pub edge_flows: Vec<u64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MinCostMaxFlow {
    pub value: u128,
    pub flow: CostFlow,
}

#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CostFlowError {
    Graph(GraphError),
    SupplyCount { expected: usize, actual: usize },
    EdgeFlowCount { expected: usize, actual: usize },
    CapacityExceeded { edge: EdgeId },
    ArithmeticOverflow,
    NegativeArcCost { edge: EdgeId },
}

impl fmt::Display for CostFlowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Graph(error) => error.fmt(f),
            Self::SupplyCount { expected, actual } => {
                write!(f, "expected {expected} vertex supplies, got {actual}")
            }
            Self::EdgeFlowCount { expected, actual } => {
                write!(f, "expected {expected} edge flows, got {actual}")
            }
            Self::CapacityExceeded { edge } => write!(f, "flow exceeds capacity on edge {edge}"),
            Self::NegativeArcCost { edge } => write!(f, "edge {edge} has a negative unit cost"),
            Self::ArithmeticOverflow => f.write_str("minimum-cost flow arithmetic exceeded i128"),
        }
    }
}

impl Error for CostFlowError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Graph(error) => Some(error),
            _ => None,
        }
    }
}

impl From<GraphError> for CostFlowError {
    fn from(error: GraphError) -> Self {
        Self::Graph(error)
    }
}

/// Cheapest directed flow meeting every vertex supply. Positive numbers send
/// flow, negative numbers receive it. Infeasible supplies yield `None`.
/// Passing all zero supplies computes a minimum-cost circulation.
pub fn minimum_cost_flow(
    graph: &Graph<ArcCost>,
    supply: &[i64],
) -> Result<Option<CostFlow>, CostFlowError> {
    feasible_amounts(graph, supply)?
        .map(|amounts| cheapest_with_same_balances(graph, amounts))
        .transpose()
}

pub(super) fn validate_supply(
    graph: &Graph<ArcCost>,
    supply: &[i64],
) -> Result<bool, CostFlowError> {
    graph.require(GraphKind::Directed)?;
    let n = graph.node_count();
    if supply.len() != n {
        return Err(CostFlowError::SupplyCount {
            expected: n,
            actual: supply.len(),
        });
    }

    Ok(supply.iter().copied().map(i128::from).sum::<i128>() == 0)
}

pub(super) fn feasible_amounts(
    graph: &Graph<ArcCost>,
    supply: &[i64],
) -> Result<Option<Vec<u64>>, CostFlowError> {
    if !validate_supply(graph, supply)? {
        return Ok(None);
    }
    let n = graph.node_count();
    let required: u128 = supply
        .iter()
        .filter(|&&amount| amount > 0)
        .map(|&amount| amount as u128)
        .sum();
    let amounts = if required == 0 {
        vec![0; graph.edge_count()]
    } else {
        let mut builder = capacities(graph, n + 2);
        for (node, &amount) in supply.iter().enumerate() {
            if amount > 0 {
                builder.add_edge(n, node, amount.unsigned_abs())?;
            } else if amount < 0 {
                builder.add_edge(node, n + 1, amount.unsigned_abs())?;
            }
        }
        let feasible = max_flow(&builder.build(), n, n + 1)?;
        if feasible.value != required {
            return Ok(None);
        }
        let mut amounts = feasible.edge_flows;
        amounts.truncate(graph.edge_count());
        amounts
    };
    Ok(Some(amounts))
}

/// Maximizes directed s–t flow value, then minimizes cost at that flow value.
pub fn min_cost_max_flow(
    graph: &Graph<ArcCost>,
    source: NodeId,
    sink: NodeId,
) -> Result<MinCostMaxFlow, CostFlowError> {
    graph.require(GraphKind::Directed)?;
    graph.validate_node(source)?;
    graph.validate_node(sink)?;
    if source == sink {
        return Err(GraphError::IdenticalTerminals.into());
    }
    let feasible = max_flow(&capacities(graph, graph.node_count()).build(), source, sink)?;
    let flow = cheapest_with_same_balances(graph, feasible.edge_flows)?;
    Ok(MinCostMaxFlow {
        value: feasible.value,
        flow,
    })
}

/// Computes cost for given edge flows. Checks length and capacities.
pub fn flow_cost(graph: &Graph<ArcCost>, amounts: &[u64]) -> Result<i128, CostFlowError> {
    if amounts.len() != graph.edge_count() {
        return Err(CostFlowError::EdgeFlowCount {
            expected: graph.edge_count(),
            actual: amounts.len(),
        });
    }
    for (id, (&amount, edge)) in amounts.iter().zip(graph.edges()).enumerate() {
        if amount > edge.weight.capacity {
            return Err(CostFlowError::CapacityExceeded { edge: id });
        }
    }

    // Consume an opposite-sign term when available, so cancellation cannot
    // overflow merely because many same-sign edges were inserted first.
    let terms = amounts
        .iter()
        .zip(graph.edges())
        .map(|(&amount, edge)| i128::from(amount) * i128::from(edge.weight.unit_cost));
    let mut positive = terms.clone().filter(|&cost| cost > 0);
    let mut negative = terms.filter(|&cost| cost < 0);
    let mut total = 0_i128;
    loop {
        let term = if total >= 0 {
            negative.next().or_else(|| positive.next())
        } else {
            positive.next().or_else(|| negative.next())
        };
        let Some(term) = term else { return Ok(total) };
        total = total
            .checked_add(term)
            .ok_or(CostFlowError::ArithmeticOverflow)?;
    }
}

fn capacities(graph: &Graph<ArcCost>, nodes: usize) -> GraphBuilder<u64> {
    let mut result =
        GraphBuilder::with_edge_capacity(nodes, GraphKind::Directed, graph.edge_count());
    for edge in graph.edges() {
        result
            .add_edge(edge.from, edge.to, edge.weight.capacity)
            .unwrap();
    }
    result
}

fn cheapest_with_same_balances(
    graph: &Graph<ArcCost>,
    amounts: Vec<u64>,
) -> Result<CostFlow, CostFlowError> {
    let n = graph.node_count();
    // At integer epsilon = 1, every simple residual cycle of length ≤ V
    // costs at least -V after scaling. Its cost is a multiple of V + 1.
    let factor = n as i128 + 1;
    let mut costs = Vec::with_capacity(graph.edge_count());
    let mut residual = Vec::with_capacity(graph.edge_count());
    let mut epsilon = 0;
    for (edge, mut amount) in graph.edges().iter().zip(amounts) {
        let cost = if edge.from == edge.to || edge.weight.capacity == 0 {
            if edge.from == edge.to && edge.weight.unit_cost != 0 {
                amount = if edge.weight.unit_cost < 0 {
                    edge.weight.capacity
                } else {
                    0
                };
            }
            0
        } else {
            i128::from(edge.weight.unit_cost) * factor
        };
        epsilon = epsilon.max(cost.abs());
        costs.push(cost);
        residual.push([edge.weight.capacity - amount, amount]);
    }
    let mut state = CostScaling {
        graph,
        costs,
        residual,
        prices: vec![0; n],
        excess: vec![0; n],
        current: vec![0; n],
        active: VecDeque::new(),
    };
    while epsilon > 1 {
        epsilon = (epsilon + 3) / 4;
        state.refine(epsilon)?;
    }
    let edge_flows: Vec<_> = state.residual.into_iter().map(|pair| pair[1]).collect();
    let total_cost = flow_cost(graph, &edge_flows)?;
    Ok(CostFlow {
        total_cost,
        edge_flows,
    })
}

struct CostScaling<'a> {
    graph: &'a Graph<ArcCost>,
    costs: Vec<i128>,
    residual: Vec<[u64; 2]>,
    prices: Vec<i128>,
    /// Inflow minus outflow relative to the input feasible flow.
    excess: Vec<i128>,
    current: Vec<usize>,
    active: VecDeque<NodeId>,
}

impl CostScaling<'_> {
    fn refine(&mut self, epsilon: i128) -> Result<(), CostFlowError> {
        self.current.fill(0);
        for (id, edge) in self.graph.edges().iter().enumerate() {
            if edge.from == edge.to || edge.weight.capacity == 0 {
                continue;
            }
            let forward = Residual {
                to: edge.to,
                edge: id,
                reversed: 0,
            };
            let cost = self.reduced_cost(edge.from, forward)?;
            if cost < 0 {
                self.push(edge.from, forward, self.residual[id][0]);
            } else if cost > 0 {
                let reverse = Residual {
                    to: edge.from,
                    edge: id,
                    reversed: 1,
                };
                self.push(edge.to, reverse, self.residual[id][1]);
            }
        }
        for (node, &amount) in self.excess.iter().enumerate() {
            if amount > 0 {
                self.active.push_back(node);
            }
        }

        while let Some(node) = self.active.pop_front() {
            while self.excess[node] > 0 {
                if let Some(arc) = residual_at(self.graph, node, self.current[node]) {
                    let capacity = self.residual[arc.edge][arc.reversed];
                    if arc.to != node && capacity > 0 && self.reduced_cost(node, arc)? < 0 {
                        let activated = self.excess[arc.to] <= 0;
                        let sent = self.excess[node].min(i128::from(capacity)) as u64;
                        self.push(node, arc, sent);
                        if activated && self.excess[arc.to] > 0 {
                            self.active.push_back(arc.to);
                        }
                    } else {
                        self.current[node] += 1;
                    }
                } else {
                    self.relabel(node, epsilon)?;
                }
            }
        }
        debug_assert!(self.excess.iter().all(|&amount| amount == 0));
        Ok(())
    }

    fn relabel(&mut self, node: NodeId, epsilon: i128) -> Result<(), CostFlowError> {
        let mut lowest = i128::MAX;
        for arc in residual_neighbors(self.graph, node) {
            if arc.to != node && self.residual[arc.edge][arc.reversed] > 0 {
                lowest = lowest.min(self.reduced_cost(node, arc)?);
            }
        }
        let adjustment = lowest
            .checked_add(epsilon)
            .ok_or(CostFlowError::ArithmeticOverflow)?;
        self.prices[node] = self.prices[node]
            .checked_sub(adjustment)
            .ok_or(CostFlowError::ArithmeticOverflow)?;
        self.current[node] = 0;
        Ok(())
    }

    fn reduced_cost(&self, from: NodeId, arc: Residual) -> Result<i128, CostFlowError> {
        let cost = self.costs[arc.edge];
        let cost = if arc.reversed == 0 { cost } else { -cost };
        self.prices[from]
            .checked_sub(self.prices[arc.to])
            .and_then(|difference| cost.checked_add(difference))
            .ok_or(CostFlowError::ArithmeticOverflow)
    }

    fn push(&mut self, node: NodeId, arc: Residual, sent: u64) {
        self.residual[arc.edge][arc.reversed] -= sent;
        self.residual[arc.edge][arc.reversed ^ 1] += sent;
        self.excess[node] -= i128::from(sent);
        self.excess[arc.to] += i128::from(sent);
    }
}
