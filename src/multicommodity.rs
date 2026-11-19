//! Maximum concurrent fractional flow using a sparse edge-based linear program.

use std::error::Error;
use std::fmt;

use microlp::{
    ComparisonOp, LinearExpr, OptimizationDirection, Problem, SolutionStatus, SolveOptions,
    SolveOutcome, Variable,
};

use crate::{EdgeId, Graph, GraphError, GraphKind, NodeId};

const SOLVER_TOLERANCE: f64 = 1e-10;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Commodity {
    pub source: NodeId,
    pub sink: NodeId,
    pub demand: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FractionalEdgeFlow {
    pub edge: EdgeId,
    pub amount: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ConcurrentFlow {
    /// Each commodity ships this multiple of its demand. The optimum can exceed one.
    pub fraction: f64,
    /// Sparse edge flows in commodity order, with each inner list sorted by edge ID.
    pub commodity_flows: Vec<Vec<FractionalEdgeFlow>>,
}

#[derive(Debug)]
#[non_exhaustive]
pub enum ConcurrentFlowError {
    Graph(GraphError),
    NoCommodities,
    ZeroDemand { commodity: usize },
    Solver(microlp::Error),
    NumericalAccuracy,
}

impl fmt::Display for ConcurrentFlowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Graph(error) => error.fmt(f),
            Self::NoCommodities => f.write_str("at least one commodity is required"),
            Self::ZeroDemand { commodity } => write!(f, "commodity {commodity} has zero demand"),
            Self::Solver(error) => write!(f, "linear programming solver failed: {error}"),
            Self::NumericalAccuracy => {
                f.write_str("the solver did not produce a numerically valid optimum")
            }
        }
    }
}

impl Error for ConcurrentFlowError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Graph(error) => Some(error),
            Self::Solver(error) => Some(error),
            _ => None,
        }
    }
}

impl From<GraphError> for ConcurrentFlowError {
    fn from(error: GraphError) -> Self {
        Self::Graph(error)
    }
}

/// Maximizes the common demand multiplier subject to shared directed capacities.
/// Traffic may split across paths. Supplies, capacities, and the solution are
/// modeled in floating point. Results are checked for feasibility to a relative
/// tolerance of 1e-7. Ill-conditioned models may return a solver or accuracy error.
pub fn maximum_concurrent_flow(
    graph: &Graph<u64>,
    commodities: &[Commodity],
) -> Result<ConcurrentFlow, ConcurrentFlowError> {
    graph.require(GraphKind::Directed)?;
    if commodities.is_empty() {
        return Err(ConcurrentFlowError::NoCommodities);
    }
    let mut largest_demand = 0;
    for (index, commodity) in commodities.iter().enumerate() {
        graph.validate_node(commodity.source)?;
        graph.validate_node(commodity.sink)?;
        if commodity.source == commodity.sink {
            return Err(GraphError::IdenticalTerminals.into());
        }
        if commodity.demand == 0 {
            return Err(ConcurrentFlowError::ZeroDemand { commodity: index });
        }
        largest_demand = largest_demand.max(commodity.demand);
    }

    let demand_scale = largest_demand as f64;
    let mut problem = Problem::new(OptimizationDirection::Maximize);
    let throughput = problem.add_var(1.0, (0.0, f64::INFINITY));
    let mut capacity_rows = vec![Vec::new(); graph.edge_count()];
    let mut variables = Vec::with_capacity(commodities.len());
    let mut queue = Vec::new();
    let mut from_source = vec![false; graph.node_count()];
    let mut to_sink = from_source.clone();
    let mut edge_vars = vec![None; graph.edge_count()];

    for &commodity in commodities {
        reachable(graph, commodity.source, false, &mut from_source, &mut queue);
        if !from_source[commodity.sink] {
            return Ok(ConcurrentFlow {
                fraction: 0.0,
                commodity_flows: vec![Vec::new(); commodities.len()],
            });
        }
        reachable(graph, commodity.sink, true, &mut to_sink, &mut queue);
        edge_vars.fill(None);
        let mut entries = Vec::new();
        for (id, &edge) in graph.edges().iter().enumerate() {
            if edge.weight > 0
                && edge.from != edge.to
                && edge.to != commodity.source
                && edge.from != commodity.sink
                && from_source[edge.from]
                && to_sink[edge.to]
            {
                let variable = problem.add_var(0.0, (0.0, f64::INFINITY));
                edge_vars[id] = Some(variable);
                capacity_rows[id].push((variable, 1.0));
                entries.push((id, variable));
            }
        }
        add_balance_rows(
            &mut problem,
            graph,
            &edge_vars,
            commodity,
            throughput,
            demand_scale,
        );
        variables.push(entries);
    }
    for (id, row) in capacity_rows.into_iter().enumerate() {
        if !row.is_empty() {
            problem.add_constraint(row, ComparisonOp::Le, graph.edges()[id].weight as f64);
        }
    }

    let mut options = SolveOptions::default();
    options.tolerances.feasibility = SOLVER_TOLERANCE;
    let outcome = problem
        .solve_with(options)
        .map_err(ConcurrentFlowError::Solver)?;
    let SolveOutcome::Solution(solution) = outcome else {
        return Err(ConcurrentFlowError::NumericalAccuracy);
    };
    if solution.status() != SolutionStatus::Optimal {
        return Err(ConcurrentFlowError::NumericalAccuracy);
    }
    let fraction = solution.var_value(throughput) / demand_scale;
    if !fraction.is_finite() || fraction <= 0.0 {
        return Err(ConcurrentFlowError::NumericalAccuracy);
    }
    let commodity_flows = variables
        .into_iter()
        .map(|entries| {
            entries
                .into_iter()
                .filter_map(|(edge, variable)| {
                    let amount = solution.var_value(variable);
                    // Remove round-off below the lower bound before checking conservation.
                    (!(-SOLVER_TOLERANCE..=0.0).contains(&amount))
                        .then_some(FractionalEdgeFlow { edge, amount })
                })
                .collect()
        })
        .collect();
    let result = ConcurrentFlow {
        fraction,
        commodity_flows,
    };
    verify(graph, commodities, &result)?;
    Ok(result)
}

fn reachable(
    graph: &Graph<u64>,
    start: NodeId,
    reverse: bool,
    reached: &mut [bool],
    queue: &mut Vec<NodeId>,
) {
    reached.fill(false);
    reached[start] = true;
    queue.clear();
    queue.push(start);
    let mut head = 0;
    while let Some(&node) = queue.get(head) {
        head += 1;
        let arcs = if reverse {
            graph.inc(node)
        } else {
            graph.out(node)
        };
        for arc in arcs {
            if !reached[arc.node] && graph.edges()[arc.edge].weight > 0 {
                reached[arc.node] = true;
                queue.push(arc.node);
            }
        }
    }
}

fn add_balance_rows(
    problem: &mut Problem,
    graph: &Graph<u64>,
    variables: &[Option<Variable>],
    commodity: Commodity,
    throughput: Variable,
    demand_scale: f64,
) {
    for node in 0..graph.node_count() {
        // This conservation equation is implied by all the others.
        if node == commodity.sink {
            continue;
        }
        let out = graph
            .out(node)
            .iter()
            .filter_map(|arc| variables[arc.edge].map(|var| (var, 1.0)));
        let inc = graph
            .inc(node)
            .iter()
            .filter_map(|arc| variables[arc.edge].map(|var| (var, -1.0)));
        let mut terms = out.chain(inc).peekable();
        if terms.peek().is_none() {
            continue;
        }
        let mut row: LinearExpr = terms.collect();
        if node == commodity.source {
            row.add(throughput, -(commodity.demand as f64 / demand_scale));
        }
        problem.add_constraint(row, ComparisonOp::Eq, 0.0);
    }
}

fn verify(
    graph: &Graph<u64>,
    commodities: &[Commodity],
    result: &ConcurrentFlow,
) -> Result<(), ConcurrentFlowError> {
    const TOLERANCE: f64 = 1e-7;
    let mut loads = vec![0.0; graph.edge_count()];
    let mut balance = vec![(0.0_f64, 0.0_f64); graph.node_count()];
    for (commodity, edges) in commodities.iter().zip(&result.commodity_flows) {
        balance.fill((0.0, 0.0));
        for &FractionalEdgeFlow { edge, amount } in edges {
            if !amount.is_finite() || amount < 0.0 {
                return Err(ConcurrentFlowError::NumericalAccuracy);
            }
            let endpoints = graph.edges()[edge];
            loads[edge] += amount;
            balance[endpoints.from].0 += amount;
            balance[endpoints.to].1 += amount;
        }
        let expected = result.fraction * commodity.demand as f64;
        balance[commodity.source].1 += expected;
        balance[commodity.sink].0 += expected;
        if balance
            .iter()
            .any(|&(out, inc)| (out - inc).abs() > TOLERANCE * out.max(inc).max(expected))
        {
            return Err(ConcurrentFlowError::NumericalAccuracy);
        }
    }
    if graph
        .edges()
        .iter()
        .zip(loads)
        .any(|(edge, used)| !used.is_finite() || used > edge.weight as f64 * (1.0 + TOLERANCE))
    {
        return Err(ConcurrentFlowError::NumericalAccuracy);
    }
    Ok(())
}
