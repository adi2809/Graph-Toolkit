use super::mean::mean_arc_ids;
use super::min_cost::{feasible_amounts, validate_supply};
use super::ratio::RatioArc;
use super::{ArcCost, CostFlow, CostFlowError, Residual, flow_cost, residual_neighbors};
use crate::{Graph, NONE, queue::IndexedHeap};

/// Goldberg–Tarjan minimum-mean residual cycle cancellation. Each pivot
/// saturates an arc on a cycle of minimum mean cost. The pivot selector
/// uses Karp's strongly connected component dynamic program.
pub fn minimum_cost_flow_minimum_mean_cycle(
    graph: &Graph<ArcCost>,
    supply: &[i64],
) -> Result<Option<CostFlow>, CostFlowError> {
    let Some(mut amounts) = feasible_amounts(graph, supply)? else {
        return Ok(None);
    };
    let mut arcs = Vec::new();
    let mut mapping = Vec::new();
    loop {
        arcs.clear();
        mapping.clear();
        for (id, (e, &flow)) in graph.edges().iter().zip(&amounts).enumerate() {
            let cost = i128::from(e.weight.unit_cost);
            for (u, v, cost, capacity, reverse) in [
                (e.from, e.to, cost, e.weight.capacity - flow, false),
                (e.to, e.from, -cost, flow, true),
            ] {
                if capacity != 0 {
                    arcs.push(RatioArc {
                        u,
                        v,
                        cost,
                        length: 1,
                    });
                    mapping.push((id, reverse, capacity));
                }
            }
        }
        let Some(cycle) = mean_arc_ids(graph.node_count(), &arcs, true)? else {
            break;
        };
        let delta = cycle.iter().map(|&i| mapping[i].2).min().unwrap();
        for i in cycle {
            let (id, reverse, _) = mapping[i];
            if reverse {
                amounts[id] -= delta;
            } else {
                amounts[id] += delta;
            }
        }
    }
    Ok(Some(CostFlow {
        total_cost: flow_cost(graph, &amounts)?,
        edge_flows: amounts,
    }))
}

struct Scaling<'a> {
    graph: &'a Graph<ArcCost>,
    flow: Vec<u64>,
    excess: Vec<i128>,
    potential: Vec<i128>,
    distance: Vec<i128>,
    previous: Vec<Option<(usize, Residual)>>,
    queue: IndexedHeap<i128>,
}

impl Scaling<'_> {
    fn capacity(&self, a: Residual) -> u64 {
        if a.reversed == 0 {
            self.graph.edges()[a.edge].weight.capacity - self.flow[a.edge]
        } else {
            self.flow[a.edge]
        }
    }
    fn price(&self, u: usize, a: Residual) -> Result<i128, CostFlowError> {
        let cost = i128::from(self.graph.edges()[a.edge].weight.unit_cost);
        let cost = if a.reversed == 0 { cost } else { -cost };
        self.potential[u]
            .checked_sub(self.potential[a.to])
            .and_then(|d| d.checked_add(cost))
            .ok_or(CostFlowError::ArithmeticOverflow)
    }
    fn push(&mut self, u: usize, a: Residual, amount: u64) {
        if a.reversed == 0 {
            self.flow[a.edge] += amount;
        } else {
            self.flow[a.edge] -= amount;
        }
        self.excess[u] -= i128::from(amount);
        self.excess[a.to] += i128::from(amount);
    }
    fn restore(&mut self, delta: u64) -> Result<(), CostFlowError> {
        for (id, e) in self.graph.edges().iter().enumerate() {
            let forward = Residual {
                edge: id,
                to: e.to,
                reversed: 0,
            };
            let (u, arc) = if self.price(e.from, forward)? < 0 {
                (e.from, forward)
            } else {
                let reverse = Residual {
                    edge: id,
                    to: e.from,
                    reversed: 1,
                };
                if self.price(e.to, reverse)? >= 0 {
                    continue;
                }
                (e.to, reverse)
            };
            let cap = self.capacity(arc);
            if cap >= delta {
                self.push(u, arc, cap);
            }
        }
        Ok(())
    }

    fn balance_path(&mut self, delta: u64) -> Result<bool, CostFlowError> {
        let threshold = i128::from(delta);
        self.distance.fill(i128::MAX);
        self.previous.fill(None);
        self.queue.clear();
        for (u, &amount) in self.excess.iter().enumerate() {
            if amount >= threshold {
                self.distance[u] = 0;
                self.queue.decrease(u, 0);
            }
        }
        let mut end = NONE;
        while let Some((d, u)) = self.queue.pop() {
            if self.excess[u] <= -threshold {
                end = u;
                break;
            }
            for a in residual_neighbors(self.graph, u) {
                if self.capacity(a) < delta {
                    continue;
                }
                let price = self.price(u, a)?;
                debug_assert!(price >= 0);
                let next = d
                    .checked_add(price)
                    .ok_or(CostFlowError::ArithmeticOverflow)?;
                if next < self.distance[a.to] {
                    self.distance[a.to] = next;
                    self.previous[a.to] = Some((u, a));
                    self.queue.decrease(a.to, next);
                }
            }
        }
        if end == NONE {
            return Ok(false);
        }
        // Truncated distance labels keep every eligible reduced cost
        // nonnegative, including arcs incident to vertices outside the search.
        let limit = self.distance[end];
        for (p, &d) in self.potential.iter_mut().zip(&self.distance) {
            *p = p
                .checked_add(d.min(limit) - limit)
                .ok_or(CostFlowError::ArithmeticOverflow)?;
        }
        let mut v = end;
        let mut amount = (-self.excess[end]).min(i128::from(u64::MAX)) as u64;
        while let Some((u, a)) = self.previous[v] {
            amount = amount.min(self.capacity(a));
            v = u;
        }
        amount = amount.min(self.excess[v].min(i128::from(u64::MAX)) as u64);
        amount -= amount % delta;
        v = end;
        while let Some((u, a)) = self.previous[v] {
            self.push(u, a, amount);
            v = u;
        }
        Ok(true)
    }
}

/// Capacity-scaling primal–dual transshipment for signed costs and finite
/// capacities. At unit scale the residual reduced costs certify optimality.
pub fn minimum_cost_flow_capacity_scaling(
    graph: &Graph<ArcCost>,
    supply: &[i64],
) -> Result<Option<CostFlow>, CostFlowError> {
    if !validate_supply(graph, supply)? {
        return Ok(None);
    }
    let excess = supply.iter().copied().map(i128::from).collect();
    scaled_amounts(graph, excess)?
        .map(|flow| priced_solution(graph, flow))
        .transpose()
}

/// Minimum-cost transshipment with a required minimum flow per original arc.
/// The edge payload supplies the upper bound and unit cost. A lower bound
/// above capacity is invalid. Incompatible vertex balances give Ok(None).
pub fn minimum_cost_flow_with_lower_bounds(
    graph: &Graph<ArcCost>,
    supply: &[i64],
    lower: &[u64],
) -> Result<Option<CostFlow>, CostFlowError> {
    let balanced = validate_supply(graph, supply)?;
    super::min_cost::validate_amounts(graph, lower)?;
    if !balanced {
        return Ok(None);
    }
    let mut excess: Vec<i128> = supply.iter().copied().map(i128::from).collect();
    let mut shifted = Graph::builder(graph.node_count(), graph.kind());
    for (e, &required) in graph.edges().iter().zip(lower) {
        excess[e.from] -= i128::from(required);
        excess[e.to] += i128::from(required);
        shifted.add_edge(
            e.from,
            e.to,
            ArcCost {
                capacity: e.weight.capacity - required,
                ..e.weight
            },
        )?;
    }
    let Some(mut amounts) = scaled_amounts(&shifted.build(), excess)? else {
        return Ok(None);
    };
    for (amount, &required) in amounts.iter_mut().zip(lower) {
        *amount += required;
    }
    Ok(Some(priced_solution(graph, amounts)?))
}

fn priced_solution(
    graph: &Graph<ArcCost>,
    edge_flows: Vec<u64>,
) -> Result<CostFlow, CostFlowError> {
    Ok(CostFlow {
        total_cost: flow_cost(graph, &edge_flows)?,
        edge_flows,
    })
}

fn scaled_amounts(
    graph: &Graph<ArcCost>,
    excess: Vec<i128>,
) -> Result<Option<Vec<u64>>, CostFlowError> {
    let n = graph.node_count();
    let mut state = Scaling {
        graph,
        flow: vec![0; graph.edge_count()],
        excess,
        potential: vec![0; n],
        distance: vec![0; n],
        previous: vec![None; n],
        queue: IndexedHeap::new(n),
    };
    if let Some(bit) = graph
        .edges()
        .iter()
        .map(|e| e.weight.capacity)
        .max()
        .and_then(u64::checked_ilog2)
    {
        for k in (0..=bit).rev() {
            let delta = 1u64 << k;
            state.restore(delta)?;
            while state.balance_path(delta)? {}
        }
    }
    if state.excess.iter().any(|&a| a != 0) {
        return Ok(None);
    }
    Ok(Some(state.flow))
}
