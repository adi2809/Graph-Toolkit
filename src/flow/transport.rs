use super::min_cost::{feasible_amounts, validate_supply};
use super::{
    ArcCost, CostFlow, CostFlowError, MinCostMaxFlow, Residual, flow_cost, residual_neighbors,
    validate_terminals,
};
use crate::queue::RadixHeap;
use crate::{Graph, GraphBuilder, GraphKind, NodeId};

/// Minimum-cost maximum s–t flow via successive shortest augmentations.
///
/// Original arcs of positive capacity must have nonnegative unit costs.
/// Reduced-cost Dijkstra handles their negative reverse residual arcs.
pub fn min_cost_max_flow_successive_shortest_path(
    graph: &Graph<ArcCost>,
    source: NodeId,
    sink: NodeId,
) -> Result<MinCostMaxFlow, CostFlowError> {
    validate_terminals(graph, source, sink)?;
    reject_negative_costs(graph)?;
    let (value, residual) = shortest_augmentations(graph, source, sink, u128::MAX)?;
    Ok(MinCostMaxFlow {
        value,
        flow: finish(graph, residual)?,
    })
}

/// Fixed-supply minimum-cost flow with nonnegative unit costs. Returns None
/// when all supplies cannot be met.
pub fn minimum_cost_flow_successive_shortest_path(
    graph: &Graph<ArcCost>,
    supply: &[i64],
) -> Result<Option<CostFlow>, CostFlowError> {
    if !validate_supply(graph, supply)? {
        return Ok(None);
    }
    reject_negative_costs(graph)?;
    let n = graph.node_count();
    let required: u128 = supply.iter().filter(|&&s| s > 0).map(|&s| s as u128).sum();
    let mut builder =
        GraphBuilder::with_edge_capacity(n + 2, GraphKind::Directed, graph.edge_count());
    for edge in graph.edges() {
        builder.add_edge(edge.from, edge.to, edge.weight).unwrap();
    }
    for (node, &amount) in supply.iter().enumerate() {
        let weight = ArcCost {
            capacity: amount.unsigned_abs(),
            unit_cost: 0,
        };
        if amount > 0 {
            builder.add_edge(n, node, weight).unwrap();
        } else if amount < 0 {
            builder.add_edge(node, n + 1, weight).unwrap();
        }
    }
    let (value, mut residual) = shortest_augmentations(&builder.build(), n, n + 1, required)?;
    if value != required {
        return Ok(None);
    }
    residual.truncate(graph.edge_count());
    Ok(Some(finish(graph, residual)?))
}

fn reject_negative_costs(graph: &Graph<ArcCost>) -> Result<(), CostFlowError> {
    for (edge, arc) in graph.edges().iter().enumerate() {
        if arc.weight.capacity != 0 && arc.weight.unit_cost < 0 {
            return Err(CostFlowError::NegativeArcCost { edge });
        }
    }
    Ok(())
}

fn shortest_augmentations(
    graph: &Graph<ArcCost>,
    source: NodeId,
    sink: NodeId,
    goal: u128,
) -> Result<(u128, Vec<[u64; 2]>), CostFlowError> {
    let n = graph.node_count();
    let mut residual: Vec<_> = graph
        .edges()
        .iter()
        .map(|e| [e.weight.capacity, 0])
        .collect();
    let mut predecessor = vec![None; n];
    let mut prices = vec![0_i128; n];
    let mut distance = vec![u128::MAX; n];
    let mut heap = RadixHeap::new(n);
    let mut value = 0;
    let mut path = Vec::new();
    while value < goal {
        distance.fill(u128::MAX);
        distance[source] = 0;
        heap.reset();
        heap.decrease(source, u128::MAX, 0);
        while let Some((cost, node)) = heap.pop(&distance) {
            for arc in residual_neighbors(graph, node) {
                if residual[arc.edge][arc.reversed] == 0 {
                    continue;
                }
                let weight = unit_cost(graph, arc) + prices[node] - prices[arc.to];
                debug_assert!(weight >= 0);
                let candidate = cost
                    .checked_add(weight as u128)
                    .ok_or(CostFlowError::ArithmeticOverflow)?;
                if candidate < distance[arc.to] {
                    heap.decrease(arc.to, distance[arc.to], candidate);
                    distance[arc.to] = candidate;
                    predecessor[arc.to] = Some((node, arc));
                }
            }
        }
        if distance[sink] == u128::MAX {
            break;
        }
        for (price, &cost) in prices.iter_mut().zip(&distance) {
            if cost != u128::MAX {
                let cost = i128::try_from(cost).map_err(|_| CostFlowError::ArithmeticOverflow)?;
                *price = price
                    .checked_add(cost)
                    .ok_or(CostFlowError::ArithmeticOverflow)?;
            }
        }
        let mut node = sink;
        path.clear();
        let mut amount = goal - value;
        while node != source {
            let (parent, arc) = predecessor[node].unwrap();
            amount = amount.min(u128::from(residual[arc.edge][arc.reversed]));
            path.push(arc);
            node = parent;
        }
        augment(&mut residual, &path, amount as u64);
        value += amount;
    }
    Ok((value, residual))
}

/// Improves a feasible fixed-supply flow by canceling negative residual
/// cycles. Supports every i64 cost and finite u64 capacity. The number of
/// iterations can be pseudopolynomial in capacities and costs.
pub fn minimum_cost_flow_cycle_canceling(
    graph: &Graph<ArcCost>,
    supply: &[i64],
) -> Result<Option<CostFlow>, CostFlowError> {
    let Some(amounts) = feasible_amounts(graph, supply)? else {
        return Ok(None);
    };
    let mut residual: Vec<_> = graph
        .edges()
        .iter()
        .zip(amounts)
        .map(|(e, amount)| [e.weight.capacity - amount, amount])
        .collect();

    let n = graph.node_count();
    let mut distance = vec![0_i128; n];
    let mut predecessor = vec![None; n];
    let mut snapshot = vec![0; n];
    let mut cycle = Vec::new();
    loop {
        distance.fill(0);
        let mut updated = None;
        for _ in 0..n {
            snapshot.copy_from_slice(&distance);
            updated = None;
            for (id, e) in graph.edges().iter().enumerate() {
                for (reverse, from, to) in [(0, e.from, e.to), (1, e.to, e.from)] {
                    if residual[id][reverse] == 0 {
                        continue;
                    }
                    let arc = Residual {
                        edge: id,
                        to,
                        reversed: reverse,
                    };
                    let next = snapshot[from] + unit_cost(graph, arc);
                    if next < distance[to] {
                        updated = Some(to);
                        distance[to] = next;
                        predecessor[to] = Some((from, arc));
                    }
                }
            }
            if updated.is_none() {
                break;
            }
        }
        let Some(mut node) = updated else {
            break;
        };
        for _ in 0..n {
            node = predecessor[node].unwrap().0;
        }
        let root = node;
        cycle.clear();
        loop {
            let (parent, arc) = predecessor[node].unwrap();
            cycle.push(arc);
            node = parent;
            if node == root {
                break;
            }
        }
        debug_assert!(cycle.iter().map(|&arc| unit_cost(graph, arc)).sum::<i128>() < 0);
        let amount = cycle
            .iter()
            .map(|arc| residual[arc.edge][arc.reversed])
            .min()
            .unwrap();
        augment(&mut residual, &cycle, amount);
    }
    Ok(Some(finish(graph, residual)?))
}

fn unit_cost(graph: &Graph<ArcCost>, arc: Residual) -> i128 {
    let cost = i128::from(graph.edges()[arc.edge].weight.unit_cost);
    if arc.reversed == 0 { cost } else { -cost }
}

fn augment(residual: &mut [[u64; 2]], path: &[Residual], amount: u64) {
    for arc in path {
        residual[arc.edge][arc.reversed] -= amount;
        residual[arc.edge][arc.reversed ^ 1] += amount;
    }
}

fn finish(graph: &Graph<ArcCost>, residual: Vec<[u64; 2]>) -> Result<CostFlow, CostFlowError> {
    let edge_flows: Vec<_> = residual.into_iter().map(|pair| pair[1]).collect();
    Ok(CostFlow {
        total_cost: flow_cost(graph, &edge_flows)?,
        edge_flows,
    })
}
