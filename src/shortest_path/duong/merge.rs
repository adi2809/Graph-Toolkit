use super::cover::Child;
use super::{Shifted, StepGraph};
use crate::connectivity::strongly_connected_components;
use crate::queue::RadixHeap;
use crate::{EdgeId, GraphBuilder, GraphKind, NONE};

/// Exact nonnegative closure alternating with synchronous negative relaxations.
/// An improvement beyond the given budget returns a concrete comparison walk.
pub(super) fn hybrid(
    graph: &StepGraph,
    initial: Vec<i128>,
    budget: usize,
) -> Result<Vec<i128>, Vec<EdgeId>> {
    let n = graph.node_count();
    // Flip the sign bit to preserve signed priority order in the radix heap.
    const SIGN: u128 = 1 << 127;
    let mut distance: Vec<_> = initial
        .into_iter()
        .map(|value| value as u128 ^ SIGN)
        .collect();
    let mut heap = RadixHeap::new(n);
    let mut queued = vec![true; n];
    let mut history = Vec::new();
    let mut trace = vec![NONE; n];
    let mut parent = vec![(NONE, NONE); n];
    let mut seeds = vec![u128::MAX; n];
    let negative: Vec<_> = graph.edges().iter().filter(|e| e.weight.cost < 0).collect();
    for (node, &cost) in distance.iter().enumerate() {
        heap.decrease(node, u128::MAX, cost);
    }

    for round in 0..=budget {
        while let Some((cost, node)) = heap.pop(&distance) {
            queued[node] = false;
            let (previous, edge) = parent[node];
            if edge != NONE {
                trace[node] = history.len();
                history.push((previous, edge));
            } else {
                trace[node] = previous;
            }
            for arc in graph.out(node) {
                let edge = graph.edges()[arc.edge].weight;
                if edge.cost >= 0 {
                    let next = cost + edge.cost as u128;
                    let old = distance[arc.node];
                    if next < old {
                        let occupied = std::mem::replace(&mut queued[arc.node], true);
                        let key = if occupied { old } else { u128::MAX };
                        heap.decrease(arc.node, key, next);
                        distance[arc.node] = next;
                        parent[arc.node] = (trace[node], edge.original);
                    }
                }
            }
        }

        for edge in &negative {
            let next = distance[edge.from] - edge.weight.cost.unsigned_abs();
            if next >= distance[edge.to] || next >= seeds[edge.to] {
                continue;
            }
            let prefix = trace[edge.from];
            if round == budget {
                let mut walk = Vec::new();
                if edge.weight.original != NONE {
                    walk.push(edge.weight.original);
                }
                let mut prefix = prefix;
                while prefix != NONE {
                    let (previous, edge) = history[prefix];
                    walk.push(edge);
                    prefix = previous;
                }
                walk.reverse();
                return Err(walk);
            }
            seeds[edge.to] = next;
            parent[edge.to] = (prefix, edge.weight.original);
        }

        heap.reset();
        let mut changed = false;
        for (node, seed) in seeds.iter_mut().enumerate() {
            if *seed != u128::MAX {
                distance[node] = std::mem::replace(seed, u128::MAX);
                heap.decrease(node, u128::MAX, distance[node]);
                queued[node] = true;
                changed = true;
            }
        }
        if !changed {
            return Ok(distance
                .into_iter()
                .map(|value| (value ^ SIGN) as i128)
                .collect());
        }
    }
    unreachable!()
}

struct CopyNode {
    child: usize,
    vertex: usize,
    price: i128,
}

/// Shift child prices and join overlapping copies using shared prefix and
/// suffix chains. Only transitions to a later child can have negative cost.
pub(super) fn ordered(
    graph: &StepGraph,
    children: Vec<Child>,
    prices: Vec<Vec<i128>>,
    ascent_limit: usize,
) -> Result<Vec<i128>, Vec<EdgeId>> {
    let n = graph.node_count();
    let mut earlier = vec![None::<i128>; n];
    let mut memberships = (0..n).map(|_| Vec::new()).collect::<Vec<_>>();
    let mut copies = Vec::new();
    for (i, (child, prices)) in children.into_iter().zip(prices).enumerate() {
        let mut shift = 0;
        for (&node, &price) in child.nodes.iter().zip(&prices) {
            for arc in graph.out(node) {
                if let Some(value) = earlier[arc.node] {
                    shift = shift.max(value - graph.edges()[arc.edge].weight.cost - price);
                }
            }
        }
        for (node, price) in child.nodes.into_iter().zip(prices) {
            let price = price + shift;
            earlier[node] = Some(earlier[node].map_or(price, |old| old.max(price)));
            memberships[node].push(copies.len());
            copies.push(CopyNode {
                child: i,
                vertex: node,
                price,
            });
        }
    }

    let mut builder = GraphBuilder::new(copies.len() * 3, GraphKind::Directed);
    let mut price = vec![0; copies.len() * 3];
    let mut add = |from, to, cost, original| {
        builder
            .add_edge(from, to, Shifted { cost, original })
            .unwrap();
    };
    for ids in &memberships {
        let mut largest = i128::MIN;
        for (position, &id) in ids.iter().enumerate() {
            let vertex = 3 * id;
            price[vertex] = copies[id].price;
            largest = largest.max(price[vertex]);
            price[vertex + 1] = largest;
            add(vertex + 1, vertex, largest - price[vertex], NONE);
            if position > 0 {
                let previous = 3 * ids[position - 1] + 1;
                add(vertex + 1, previous, largest - price[previous], NONE);
            }
        }
        largest = i128::MIN;
        for (position, &id) in ids.iter().enumerate().rev() {
            let vertex = 3 * id;
            largest = largest.max(price[vertex]);
            price[vertex + 2] = largest;
            add(vertex + 2, vertex, largest - price[vertex], NONE);
            if let Some(&next) = ids.get(position + 1) {
                let next = 3 * next + 2;
                add(vertex + 2, next, largest - price[next], NONE);
            }
        }
    }
    // Merge the sorted membership lists for each logical arc. The two
    // connector indices advance monotonically over its source copies.
    for edge in graph.edges() {
        let to = &memberships[edge.to];
        let mut split = 0;
        for &id in &memberships[edge.from] {
            let copy = &copies[id];
            debug_assert_eq!(copy.vertex, edge.from);
            while split < to.len() && copies[to[split]].child <= copy.child {
                split += 1;
            }
            let from = 3 * id;
            let cost = edge.weight.cost + price[from];
            if split > 0 {
                let target = 3 * to[split - 1] + 1;
                debug_assert!(cost >= price[target]);
                add(from, target, cost - price[target], edge.weight.original);
            }
            if let Some(&target) = to.get(split) {
                let target = 3 * target + 2;
                add(from, target, cost - price[target], edge.weight.original);
            }
        }
    }
    let auxiliary = builder.build();
    let initial = price.iter().map(|&value| -value).collect();
    let distances = hybrid(
        &auxiliary,
        initial,
        (ascent_limit + 3).min(auxiliary.node_count()),
    )?;
    Ok(memberships
        .into_iter()
        .map(|ids| {
            let first = 3 * ids[0];
            distances[first] + price[first]
        })
        .collect())
}

/// Seed the relaxation with exact distances using only nonpositive arcs.
/// The seeds correspond to real walks and can skip many hybrid rounds.
pub(super) fn nonpositive_seed(graph: &StepGraph) -> Vec<i128> {
    let n = graph.node_count();
    let mut builder = GraphBuilder::new(n, GraphKind::Directed);
    for edge in graph.edges().iter().filter(|edge| edge.weight.cost <= 0) {
        builder
            .add_edge(edge.from, edge.to, edge.weight.cost)
            .unwrap();
    }
    let nonpositive = builder.build();
    let components = strongly_connected_components(&nonpositive);
    let mut groups = vec![Vec::new(); components.count];
    for edge in nonpositive.edges() {
        let from = components.component_of[edge.from];
        let to = components.component_of[edge.to];
        if from != to {
            groups[from].push((to, edge.weight));
        } else if edge.weight < 0 {
            return vec![0; n];
        }
    }
    let mut prices = vec![0; components.count];
    for (from, outgoing) in groups.into_iter().enumerate() {
        for (to, cost) in outgoing {
            prices[to] = prices[to].min(prices[from] + cost);
        }
    }
    components
        .component_of
        .into_iter()
        .map(|id| prices[id])
        .collect()
}
