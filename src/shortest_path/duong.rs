//! Goldberg scaling with Duong's local padded decomposition and ordered merge.

mod cover;
mod merge;

use super::signed::{
    NegativeCycle, SignedPathError, SignedShortestPaths, UNREACHABLE, bellman_ford,
    relax_by_component,
};
use super::{INF, NO_LINK, dijkstra_core};
use crate::queue::RadixHeap;
use crate::{EdgeId, Graph, GraphBuilder, GraphKind, NONE, NodeId};

#[derive(Clone, Copy)]
struct Shifted {
    cost: i128,
    original: EdgeId,
}

type StepGraph = Graph<Shifted>;

/// Exact deterministic signed single-source paths using Duong's local padded
/// decomposition (Lemma 4.9), ordered merge and scaling with cycle witnesses.
/// A bounded SCC-aware Bellman–Ford pass handles easy instances upfront.
/// Only negative cycles reachable from the source are reported.
///
/// The size-sensitive peeling regime of the paper starts at L = 256, beyond
/// this crate's index widths. This implementation uses a radix heap and does
/// not claim the paper's word-RAM bound.
pub fn duong_shortest_paths(
    graph: &Graph<i64>,
    source: NodeId,
) -> Result<SignedShortestPaths, SignedPathError> {
    graph.validate_node(source)?;
    let n = graph.node_count();
    let signed = graph.edges().iter().any(|edge| edge.weight < 0);
    if signed {
        let mut distance = vec![UNREACHABLE; n];
        distance[source] = 0;
        if let Some(result) = relax_by_component(graph, distance, 48) {
            return Ok(result?);
        }
    }
    scaled_shortest_paths(graph, source, signed)
}

fn scaled_shortest_paths(
    graph: &Graph<i64>,
    source: NodeId,
    signed: bool,
) -> Result<SignedShortestPaths, SignedPathError> {
    let n = graph.node_count();
    let mut predecessor = vec![NO_LINK; n];
    let mut distance = vec![INF; n];
    let mut queue = RadixHeap::new(n);
    let mut potential = vec![0_i128; n];

    if signed {
        let mut nodes = vec![source];
        let mut local = vec![NONE; n];
        local[source] = 0;
        let mut head = 0;
        while let Some(&node) = nodes.get(head) {
            head += 1;
            for arc in graph.out(node) {
                if local[arc.node] == NONE {
                    local[arc.node] = nodes.len();
                    nodes.push(arc.node);
                }
            }
        }

        let mut edges = Vec::new();
        let mut low = 0_i128;
        let mut magnitude = 0_i128;
        let factor = nodes.len() as i128 + 1;
        for (from, &node) in nodes.iter().enumerate() {
            for arc in graph.out(node) {
                let weight = i128::from(graph.edges()[arc.edge].weight);
                if weight < 0 && (arc.node == node || graph.kind() == GraphKind::Undirected) {
                    let (nodes, edges, cost) = if arc.node == node {
                        (vec![node, node], vec![arc.edge], weight)
                    } else {
                        (vec![node, arc.node, node], vec![arc.edge; 2], weight * 2)
                    };
                    return Err(NegativeCycle { nodes, edges, cost }.into());
                }
                let cost = weight * factor;
                low = low.min(cost);
                magnitude = magnitude.max(cost.abs());
                edges.push((
                    from,
                    local[arc.node],
                    Shifted {
                        cost,
                        original: arc.edge,
                    },
                ));
            }
        }
        if low < 0 {
            // Keep exact signed arithmetic for the wider intermediate distances.
            // On exceptionally large numerical instances use the existing exact solver.
            let fits = magnitude
                .checked_mul(factor)
                .and_then(|value| value.checked_mul(edges.len() as i128 + 1))
                .is_some_and(|value| value <= i128::MAX / 1024);
            if !fits {
                return bellman_ford(graph, source);
            }
            let mut prices = vec![0; nodes.len()];
            let total_mass = 2 * edges.len() + nodes.len();
            let log_mass = (total_mass - 1).ilog2() as usize + 1;
            let mut mass = vec![1; nodes.len()];
            for &(from, to, _) in &edges {
                mass[from] += 1;
                mass[to] += 1;
            }
            while low < -1 {
                let shift = (-low).div_euclid(2) + (-low).rem_euclid(2);
                let mut builder =
                    GraphBuilder::with_edge_capacity(nodes.len(), GraphKind::Directed, edges.len());
                for &(from, to, mut edge) in &edges {
                    edge.cost += prices[from] - prices[to] + shift;
                    builder.add_edge(from, to, edge).unwrap();
                }
                let step = builder.build();
                // A simple path cannot weigh less than -(|V|-1)*2*shift.
                let diameter = (2 * nodes.len() as u128 * shift as u128).next_power_of_two();
                let solver = ScalingStep {
                    root: &step,
                    graph,
                    nodes: &nodes,
                    shift,
                    diameter,
                    log_mass,
                };
                let change = solver.solve(&step, &mass, diameter)?;
                for (price, delta) in prices.iter_mut().zip(change) {
                    *price += delta;
                }
                low = edges
                    .iter()
                    .map(|&(u, v, e)| e.cost + prices[u] - prices[v])
                    .min()
                    .unwrap();
                debug_assert!(low >= -shift);
            }

            // Every scaled reduced arc weighs at least -1. Choose a common
            // offset so no price has residue factor-1 before rounding down.
            let mut used = vec![false; nodes.len() + 1];
            for &price in &prices {
                used[price.rem_euclid(factor) as usize] = true;
            }
            let gap = used.iter().position(|&taken| !taken).unwrap() as i128;
            let offset = factor - 1 - gap;
            for (&node, price) in nodes.iter().zip(prices) {
                potential[node] = (price + offset).div_euclid(factor);
            }
        }
    }

    // The original graph can include unreachable arcs with negative weights.
    // The Dijkstra kernel invokes this callback only for reachable tails.
    dijkstra_core(
        graph,
        source,
        false,
        |arc| {
            let edge = graph.edges()[arc.edge];
            let reduced = i128::from(edge.weight) + potential[edge.from] - potential[edge.to];
            debug_assert!(reduced >= 0);
            reduced as u128
        },
        &mut distance,
        &mut queue,
        |node, parent, edge| predecessor[node] = (parent, edge),
    );
    let distance = distance
        .into_iter()
        .zip(&potential)
        .map(|(cost, price)| {
            if cost == INF {
                UNREACHABLE
            } else {
                cost as i128 + price - potential[source]
            }
        })
        .collect();
    Ok(SignedShortestPaths {
        distance,
        predecessor,
    })
}

struct ScalingStep<'a> {
    root: &'a StepGraph,
    graph: &'a Graph<i64>,
    nodes: &'a [NodeId],
    shift: i128,
    diameter: u128,
    log_mass: usize,
}

impl ScalingStep<'_> {
    fn solve(
        &self,
        graph: &StepGraph,
        mass: &[usize],
        diameter: u128,
    ) -> Result<Vec<i128>, NegativeCycle> {
        let n = graph.node_count();
        let Some(edge) = graph.edges().iter().find(|e| e.weight.cost < 0) else {
            return Ok(vec![0; n]);
        };
        if diameter < self.shift as u128 {
            return Err(self.witness(vec![edge.weight.original], diameter));
        }
        // A bounded hybrid pass often solves an entire child without partitioning.
        // Small children need at most n-1 negative arcs on any simple path.
        let initial = if n > 4 && diameter == self.diameter && std::ptr::eq(graph, self.root) {
            merge::nonpositive_seed(graph)
        } else {
            vec![0; n]
        };
        match merge::hybrid(graph, initial, n.min(4)) {
            Ok(prices) => return Ok(prices),
            Err(walk) if n <= 4 => return Err(self.witness(walk, self.diameter)),
            Err(_) => {}
        }

        let (children, ascent_limit) = cover::decompose(graph, mass, diameter, self.log_mass);
        let mut local = vec![NONE; n];
        let mut prices = Vec::with_capacity(children.len());
        for child in &children {
            for (id, &node) in child.nodes.iter().enumerate() {
                local[node] = id;
            }
            let mut builder = GraphBuilder::new(child.nodes.len(), GraphKind::Directed);
            for &node in &child.nodes {
                for arc in graph.out(node) {
                    if local[arc.node] != NONE {
                        builder
                            .add_edge(local[node], local[arc.node], graph.edges()[arc.edge].weight)
                            .unwrap();
                    }
                }
            }
            let child_mass: Vec<_> = child.nodes.iter().map(|&u| mass[u]).collect();
            prices.push(self.solve(&builder.build(), &child_mass, child.diameter)?);
            for &node in &child.nodes {
                local[node] = NONE;
            }
        }
        merge::ordered(graph, children, prices, ascent_limit)
            .map_err(|walk| self.witness(walk, diameter))
    }

    fn witness(&self, mut edges: Vec<EdgeId>, diameter: u128) -> NegativeCycle {
        // Remove a negative simple cycle from the walk, closing it in the
        // root clipped graph first if its weak diameter was certified.
        let start = self.graph.edges()[edges[0]].from;
        let end = self.graph.edges()[*edges.last().unwrap()].to;
        if start != end && diameter < self.diameter {
            let (source, target) = self
                .nodes
                .iter()
                .enumerate()
                .fold((NONE, NONE), |(s, t), (i, &v)| {
                    (if v == end { i } else { s }, if v == start { i } else { t })
                });
            let n = self.root.node_count();
            let mut distance = vec![INF; n];
            let mut predecessor = vec![NO_LINK; n];
            dijkstra_core(
                self.root,
                source,
                false,
                |arc| self.root.edges()[arc.edge].weight.cost.max(0) as u128,
                &mut distance,
                &mut RadixHeap::new(n),
                |node, parent, edge| predecessor[node] = (parent, edge),
            );
            debug_assert!(distance[target] <= diameter);
            let mut current = target;
            let end_of_walk = edges.len();
            while current != source {
                let (parent, edge) = predecessor[current];
                edges.push(self.root.edges()[edge].weight.original);
                current = parent;
            }
            edges[end_of_walk..].reverse();
        }
        let mut position = vec![NONE; self.graph.node_count()];
        let mut stack = vec![start];
        let mut simple = Vec::new();
        let mut prefix = vec![0_i128];
        position[start] = 0;
        for id in edges {
            let edge = self.graph.edges()[id];
            let total = prefix.last().unwrap() + i128::from(edge.weight);
            simple.push(id);
            let previous = position[edge.to];
            if previous == NONE {
                position[edge.to] = stack.len();
                stack.push(edge.to);
                prefix.push(total);
            } else {
                let cost = total - prefix[previous];
                if cost < 0 {
                    stack.push(edge.to);
                    return NegativeCycle {
                        cost,
                        nodes: stack.split_off(previous),
                        edges: simple.split_off(previous),
                    };
                }
                for node in stack.drain(previous + 1..) {
                    position[node] = NONE;
                }
                simple.truncate(previous);
                prefix.truncate(previous + 1);
            }
        }
        unreachable!("Duong certificate must contain a negative cycle");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scaling_without_startup_pass_matches_bellman_ford() {
        let mut random = 4390587_u64;
        let mut sample = |bound: usize| {
            random ^= random << 13;
            random ^= random >> 7;
            random ^= random << 17;
            random as usize % bound
        };
        for run in 0..220 {
            let n = 7 + sample(33);
            let mut builder = GraphBuilder::new(n, GraphKind::Directed);
            let price: Vec<_> = (0..n).map(|_| sample(281) as i64 - 140).collect();
            for from in 0..n {
                let to = (from + 1) % n;
                let cost = if run % 4 == 0 { -1 } else { sample(4) as i64 };
                builder
                    .add_edge(from, to, cost + price[to] - price[from])
                    .unwrap();
                for _ in 0..run % 5 {
                    let target = sample(n);
                    let cost = sample(89) as i64;
                    builder
                        .add_edge(from, target, cost + price[target] - price[from])
                        .unwrap();
                }
            }
            let graph = builder.build();
            match (
                bellman_ford(&graph, 0),
                scaled_shortest_paths(&graph, 0, true),
            ) {
                (Ok(expected), Ok(actual)) => {
                    assert!(expected.distances().eq(actual.distances()), "case {run}")
                }
                (
                    Err(SignedPathError::NegativeCycle(_)),
                    Err(SignedPathError::NegativeCycle(cycle)),
                ) => {
                    assert_eq!(cycle.nodes.first(), cycle.nodes.last());
                    assert_eq!(cycle.nodes.len(), cycle.edges.len() + 1);
                    let mut cost = 0_i128;
                    for (&id, pair) in cycle.edges.iter().zip(cycle.nodes.windows(2)) {
                        let edge = graph.edges()[id];
                        assert_eq!([edge.from, edge.to], pair);
                        cost += i128::from(edge.weight);
                    }
                    assert_eq!(cycle.cost, cost);
                    assert!(cost < 0);
                }
                pair => panic!("case {run}: {pair:?}"),
            }
        }
    }
}
