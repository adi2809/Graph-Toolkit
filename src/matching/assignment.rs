//! Sparse rectangular bipartite assignment via reduced-cost augmenting paths.

use std::collections::VecDeque;

use crate::queue::IndexedHeap;
use crate::{EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

use super::{GeneralMatching, VACANT, WeightedMatching};

/// Colors an undirected graph in O(V + E) time. The lowest vertex of each
/// component gets color `false`. An odd cycle or loop yields an error.
pub fn bipartite_coloring<W>(graph: &Graph<W>) -> Result<Vec<bool>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let mut colors = vec![None; graph.node_count()];
    let mut queue = VecDeque::new();
    for root in 0..graph.node_count() {
        if colors[root].is_some() {
            continue;
        }
        colors[root] = Some(false);
        queue.push_back((root, false));
        while let Some((node, color)) = queue.pop_front() {
            for arc in graph.out(node) {
                match colors[arc.node] {
                    Some(other) if color == other => {
                        return Err(GraphError::InvalidBipartition { edge: arc.edge });
                    }
                    None => {
                        colors[arc.node] = Some(!color);
                        queue.push_back((arc.node, !color));
                    }
                    _ => {}
                }
            }
        }
    }
    Ok(colors.into_iter().map(Option::unwrap).collect())
}

/// Cheapest maximum-cardinality bipartite matching with a supplied left side.
/// Every vertex outside `left` belongs to the right side. Duplicate IDs in
/// `left` are accepted. Edges must cross the partition, including zero-cost edges.
///
/// Uses O(V + E) space and O(k(V + E) log V) time for k matched edges.
pub fn minimum_cost_bipartite_matching(
    graph: &Graph<i64>,
    left: &[NodeId],
) -> Result<WeightedMatching, GraphError> {
    let (assignment, _) = assignment(graph, left)?;
    Ok(assignment)
}

/// Cheapest assignment covering the smaller side, or `None` if infeasible.
pub fn minimum_cost_full_bipartite_matching(
    graph: &Graph<i64>,
    left: &[NodeId],
) -> Result<Option<WeightedMatching>, GraphError> {
    let (assignment, full_size) = assignment(graph, left)?;
    Ok((assignment.matching.len() == full_size).then_some(assignment))
}

fn assignment(
    graph: &Graph<i64>,
    left: &[NodeId],
) -> Result<(WeightedMatching, usize), GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    let mut on_left = vec![false; n];
    for &node in left {
        graph.validate_node(node)?;
        on_left[node] = true;
    }
    let mut left = Vec::new();
    for (node, &included) in on_left.iter().enumerate() {
        if included {
            left.push(node);
        }
    }
    let full_size = left.len().min(n - left.len());

    // An implicit source reaches every left vertex at cost zero in the initial
    // residual DAG. The potentials are its shortest distances.
    let mut potential = vec![i128::MAX; n + 1];
    for (id, edge) in graph.edges().iter().enumerate() {
        if on_left[edge.from] == on_left[edge.to] {
            return Err(GraphError::InvalidBipartition { edge: id });
        }
        let right = if on_left[edge.from] {
            edge.to
        } else {
            edge.from
        };
        potential[right] = potential[right].min(i128::from(edge.weight));
    }
    potential[n] = potential[..n].iter().copied().min().unwrap_or(0);
    for value in &mut potential {
        if *value == i128::MAX {
            *value = 0;
        }
    }

    let matching = GeneralMatching {
        mates: vec![VACANT; n],
        len: 0,
    };
    let result = Assignment {
        graph,
        on_left,
        left,
        potential,
        matching,
        distance: vec![i128::MAX; n + 1],
        predecessor: vec![NONE; n + 1],
        queue: IndexedHeap::new(n + 1),
    }
    .solve(full_size);
    Ok((result, full_size))
}

struct Assignment<'a> {
    graph: &'a Graph<i64>,
    left: Vec<NodeId>,
    on_left: Vec<bool>,
    potential: Vec<i128>,
    matching: GeneralMatching,
    distance: Vec<i128>,
    predecessor: Vec<EdgeId>,
    queue: IndexedHeap<i128>,
}

impl Assignment<'_> {
    fn solve(mut self, full_size: usize) -> WeightedMatching {
        while self.matching.len < full_size && self.augment() {}
        let total_weight = self
            .matching
            .edges()
            .map(|id| i128::from(self.graph.edges()[id].weight))
            .sum();
        WeightedMatching {
            matching: self.matching,
            total_weight,
        }
    }

    fn augment(&mut self) -> bool {
        let sink = self.graph.node_count();
        self.distance.fill(i128::MAX);
        self.queue.clear();
        for &node in &self.left {
            if self.matching.mates[node] == VACANT && !self.graph.out(node).is_empty() {
                self.distance[node] = 0;
                self.queue.decrease(node, 0);
            }
        }
        while let Some((distance, node)) = self.queue.pop() {
            if node == sink {
                // Truncation preserves all reduced-cost inequalities after early exit.
                for (potential, &label) in self.potential.iter_mut().zip(&self.distance) {
                    *potential += label.min(distance);
                }
                self.flip(self.predecessor[sink]);
                return true;
            }
            let matched = self.matching.mates[node];
            if self.on_left[node] {
                for arc in self.graph.out(node) {
                    if arc.edge != matched.edge {
                        let cost = i128::from(self.graph.edges()[arc.edge].weight);
                        self.relax(node, arc.node, cost, arc.edge);
                    }
                }
            } else if matched != VACANT {
                let cost = -i128::from(self.graph.edges()[matched.edge].weight);
                self.relax(node, matched.node, cost, matched.edge);
            } else {
                self.relax(node, sink, 0, node);
            }
        }
        false
    }

    fn relax(&mut self, a: NodeId, b: NodeId, cost: i128, predecessor: usize) {
        let reduced = cost + self.potential[a] - self.potential[b];
        debug_assert!(reduced >= 0);
        let candidate = self.distance[a] + reduced;
        if candidate < self.distance[b] {
            self.distance[b] = candidate;
            self.predecessor[b] = predecessor;
            self.queue.decrease(b, candidate);
        }
    }

    fn flip(&mut self, mut right: NodeId) {
        loop {
            let id = self.predecessor[right];
            let edge = &self.graph.edges()[id];
            let left = if edge.from == right {
                edge.to
            } else {
                edge.from
            };
            let previous = self.matching.mates[left];
            self.matching.pair(left, right, id);
            if previous == VACANT {
                break;
            }
            right = previous.node;
        }
        self.matching.len += 1;
    }
}
