//! Gabow's phased maximum-cardinality matching. The canonical-path formulation
//! of the level search follows Mehlhorn and Nobahari (2026), arXiv:2603.22909.

use crate::{Arc, EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

use super::tree_union::TreeUnion;
use super::{GeneralMatching, VACANT, greedy_matching};

/// Exact general-graph matching in O(V + E√V) time and O(V + E) space.
/// Augments a maximal set of vertex-disjoint shortest augmenting paths per phase.
/// Parallel edges retain their original IDs. Edge weights and self-loops are ignored.
pub fn maximum_general_matching_gabow<W>(graph: &Graph<W>) -> Result<GeneralMatching, GraphError> {
    graph.require(GraphKind::Undirected)?;
    if graph.node_count() <= graph.edge_count().saturating_mul(2) {
        return Ok(solve(graph));
    }

    // Remove isolated vertices from the per-phase state when V dominates E.
    let mut originals = Vec::new();
    let mut compact = vec![NONE; graph.node_count()];
    for (node, id) in compact.iter_mut().enumerate() {
        if graph.out(node).iter().any(|arc| arc.node != node) {
            *id = originals.len();
            originals.push(node);
        }
    }
    let mut builder = Graph::builder(originals.len(), GraphKind::Undirected);
    for (id, edge) in graph.edges().iter().enumerate() {
        if edge.from != edge.to {
            builder
                .add_edge(compact[edge.from], compact[edge.to], id)
                .unwrap();
        }
    }
    let interior = builder.build();
    let matching = solve(&interior);
    let mut result = GeneralMatching {
        mates: vec![VACANT; graph.node_count()],
        len: matching.len,
    };
    for id in matching.edges() {
        let edge = interior.edges()[id];
        result.pair(originals[edge.from], originals[edge.to], edge.weight);
    }
    Ok(result)
}

fn solve<W>(graph: &Graph<W>) -> GeneralMatching {
    let mut matching = greedy_matching(graph);
    let upper_bound = (graph.node_count() / 2).min(2 * matching.len);
    if matching.len < upper_bound {
        augment_length_three(graph, &mut matching, upper_bound);
    }
    while matching.len < upper_bound {
        let mut search = LevelSearch::new(graph, &matching);
        let Some(delta) = search.shortest_augmenting_level() else {
            break;
        };
        let representative = search.representatives(delta);
        let auxiliary = search.auxiliary_graph(delta, &representative);
        let mates: Vec<_> = matching.mates.iter().map(|arc| arc.node).collect();
        let mate = auxiliary_mates(&auxiliary, graph, &mates);
        let (mut replacement, count) = BlockingSearch::new(&auxiliary, &mate).paths();
        assert!(
            count > 0,
            "the shortest augmenting level must contain an augmenting path"
        );

        // Each original blossom is visited by at most one path in this phase.
        let mut stack = Vec::new();
        // Only the non-matching edges in H need to be expanded here. The new
        // edges produced inside old blossoms already belong to the original graph.
        let boundary_edges = replacement.len();
        for cursor in 0..boundary_edges {
            let edge = &graph.edges()[replacement[cursor]];
            for node in [edge.from, edge.to] {
                search.path.collect(
                    node,
                    representative[node],
                    &mates,
                    &mut stack,
                    &mut replacement,
                );
            }
        }
        for id in replacement {
            let edge = &graph.edges()[id];
            matching.pair(edge.from, edge.to, id);
        }
        matching.len += count;
    }
    matching
}

/// The greedy matching is maximal, so its shortest possible augmentation
/// uses three edges. Process each original matched edge at most once.
fn augment_length_three<W>(graph: &Graph<W>, matching: &mut GeneralMatching, upper_bound: usize) {
    let initial_pairs: Vec<_> = matching
        .mates
        .iter()
        .enumerate()
        .filter_map(|(a, mate)| (mate.node != NONE && a < mate.node).then_some((a, mate.node)))
        .collect();
    for (a, b) in initial_pairs {
        // Augmentations only shrink the set of free vertices. A new matched
        // edge touches a formerly free vertex, which has no free neighbors.
        let Some(mut left) = free_neighbor(graph, matching, a, NONE) else {
            continue;
        };
        let Some(mut right) = free_neighbor(graph, matching, b, NONE) else {
            continue;
        };
        if left.node == right.node {
            if let Some(other) = free_neighbor(graph, matching, a, right.node) {
                left = other;
            } else if let Some(other) = free_neighbor(graph, matching, b, left.node) {
                right = other;
            } else {
                continue;
            }
        }
        matching.pair(a, left.node, left.edge);
        matching.pair(b, right.node, right.edge);
        matching.len += 1;
        if matching.len == upper_bound {
            break;
        }
    }
}

fn free_neighbor<W>(
    graph: &Graph<W>,
    matching: &GeneralMatching,
    vertex: NodeId,
    exclude: NodeId,
) -> Option<Arc> {
    graph
        .out(vertex)
        .iter()
        .copied()
        .find(|arc| arc.node != exclude && matching.mates[arc.node] == VACANT)
}

fn auxiliary_mates<W>(
    auxiliary: &Graph<EdgeId>,
    graph: &Graph<W>,
    mates: &[NodeId],
) -> Vec<NodeId> {
    let mut result = vec![NONE; graph.node_count()];
    for edge in auxiliary.edges() {
        let original = &graph.edges()[edge.weight];
        if mates[original.from] == original.to {
            result[edge.from] = edge.to;
            result[edge.to] = edge.from;
        }
    }
    result
}

#[derive(Clone, Copy)]
struct Bridge {
    // The first endpoint lies on the arm on which this odd vertex was found.
    from: NodeId,
    to: NodeId,
    edge: EdgeId,
}

struct CanonicalPaths {
    // Stores the non-matching tree edge entering each vertex born odd.
    incoming: Vec<Arc>,
    bridge: Vec<Option<Bridge>>,
}

impl CanonicalPaths {
    fn new(n: usize) -> Self {
        Self {
            incoming: vec![VACANT; n],
            bridge: vec![None; n],
        }
    }

    /// Collects the replacement non-matching edges along an even alternating path.
    /// Order is immaterial because the lifted paths are vertex-disjoint.
    fn collect(
        &self,
        start: NodeId,
        end: NodeId,
        mate: &[NodeId],
        stack: &mut Vec<(NodeId, NodeId)>,
        edges: &mut Vec<EdgeId>,
    ) {
        stack.push((start, end));
        while let Some((mut node, stop)) = stack.pop() {
            while node != stop {
                let partner = mate[node];
                if let Some(bridge) = self.bridge[node] {
                    edges.push(bridge.edge);
                    stack.push((bridge.from, partner));
                    node = bridge.to;
                } else {
                    let parent = self.incoming[partner];
                    edges.push(parent.edge);
                    node = parent.node;
                }
            }
        }
    }
}

struct LevelSearch<'a, W> {
    graph: &'a Graph<W>,
    matching: &'a GeneralMatching,
    sets: TreeUnion,
    /// Creation order is a topological order of the original alternating forest.
    discovered: Vec<NodeId>,
    root: Vec<NodeId>,
    even: Vec<usize>,
    odd: Vec<usize>,
    linked_in: Vec<usize>,
    path: CanonicalPaths,
    growth: Vec<Vec<NodeId>>,
    bridges: Vec<Vec<Bridge>>,
    seen: Vec<usize>,
    epoch: usize,
}

impl<'a, W> LevelSearch<'a, W> {
    fn new(graph: &'a Graph<W>, matching: &'a GeneralMatching) -> Self {
        let n = graph.node_count();
        let buckets = n / 2 + 1;
        Self {
            graph,
            matching,
            sets: TreeUnion::new(n),
            discovered: Vec::new(),
            root: vec![NONE; n],
            even: vec![NONE; n],
            odd: vec![NONE; n],
            linked_in: vec![NONE; n],
            path: CanonicalPaths::new(n),
            growth: vec![Vec::new(); buckets],
            bridges: vec![Vec::new(); buckets],
            seen: vec![0; n],
            epoch: 0,
        }
    }

    fn shortest_augmenting_level(&mut self) -> Option<usize> {
        for (node, mate) in self.matching.mates.iter().enumerate() {
            if *mate == VACANT {
                self.root[node] = node;
                self.discovered.push(node);
                self.sets.add_root(node);
                self.make_even(node, 0, 0);
            }
        }

        for delta in 1..self.growth.len() {
            while let Some(node) = self.growth[delta].pop() {
                self.grow_from(node, delta);
            }
            while let Some(bridge) = self.bridges[delta].pop() {
                let a = self.sets.find(bridge.from);
                let b = self.sets.find(bridge.to);
                if a == b {
                    continue;
                }
                if self.root[a] != self.root[b] {
                    return Some(delta);
                }
                let base = self.common_base(a, b);
                self.shrink_arm(a, base, bridge, delta);
                self.shrink_arm(
                    b,
                    base,
                    Bridge {
                        from: bridge.to,
                        to: bridge.from,
                        edge: bridge.edge,
                    },
                    delta,
                );
            }
        }
        None
    }

    fn make_even(&mut self, node: NodeId, level: usize, delta: usize) {
        self.even[node] = level;
        if let Some(bucket) = self.growth.get_mut(level + 2) {
            bucket.push(node);
        }
        for arc in self.graph.out(node) {
            let next = arc.node;
            if self.even[next] == NONE || node == next || self.matching.mates[node].node == next {
                continue;
            }
            let when = (level + self.even[next]) / 2 + 1;
            debug_assert!(when >= delta);
            if let Some(bucket) = self.bridges.get_mut(when) {
                bucket.push(Bridge {
                    from: node,
                    to: next,
                    edge: arc.edge,
                });
            }
        }
    }

    fn grow_from(&mut self, node: NodeId, delta: usize) {
        for &arc in self.graph.out(node) {
            let odd = arc.node;
            if self.root[odd] != NONE {
                continue;
            }
            let even = self.matching.mates[odd].node;
            self.root[odd] = self.root[node];
            self.root[even] = self.root[node];
            self.odd[odd] = delta - 1;
            self.path.incoming[odd] = Arc {
                node,
                edge: arc.edge,
            };
            self.discovered.extend([odd, even]);
            self.sets.grow(node, odd);
            self.sets.grow(odd, even);
            self.make_even(even, delta, delta);
        }
    }

    fn common_base(&mut self, mut a: NodeId, mut b: NodeId) -> NodeId {
        self.epoch += 1;
        // Lockstep walks charge to the contractions on the two arms, even
        // when one side reaches the common base much earlier than the other.
        loop {
            if a != NONE {
                if self.seen[a] == self.epoch {
                    return a;
                }
                self.seen[a] = self.epoch;
                let odd = self.matching.mates[a].node;
                a = if odd == NONE {
                    NONE
                } else {
                    self.sets.find(self.path.incoming[odd].node)
                };
            }
            std::mem::swap(&mut a, &mut b);
        }
    }

    fn shrink_arm(&mut self, mut vertex: NodeId, base: NodeId, bridge: Bridge, delta: usize) {
        while vertex != base {
            let odd = self.matching.mates[vertex].node;
            self.sets.link(vertex);
            self.sets.link(odd);
            self.linked_in[vertex] = delta;
            self.linked_in[odd] = delta;
            self.path.bridge[odd] = Some(bridge);
            self.make_even(odd, 2 * delta - 1 - self.odd[odd], delta);
            vertex = self.sets.find(self.path.incoming[odd].node);
        }
    }

    fn representatives(&self, delta: usize) -> Vec<NodeId> {
        // Only contractions completed before the breakthrough belong in H.
        let mut rep: Vec<_> = (0..self.graph.node_count()).collect();
        for &node in &self.discovered {
            if self.linked_in[node] < delta {
                let parent = if self.odd[node] == NONE {
                    self.matching.mates[node].node
                } else {
                    self.path.incoming[node].node
                };
                rep[node] = rep[parent];
            }
        }
        rep
    }

    fn level_before_breakthrough(&self, node: NodeId, delta: usize) -> Level {
        if self.odd[node] != NONE && self.linked_in[node] >= delta {
            Level::Odd(self.odd[node])
        } else if self.even[node] != NONE {
            Level::Even(self.even[node])
        } else {
            Level::Unseen
        }
    }

    fn auxiliary_graph(&self, delta: usize, rep: &[NodeId]) -> Graph<EdgeId> {
        let mut builder = Graph::builder(self.graph.node_count(), GraphKind::Undirected);
        for (id, edge) in self.graph.edges().iter().enumerate() {
            let (a, b) = (edge.from, edge.to);
            if rep[a] == rep[b] {
                continue;
            }
            let paired = self.matching.mates[a].node == b;
            let matched = paired && self.matching.mates[a].edge == id;
            if paired && !matched {
                continue;
            }
            let left = self.level_before_breakthrough(a, delta);
            let right = self.level_before_breakthrough(b, delta);
            if admissible(left, right, matched, delta) {
                builder.add_edge(rep[a], rep[b], id).unwrap();
            }
        }
        builder.build()
    }
}

#[derive(Clone, Copy)]
enum Level {
    Even(usize),
    Odd(usize),
    Unseen,
}

fn admissible(a: Level, b: Level, matched: bool, delta: usize) -> bool {
    match (a, b) {
        (Level::Even(a), Level::Even(b)) => !matched && a + b + 2 == 2 * delta,
        (Level::Even(even), Level::Odd(odd)) | (Level::Odd(odd), Level::Even(even)) => {
            if matched {
                even == odd + 1
            } else {
                odd == even + 1
            }
        }
        (Level::Even(even), Level::Unseen) | (Level::Unseen, Level::Even(even)) => {
            !matched && even + 2 == delta
        }
        (Level::Unseen, Level::Unseen) => matched,
        _ => false,
    }
}

struct BlockingSearch<'a> {
    graph: &'a Graph<EdgeId>,
    mate: &'a [NodeId],
    sets: TreeUnion,
    visited: Vec<bool>,
    even_time: Vec<usize>,
    clock: usize,
    path: CanonicalPaths,
    dfs: Vec<(NodeId, usize)>,
    lift: Vec<(NodeId, NodeId)>,
}

impl<'a> BlockingSearch<'a> {
    fn new(graph: &'a Graph<EdgeId>, mate: &'a [NodeId]) -> Self {
        let n = graph.node_count();
        Self {
            graph,
            mate,
            sets: TreeUnion::new(n),
            visited: vec![false; n],
            even_time: vec![0; n],
            clock: 0,
            path: CanonicalPaths::new(n),
            dfs: Vec::new(),
            lift: Vec::new(),
        }
    }

    fn paths(mut self) -> (Vec<EdgeId>, usize) {
        let mut replacement = Vec::new();
        let mut count = 0;
        for root in 0..self.graph.node_count() {
            if self.visited[root] || self.mate[root] != NONE || self.graph.out(root).is_empty() {
                continue;
            }
            self.sets.add_root(root);
            self.push_even(root);
            if let Some((node, edge)) = self.find_path() {
                self.path
                    .collect(node, root, self.mate, &mut self.lift, &mut replacement);
                replacement.push(edge);
                count += 1;
                self.dfs.clear();
            }
        }
        (replacement, count)
    }

    fn push_even(&mut self, node: NodeId) {
        self.clock += 1;
        self.visited[node] = true;
        self.even_time[node] = self.clock;
        self.dfs.push((node, 0));
    }

    fn find_path(&mut self) -> Option<(NodeId, EdgeId)> {
        while let Some((node, cursor)) = self.dfs.last_mut() {
            let node = *node;
            let Some(&arc) = self.graph.out(node).get(*cursor) else {
                self.dfs.pop();
                continue;
            };
            *cursor += 1;
            let next = arc.node;
            if self.mate[node] == next {
                continue;
            }

            let edge = self.graph.edges()[arc.edge].weight;
            if !self.visited[next] {
                self.visited[next] = true;
                self.sets.grow(node, next);
                if self.mate[next] == NONE {
                    return Some((node, edge));
                }
                self.path.incoming[next] = Arc { node, edge };
                self.sets.grow(next, self.mate[next]);
                self.push_even(self.mate[next]);
                continue;
            }
            // An edge to a proper descendant is a forward edge in the DFS
            // forest. Other visited vertices need no further exploration.
            let base = self.sets.find(node);
            let next_base = self.sets.find(next);
            if self.even_time[base] < self.even_time[next_base] {
                self.contract(
                    base,
                    next_base,
                    Bridge {
                        from: next,
                        to: node,
                        edge,
                    },
                );
            }
        }
        None
    }

    fn contract(&mut self, base: NodeId, mut vertex: NodeId, bridge: Bridge) {
        while vertex != base {
            let odd = self.mate[vertex];
            self.sets.link(vertex);
            self.sets.link(odd);
            self.path.bridge[odd] = Some(bridge);
            // Descendant-to-ancestor insertion makes the closest newly even
            // vertex the first one visited by the explicit DFS stack.
            self.dfs.push((odd, 0));
            vertex = self.sets.find(self.path.incoming[odd].node);
        }
    }
}
