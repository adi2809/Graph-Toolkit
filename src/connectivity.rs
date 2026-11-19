//! Graph structure and minimum spanning forests.

use std::collections::VecDeque;

use crate::traversal::depth_first;
use crate::{Arc, EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

mod analysis;
mod cycles;
mod dag;
mod disjoint;
mod dominance;
mod dynamic;
mod spanning;
pub use analysis::{
    BiconnectedComponents, Condensation, EdgeCut, biconnected_components, condensation,
    connected_components, edge_connectivity, is_connected, is_strongly_connected,
    is_weakly_connected, minimum_st_edge_cut, st_connected,
};
pub use cycles::{
    EdgeWalk, cycle_basis, eulerian_circuit, eulerian_path, has_eulerian_path, is_eulerian,
};
pub use dag::{
    is_directed_acyclic_graph, lexicographical_topological_sort, reverse_topological_sort,
    topological_generations, transitive_closure_dag, transitive_reduction,
};
pub use disjoint::{edge_disjoint_paths, minimum_st_vertex_cut, node_disjoint_paths};
pub use dominance::{dominance_frontiers, immediate_dominators};
pub use dynamic::{
    ConnectivityAnswer, ConnectivityEvent, DynamicConnectivityError, offline_dynamic_connectivity,
};
pub use spanning::{
    BottleneckTree, SpanningTree, minimum_bottleneck_arborescence,
    minimum_bottleneck_spanning_tree, minimum_spanning_tree,
};

#[derive(Clone, Debug)]
pub struct StrongComponents {
    /// Component ID for each vertex. IDs are in topological order in the condensation DAG.
    pub component_of: Vec<usize>,
    pub count: usize,
}

/// Iterative Kosaraju decomposition using the graph's precomputed transpose.
pub fn strongly_connected_components<W>(graph: &Graph<W>) -> StrongComponents {
    let n = graph.node_count();
    let mut visited = vec![false; n];
    let mut stack = Vec::new();
    let mut finish = Vec::with_capacity(n);
    for start in 0..n {
        if !visited[start] {
            depth_first(
                graph,
                start,
                false,
                &mut visited,
                &mut stack,
                |node, exiting| {
                    if exiting {
                        finish.push(node);
                    }
                },
            );
        }
    }

    visited.fill(false);
    let mut component_of = vec![NONE; n];
    let mut count = 0;
    for start in finish.into_iter().rev() {
        if !visited[start] {
            depth_first(graph, start, true, &mut visited, &mut stack, |node, _| {
                component_of[node] = count;
            });
            count += 1;
        }
    }
    StrongComponents {
        component_of,
        count,
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct WeakComponents {
    /// Component ID for each vertex. IDs follow the lowest vertex in each component.
    pub component_of: Vec<usize>,
    pub count: usize,
}

/// Connected components when edge orientation is ignored, in O(V + E) time.
pub fn weakly_connected_components<W>(graph: &Graph<W>) -> WeakComponents {
    let n = graph.node_count();
    let mut result = WeakComponents {
        component_of: vec![NONE; n],
        count: 0,
    };
    let mut queue = VecDeque::new();
    for root in 0..n {
        if result.component_of[root] != NONE {
            continue;
        }
        result.component_of[root] = result.count;
        queue.push_back(root);
        while let Some(node) = queue.pop_front() {
            for arc in incident(graph, node) {
                if result.component_of[arc.node] == NONE {
                    result.component_of[arc.node] = result.count;
                    queue.push_back(arc.node);
                }
            }
        }
        result.count += 1;
    }
    result
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct UndirectedCuts {
    /// Logical bridge edge IDs in insertion order.
    pub bridges: Vec<EdgeId>,
    /// Cut vertices in ascending ID order.
    pub articulation_points: Vec<NodeId>,
}

/// Finds all bridges and articulation points in one iterative O(V + E) low-link DFS.
/// A single parallel edge never disconnects its endpoints, and loops are ignored.
pub fn bridges_and_articulation_points<W>(graph: &Graph<W>) -> Result<UndirectedCuts, GraphError> {
    graph.require(GraphKind::Undirected)?;
    Ok(undirected_low_links(graph, false).0)
}

fn undirected_low_links<W>(
    graph: &Graph<W>,
    collect_blocks: bool,
) -> (UndirectedCuts, Vec<Vec<EdgeId>>) {
    let n = graph.node_count();
    let mut pending = Vec::new();
    let mut blocks = Vec::new();
    let mut seen = vec![NONE; n];
    let mut low = vec![0; n];
    let mut bridges = vec![false; graph.edge_count()];
    let mut cuts = vec![false; n];
    let mut stack = Vec::new();
    let mut time = 0;
    for root in 0..n {
        if seen[root] != NONE {
            continue;
        }
        seen[root] = time;
        low[root] = time;
        time += 1;
        stack.push(CutFrame::new(root, NONE));
        while let Some(frame) = stack.last_mut() {
            let node = frame.node;
            if let Some(arc) = graph.out(node).get(frame.cursor) {
                frame.cursor += 1;
                // Skip the exact tree edge. A parallel edge is a back edge.
                if arc.edge == frame.parent_edge {
                    continue;
                }
                if collect_blocks {
                    if arc.node == node {
                        blocks.push(vec![arc.edge]);
                    } else if seen[arc.node] < seen[node] || seen[arc.node] == NONE {
                        pending.push(arc.edge);
                    }
                }
                if seen[arc.node] != NONE {
                    low[node] = low[node].min(seen[arc.node]);
                } else {
                    let child = arc.node;
                    frame.children += 1;
                    seen[child] = time;
                    low[child] = time;
                    time += 1;
                    stack.push(CutFrame::new(child, arc.edge));
                }
                continue;
            }
            let child = stack.pop().unwrap();
            if let Some(parent) = stack.last() {
                let p = parent.node;
                low[p] = low[p].min(low[node]);
                if low[node] > seen[p] {
                    bridges[child.parent_edge] = true;
                }
                if low[node] >= seen[p] {
                    if parent.parent_edge != NONE {
                        cuts[p] = true;
                    }
                    if collect_blocks {
                        let mut block = Vec::new();
                        loop {
                            let edge = pending.pop().unwrap();
                            block.push(edge);
                            if edge == child.parent_edge {
                                break;
                            }
                        }
                        blocks.push(block);
                    }
                }
            } else {
                cuts[node] = child.children > 1;
            }
        }
    }
    (
        UndirectedCuts {
            bridges: marked_ids(bridges),
            articulation_points: marked_ids(cuts),
        },
        blocks,
    )
}

struct CutFrame {
    node: NodeId,
    parent_edge: EdgeId,
    cursor: usize,
    children: usize,
}

impl CutFrame {
    fn new(node: NodeId, parent_edge: EdgeId) -> Self {
        Self {
            node,
            parent_edge,
            cursor: 0,
            children: 0,
        }
    }
}

fn marked_ids(marked: Vec<bool>) -> Vec<usize> {
    marked
        .into_iter()
        .enumerate()
        .filter_map(|(id, yes)| yes.then_some(id))
        .collect()
}

/// Batagelj–Zaveršnik core decomposition in O(V + E) time.
/// A vertex's core number is the largest k for which it belongs to a k-core.
/// Parallel edges contribute separately, self-loops are ignored, and directed
/// degree is the sum of in-degree and out-degree.
pub fn core_numbers<W>(graph: &Graph<W>) -> Vec<usize> {
    let n = graph.node_count();
    let mut degree: Vec<_> = (0..n)
        .map(|node| incident(graph, node).filter(|arc| arc.node != node).count())
        .collect();
    let mut bins = vec![0; degree.iter().copied().max().unwrap_or(0) + 1];
    for &d in &degree {
        bins[d] += 1;
    }
    let mut start = 0;
    for bin in &mut bins {
        let count = *bin;
        *bin = start;
        start += count;
    }
    let mut cursor = bins.clone();
    let mut position = vec![0; n];
    let mut order = vec![0; n];
    for (node, &d) in degree.iter().enumerate() {
        let slot = cursor[d];
        position[node] = slot;
        order[slot] = node;
        cursor[d] += 1;
    }
    for index in 0..n {
        let node = order[index];
        for arc in incident(graph, node) {
            let neighbor = arc.node;
            let d = degree[neighbor];
            if d > degree[node] {
                let first = bins[d];
                let at = position[neighbor];
                position[order[first]] = at;
                position[neighbor] = first;
                order.swap(first, at);
                bins[d] += 1;
                degree[neighbor] -= 1;
            }
        }
    }
    degree
}

fn incident<W>(graph: &Graph<W>, node: NodeId) -> impl Iterator<Item = &Arc> {
    let incoming = if graph.kind() == GraphKind::Directed {
        graph.inc(node)
    } else {
        &[]
    };
    graph.out(node).iter().chain(incoming)
}

/// Kahn's O(V + E) topological ordering. A self-loop is a cycle.
pub fn topological_sort<W>(graph: &Graph<W>) -> Result<Vec<NodeId>, GraphError> {
    graph.require(GraphKind::Directed)?;
    let n = graph.node_count();
    let mut remaining: Vec<_> = (0..n).map(|node| graph.inc(node).len()).collect();
    let mut queue: VecDeque<_> = (0..n).filter(|&node| remaining[node] == 0).collect();
    let mut order = Vec::with_capacity(n);
    while let Some(node) = queue.pop_front() {
        order.push(node);
        for arc in graph.out(node) {
            remaining[arc.node] -= 1;
            if remaining[arc.node] == 0 {
                queue.push_back(arc.node);
            }
        }
    }
    if order.len() != n {
        return Err(GraphError::Cycle);
    }
    Ok(order)
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SpanningForest {
    pub edges: Vec<EdgeId>,
    pub total_weight: u128,
    pub component_count: usize,
}

/// Kruskal's minimum spanning forest with path halving and union by size.
pub fn minimum_spanning_forest(graph: &Graph<u64>) -> Result<SpanningForest, GraphError> {
    let (edges, component_count) = kruskal(graph)?;
    let total_weight = edges
        .iter()
        .map(|&id| u128::from(graph.edges()[id].weight))
        .sum();
    Ok(SpanningForest {
        edges,
        total_weight,
        component_count,
    })
}

fn kruskal<W: Ord + Copy>(graph: &Graph<W>) -> Result<(Vec<EdgeId>, usize), GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    let weights = graph.edges();
    let mut order: Vec<_> = weights
        .iter()
        .enumerate()
        .filter_map(|(id, edge)| (edge.from != edge.to).then_some(id))
        .collect();
    order.sort_unstable_by_key(|&id| (weights[id].weight, id));
    let mut sets = DisjointSet::new(n);
    let mut selected = Vec::new();
    let mut components = n;
    for id in order {
        let edge = weights[id];
        if sets.union(edge.from, edge.to).is_some() {
            selected.push(id);
            components -= 1;
            if components == 1 {
                break;
            }
        }
    }
    Ok((selected, components))
}

/// Dense union–find with union by size and iterative path halving.
/// Stores one word per element. Operations take amortized O(α(V)) time.
/// Vertex arguments are indexed like slices and panic if out of bounds.
#[derive(Clone, Debug)]
pub struct DisjointSet {
    // A root stores the negative size, all other entries store parent IDs.
    parent: Vec<isize>,
    count: usize,
}

impl DisjointSet {
    pub fn new(n: usize) -> Self {
        Self {
            parent: vec![-1; n],
            count: n,
        }
    }

    pub fn len(&self) -> usize {
        self.parent.len()
    }

    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    pub fn component_count(&self) -> usize {
        self.count
    }

    pub fn find(&mut self, mut node: usize) -> usize {
        loop {
            let parent = self.parent[node];
            if parent < 0 {
                return node;
            }
            let grandparent = self.parent[parent as usize];
            if grandparent < 0 {
                return parent as usize;
            }
            self.parent[node] = grandparent;
            node = grandparent as usize;
        }
    }

    pub fn connected(&mut self, a: usize, b: usize) -> bool {
        self.find(a) == self.find(b)
    }

    pub fn component_size(&mut self, vertex: usize) -> usize {
        let root = self.find(vertex);
        (-self.parent[root]) as usize
    }

    /// Merges two sets and returns the surviving representative. Returns
    /// `None` when both vertices already belong to the same set.
    pub fn union(&mut self, a: usize, b: usize) -> Option<usize> {
        let (mut a, mut b) = (self.find(a), self.find(b));
        if a == b {
            return None;
        }
        if self.parent[a] > self.parent[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.parent[a] += self.parent[b];
        self.parent[b] = a as isize;
        self.count -= 1;
        Some(a)
    }
    /// Internal blossom search resets every member of each affected set.
    pub(crate) fn reset(&mut self, node: usize) {
        self.count += usize::from(self.parent[node] >= 0);
        self.parent[node] = -1;
    }
}
