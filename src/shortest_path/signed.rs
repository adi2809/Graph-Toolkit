//! Signed shortest paths with original edge IDs and cycle witnesses.

use std::error::Error;
use std::fmt;

use crate::connectivity::{strongly_connected_components, topological_sort};
use crate::queue::RadixHeap;
use crate::{EdgeId, Graph, GraphError, NodeId};

use super::{INF, NO_LINK, dijkstra_core, trace_path};

pub(super) const UNREACHABLE: i128 = i128::MAX;
type Previous = (NodeId, EdgeId);

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SignedPath {
    pub cost: i128,
    pub nodes: Vec<NodeId>,
    pub edges: Vec<EdgeId>,
}

/// A closed, oriented walk of negative total cost. In an undirected graph a
/// negative edge can appear twice, once in each direction.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NegativeCycle {
    pub cost: i128,
    pub nodes: Vec<NodeId>,
    pub edges: Vec<EdgeId>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum SignedPathError {
    Graph(GraphError),
    NegativeCycle(NegativeCycle),
    MatrixAllocation { node_count: usize },
}

impl fmt::Display for SignedPathError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Graph(error) => error.fmt(f),
            Self::NegativeCycle(cycle) => write!(f, "negative cycle of cost {}", cycle.cost),
            Self::MatrixAllocation { node_count } => {
                write!(f, "could not allocate matrices for {node_count} vertices")
            }
        }
    }
}

impl Error for SignedPathError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Graph(error) => Some(error),
            _ => None,
        }
    }
}

impl From<GraphError> for SignedPathError {
    fn from(error: GraphError) -> Self {
        Self::Graph(error)
    }
}

#[derive(Clone, Debug)]
pub struct SignedShortestPaths {
    pub(super) distance: Vec<i128>,
    pub(super) predecessor: Vec<Previous>,
}

impl SignedShortestPaths {
    pub fn distance(&self, node: NodeId) -> Option<i128> {
        finite(*self.distance.get(node)?)
    }

    /// Iterates over all distances in vertex order.
    pub fn distances(&self) -> impl ExactSizeIterator<Item = Option<i128>> + '_ {
        self.distance.iter().copied().map(finite)
    }

    pub fn path_to(&self, node: NodeId) -> Option<SignedPath> {
        Some(path(self.distance(node)?, &self.predecessor, node))
    }
}

/// Row-major distances and predecessors. Both matrices use O(V²) space.
#[derive(Clone, Debug)]
pub struct AllPairsShortestPaths {
    n: usize,
    distance: Vec<i128>,
    predecessor: Vec<Previous>,
}

impl AllPairsShortestPaths {
    pub fn node_count(&self) -> usize {
        self.n
    }

    pub fn distance(&self, from: NodeId, to: NodeId) -> Option<i128> {
        if from >= self.n || to >= self.n {
            return None;
        }
        finite(self.distance[from * self.n + to])
    }

    /// Iterates over one source row in target vertex order.
    pub fn distances_from(
        &self,
        from: NodeId,
    ) -> Option<impl ExactSizeIterator<Item = Option<i128>> + '_> {
        if from >= self.n {
            return None;
        }
        let start = from * self.n;
        Some(
            self.distance[start..start + self.n]
                .iter()
                .copied()
                .map(finite),
        )
    }

    pub fn path(&self, from: NodeId, to: NodeId) -> Option<SignedPath> {
        let cost = self.distance(from, to)?;
        let start = from * self.n;
        Some(path(cost, &self.predecessor[start..start + self.n], to))
    }

    fn new(n: usize) -> Result<Self, SignedPathError> {
        let allocation_error = SignedPathError::MatrixAllocation { node_count: n };
        let len = n.checked_mul(n).ok_or_else(|| allocation_error.clone())?;
        let mut distance = Vec::new();
        let mut predecessor = Vec::new();
        distance
            .try_reserve_exact(len)
            .map_err(|_| allocation_error.clone())?;
        predecessor
            .try_reserve_exact(len)
            .map_err(|_| allocation_error)?;
        distance.resize(len, UNREACHABLE);
        predecessor.resize(len, NO_LINK);
        Ok(Self {
            n,
            distance,
            predecessor,
        })
    }
}

fn finite(value: i128) -> Option<i128> {
    (value != UNREACHABLE).then_some(value)
}

fn path(cost: i128, predecessor: &[Previous], node: NodeId) -> SignedPath {
    let (nodes, edges) = trace_path(predecessor, node);
    SignedPath { cost, nodes, edges }
}

/// Single-source Bellman–Ford. Reports a negative cycle reachable from `source`.
/// SCCs are processed in topological order, with active-frontier relaxation
/// inside each component. Runs in O(V + E) on DAGs and O(V + VE) in the worst case.
pub fn bellman_ford(
    graph: &Graph<i64>,
    source: NodeId,
) -> Result<SignedShortestPaths, SignedPathError> {
    graph.validate_node(source)?;
    let mut distance = vec![UNREACHABLE; graph.node_count()];
    distance[source] = 0;
    Ok(relax_by_component(graph, distance, usize::MAX).unwrap()?)
}

/// Finds a negative-cycle witness anywhere in the graph.
pub fn find_negative_cycle(graph: &Graph<i64>) -> Option<NegativeCycle> {
    if graph.edges().iter().all(|edge| edge.weight >= 0) {
        return None;
    }
    relax_by_component(graph, vec![0; graph.node_count()], usize::MAX)
        .unwrap()
        .err()
}

impl From<NegativeCycle> for SignedPathError {
    fn from(cycle: NegativeCycle) -> Self {
        Self::NegativeCycle(cycle)
    }
}

pub(super) fn relax_by_component(
    graph: &Graph<i64>,
    mut distance: Vec<i128>,
    limit: usize,
) -> Option<Result<SignedShortestPaths, NegativeCycle>> {
    let n = graph.node_count();
    let components = strongly_connected_components(graph);
    let mut offsets = vec![0; components.count + 1];
    for &id in &components.component_of {
        offsets[id + 1] += 1;
    }
    for id in 0..components.count {
        offsets[id + 1] += offsets[id];
    }
    let mut nodes = vec![0; n];
    let mut cursor = offsets[..components.count].to_vec();
    for (node, &id) in components.component_of.iter().enumerate() {
        nodes[cursor[id]] = node;
        cursor[id] += 1;
    }

    let mut predecessor = vec![NO_LINK; n];
    let mut pending = vec![false; n];
    let mut current = Vec::new();
    let mut next = Vec::new();
    for (id, range) in offsets.windows(2).enumerate() {
        let members = &nodes[range[0]..range[1]];
        current.clear();
        current.extend(
            members
                .iter()
                .filter_map(|&node| finite(distance[node]).map(|cost| (node, cost))),
        );
        if current.is_empty() {
            continue;
        }
        for round in 0..members.len() {
            // Each scan uses its frontier label from the preceding round.
            for &(node, cost) in &current {
                for arc in graph.out(node) {
                    if components.component_of[arc.node] != id {
                        continue;
                    }
                    let candidate = cost + i128::from(graph.edges()[arc.edge].weight);
                    if candidate < distance[arc.node] {
                        distance[arc.node] = candidate;
                        predecessor[arc.node] = (node, arc.edge);
                        if !pending[arc.node] {
                            pending[arc.node] = true;
                            next.push(arc.node);
                        }
                    }
                }
            }
            let Some(&updated) = next.last() else {
                break;
            };
            if round + 1 == members.len() {
                // If the final predecessor chain were simple, it would use
                // at most size − 1 internal edges and could not improve now.
                return Some(Err(cycle_from_predecessors(
                    graph,
                    &predecessor,
                    updated,
                    members.len(),
                )));
            }
            if round + 1 == limit {
                return None;
            }
            current.clear();
            for node in next.drain(..) {
                pending[node] = false;
                current.push((node, distance[node]));
            }
        }

        // The condensation DAG cannot revisit this component. Cross edges
        // need exactly one scan after its internal shortest paths are final.
        for &node in members {
            if distance[node] == UNREACHABLE {
                continue;
            }
            for arc in graph.out(node) {
                if components.component_of[arc.node] != id {
                    let candidate = distance[node] + i128::from(graph.edges()[arc.edge].weight);
                    if candidate < distance[arc.node] {
                        distance[arc.node] = candidate;
                        predecessor[arc.node] = (node, arc.edge);
                    }
                }
            }
        }
    }
    Some(Ok(SignedShortestPaths {
        distance,
        predecessor,
    }))
}

fn cycle_from_predecessors(
    graph: &Graph<i64>,
    predecessor: &[Previous],
    mut vertex: NodeId,
    size: usize,
) -> NegativeCycle {
    for _ in 0..size {
        vertex = predecessor[vertex].0;
    }
    let start = vertex;
    let mut nodes = vec![start];
    let mut edges = Vec::new();
    let mut cost = 0;
    loop {
        let (previous, edge) = predecessor[vertex];
        edges.push(edge);
        nodes.push(previous);
        cost += i128::from(graph.edges()[edge].weight);
        if previous == start {
            // On a zero-cost predecessor cycle, every edge would have
            // used the final label of a predecessor from an earlier round.
            // The strict ordering of those rounds rules out such a cycle.
            debug_assert!(cost < 0);
            nodes.reverse();
            edges.reverse();
            return NegativeCycle { cost, nodes, edges };
        }
        vertex = previous;
    }
}

/// Johnson all-pairs shortest paths. Bellman–Ford computes feasible potentials,
/// then each source uses indexed radix-heap Dijkstra on nonnegative reduced costs.
/// Any negative cycle is reported before allocating the result matrices.
pub fn johnson(graph: &Graph<i64>) -> Result<AllPairsShortestPaths, SignedPathError> {
    let n = graph.node_count();
    let potential = if graph.edges().iter().any(|edge| edge.weight < 0) {
        relax_by_component(graph, vec![0; n], usize::MAX)
            .unwrap()?
            .distance
    } else {
        vec![0; n]
    };
    let reduced: Vec<_> = graph
        .edges()
        .iter()
        .map(|edge| {
            let weight = i128::from(edge.weight) + potential[edge.from] - potential[edge.to];
            debug_assert!(weight >= 0);
            weight as u128
        })
        .collect();
    let mut result = AllPairsShortestPaths::new(n)?;
    let mut distance = vec![INF; n];
    let mut queue = RadixHeap::new(n);
    for source in 0..n {
        let start = source * n;
        let parents = &mut result.predecessor[start..start + n];
        dijkstra_core(
            graph,
            source,
            false,
            |arc| reduced[arc.edge],
            &mut distance,
            &mut queue,
            |node, previous, edge| parents[node] = (previous, edge),
        );
        for (target, &cost) in distance.iter().enumerate() {
            if cost != INF {
                result.distance[start + target] =
                    cost as i128 - potential[source] + potential[target];
            }
        }
    }
    Ok(result)
}

/// Cache-blocked Floyd–Warshall. Reports any negative cycle, including one in
/// a disconnected component. The distance kernel takes O(V³ + E) time.
pub fn floyd_warshall(graph: &Graph<i64>) -> Result<AllPairsShortestPaths, SignedPathError> {
    const BLOCK: usize = 32;
    let n = graph.node_count();
    let mut result = AllPairsShortestPaths::new(n)?;
    let mut hops = Vec::new();
    hops.try_reserve_exact(result.distance.len())
        .map_err(|_| SignedPathError::MatrixAllocation { node_count: n })?;
    hops.resize(result.distance.len(), u32::MAX);
    for node in 0..n {
        result.distance[node * n + node] = 0;
        hops[node * n + node] = 0;
        for arc in graph.out(node) {
            let entry = node * n + arc.node;
            let cost = i128::from(graph.edges()[arc.edge].weight);
            if cost < result.distance[entry] {
                result.distance[entry] = cost;
                result.predecessor[entry] = (node, arc.edge);
                hops[entry] = 1;
            }
        }
    }
    if negative_diagonal(&result, 0, n) {
        return Err(floyd_cycle(graph));
    }
    for via in (0..n).step_by(BLOCK) {
        update_tile(&mut result, &mut hops, via, via, via, BLOCK);
        if negative_diagonal(&result, via, (via + BLOCK).min(n)) {
            return Err(floyd_cycle(graph));
        }
        // Close the pivot tile, then update its row and column before
        // using either to relax tiles outside the pivot row and column.
        for other in (0..n).step_by(BLOCK).filter(|&index| index != via) {
            update_tile(&mut result, &mut hops, via, via, other, BLOCK);
            update_tile(&mut result, &mut hops, via, other, via, BLOCK);
        }
        for from in (0..n).step_by(BLOCK).filter(|&index| index != via) {
            for to in (0..n).step_by(BLOCK).filter(|&index| index != via) {
                update_tile(&mut result, &mut hops, via, from, to, BLOCK);
            }
            if negative_diagonal(&result, from, (from + BLOCK).min(n)) {
                return Err(floyd_cycle(graph));
            }
        }
    }
    Ok(result)
}

fn negative_diagonal(data: &AllPairsShortestPaths, start: usize, end: usize) -> bool {
    (start..end).any(|node| data.distance[node * data.n + node] < 0)
}

fn floyd_cycle(graph: &Graph<i64>) -> SignedPathError {
    find_negative_cycle(graph)
        .expect("a negative diagonal certifies a negative cycle")
        .into()
}

fn update_tile(
    result: &mut AllPairsShortestPaths,
    hops: &mut [u32],
    via: usize,
    from: usize,
    to: usize,
    block: usize,
) {
    let n = result.n;
    for k in via..(via + block).min(n) {
        let row_k = k * n;
        for i in from..(from + block).min(n) {
            let row_i = i * n;
            let left = result.distance[row_i + k];
            if left == UNREACHABLE {
                continue;
            }
            let left_hops = hops[row_i + k];
            for j in to..(to + block).min(n) {
                let right = result.distance[row_k + j];
                if right != UNREACHABLE {
                    // With a negative cycle, labels can grow in magnitude
                    // before the diagonal test. Saturation preserves its sign.
                    let candidate = left.saturating_add(right);
                    let entry = row_i + j;
                    if candidate <= result.distance[entry] {
                        // Among equal-cost walks use the shortest in edge count.
                        // This keeps predecessor rows acyclic around zero cycles.
                        let count = left_hops.saturating_add(hops[row_k + j]);
                        if candidate < result.distance[entry] || count < hops[entry] {
                            result.distance[entry] = candidate;
                            result.predecessor[entry] = result.predecessor[row_k + j];
                            hops[entry] = count;
                        }
                    }
                }
            }
        }
    }
}

/// Reuses a topological order across signed routing queries on the same DAG.
#[derive(Clone, Debug)]
pub struct DagShortestPathIndex<'g> {
    graph: &'g Graph<i64>,
    order: Vec<NodeId>,
}

impl<'g> DagShortestPathIndex<'g> {
    /// Validates acyclicity and caches a topological order in O(V + E).
    pub fn new(graph: &'g Graph<i64>) -> Result<Self, GraphError> {
        Ok(Self {
            graph,
            order: topological_sort(graph)?,
        })
    }

    /// Exact shortest routes from the closest source in O(V + E). A source
    /// can itself have a negative-cost route from another source.
    pub fn from_sources(
        &self,
        sources: impl IntoIterator<Item = NodeId>,
    ) -> Result<SignedShortestPaths, GraphError> {
        let mut distance = vec![UNREACHABLE; self.order.len()];
        for source in sources {
            self.graph.validate_node(source)?;
            distance[source] = 0;
        }
        Ok(self.relax(distance, None))
    }

    /// Exact minimum-cost route between two vertices. Stops when the target
    /// is reached in topological order.
    pub fn between(
        &self,
        source: NodeId,
        target: NodeId,
    ) -> Result<Option<SignedPath>, GraphError> {
        self.graph.validate_node(source)?;
        self.graph.validate_node(target)?;
        let mut distance = vec![UNREACHABLE; self.order.len()];
        distance[source] = 0;
        Ok(self.relax(distance, Some(target)).path_to(target))
    }

    /// Globally cheapest path. Single vertices have cost zero. An empty DAG
    /// has no path.
    pub fn minimum_path(&self) -> Option<SignedPath> {
        let paths = self.relax(vec![0; self.order.len()], None);
        let (end, _) = paths
            .distance
            .iter()
            .enumerate()
            .min_by_key(|&(node, &cost)| (cost, node))?;
        paths.path_to(end)
    }

    fn relax(&self, mut distance: Vec<i128>, target: Option<NodeId>) -> SignedShortestPaths {
        let mut predecessor = vec![NO_LINK; self.order.len()];
        for &node in &self.order {
            if Some(node) == target {
                break;
            }
            if distance[node] == UNREACHABLE {
                continue;
            }
            for arc in self.graph.out(node) {
                let candidate = distance[node] + i128::from(self.graph.edges()[arc.edge].weight);
                if candidate < distance[arc.node] {
                    distance[arc.node] = candidate;
                    predecessor[arc.node] = (node, arc.edge);
                }
            }
        }
        SignedShortestPaths {
            distance,
            predecessor,
        }
    }
}

/// Single-source signed shortest paths in a DAG in O(V + E) time.
pub fn dag_shortest_paths(
    graph: &Graph<i64>,
    source: NodeId,
) -> Result<SignedShortestPaths, GraphError> {
    graph.validate_node(source)?;
    DagShortestPathIndex::new(graph)?.from_sources([source])
}

/// Minimum-cost route between two vertices of a signed DAG.
pub fn dag_minimum_cost_path(
    graph: &Graph<i64>,
    source: NodeId,
    target: NodeId,
) -> Result<Option<SignedPath>, GraphError> {
    DagShortestPathIndex::new(graph)?.between(source, target)
}

/// Maximum-total-weight path in a directed acyclic graph in O(V + E) time.
/// A single vertex is a path of weight zero. An empty graph returns `None`.
pub fn dag_longest_path(graph: &Graph<i64>) -> Result<Option<SignedPath>, GraphError> {
    let order = topological_sort(graph)?;
    let mut distance = vec![0; graph.node_count()];
    let mut predecessor = vec![NO_LINK; graph.node_count()];
    for node in order {
        for arc in graph.out(node) {
            let candidate = distance[node] + i128::from(graph.edges()[arc.edge].weight);
            if candidate > distance[arc.node] {
                distance[arc.node] = candidate;
                predecessor[arc.node] = (node, arc.edge);
            }
        }
    }
    Ok(distance
        .iter()
        .copied()
        .enumerate()
        .max_by_key(|&(node, cost)| (cost, std::cmp::Reverse(node)))
        .map(|(node, cost)| path(cost, &predecessor, node)))
}
