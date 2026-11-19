//! Shortest paths on signed and nonnegative integer weights, and bottleneck routes.

use std::cmp::Reverse;
use std::collections::VecDeque;

use crate::queue::{IndexedHeap, RadixHeap};
use crate::{Arc, EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

mod bottleneck;
mod duong;
mod signed;
mod yen;
pub(crate) use bottleneck::bottleneck_search;
pub use bottleneck::{BottleneckPath, minimum_bottleneck_path};
pub use duong::duong_shortest_paths;
pub use signed::{
    AllPairsShortestPaths, DagShortestPathIndex, NegativeCycle, SignedPath, SignedPathError,
    SignedShortestPaths, bellman_ford, dag_longest_path, dag_minimum_cost_path, dag_shortest_paths,
    find_negative_cycle, floyd_warshall, johnson,
};
pub use yen::k_shortest_simple_paths;

/// Distances use a wider type than edge weights to accommodate long paths.
pub type Cost = u128;
const INF: Cost = Cost::MAX;
const NO_LINK: (NodeId, EdgeId) = (NONE, NONE);

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Path {
    pub cost: Cost,
    pub nodes: Vec<NodeId>,
    /// Logical edge IDs, in traversal order, so parallel edges are unambiguous.
    pub edges: Vec<EdgeId>,
}

#[derive(Clone, Debug)]
pub struct ShortestPaths {
    distance: Vec<Cost>,
    parent: Vec<(NodeId, EdgeId)>,
}

impl ShortestPaths {
    fn new(n: usize) -> Self {
        Self {
            distance: vec![INF; n],
            parent: vec![NO_LINK; n],
        }
    }

    pub fn distance(&self, node: NodeId) -> Option<Cost> {
        self.distance.get(node).copied().filter(|&cost| cost != INF)
    }

    pub fn path_to(&self, node: NodeId) -> Option<Path> {
        let cost = self.distance(node)?;
        let (nodes, edges) = trace_path(&self.parent, node);
        Some(Path { cost, nodes, edges })
    }
}

fn trace_path(parent: &[(NodeId, EdgeId)], mut node: NodeId) -> (Vec<NodeId>, Vec<EdgeId>) {
    let mut nodes = vec![node];
    let mut edges = Vec::new();
    while parent[node] != NO_LINK {
        debug_assert!(
            nodes.len() < parent.len(),
            "a shortest path cannot repeat vertices"
        );
        let (previous, edge) = parent[node];
        node = previous;
        nodes.push(node);
        edges.push(edge);
    }
    nodes.reverse();
    edges.reverse();
    (nodes, edges)
}

/// Single-source Dijkstra with an indexed monotone radix heap and O(V) auxiliary space.
pub fn dijkstra(graph: &Graph<u64>, source: NodeId) -> Result<ShortestPaths, GraphError> {
    graph.validate_node(source)?;
    let mut parent = vec![NO_LINK; graph.node_count()];
    let distance = dijkstra_search(graph, source, false, |node, previous, edge| {
        parent[node] = (previous, edge);
    });
    Ok(ShortestPaths { distance, parent })
}

fn dijkstra_search(
    graph: &Graph<u64>,
    source: NodeId,
    reverse: bool,
    relax: impl FnMut(NodeId, NodeId, EdgeId),
) -> Vec<Cost> {
    let mut distance = vec![INF; graph.node_count()];
    let mut queue = RadixHeap::new(graph.node_count());
    dijkstra_core(
        graph,
        source,
        reverse,
        |arc| Cost::from(graph.edges()[arc.edge].weight),
        &mut distance,
        &mut queue,
        relax,
    );
    distance
}

fn dijkstra_core<W>(
    graph: &Graph<W>,
    source: NodeId,
    reverse: bool,
    weight: impl Fn(Arc) -> Cost,
    distance: &mut [Cost],
    queue: &mut RadixHeap,
    mut relax: impl FnMut(NodeId, NodeId, EdgeId),
) {
    distance.fill(INF);
    distance[source] = 0;
    queue.reset();
    queue.decrease(source, INF, 0);
    while let Some((cost, node)) = queue.pop(distance) {
        let arcs = if reverse {
            graph.inc(node)
        } else {
            graph.out(node)
        };
        for arc in arcs {
            let candidate = cost + weight(*arc);
            if candidate < distance[arc.node] {
                queue.decrease(arc.node, distance[arc.node], candidate);
                distance[arc.node] = candidate;
                relax(arc.node, node, arc.edge);
            }
        }
    }
}

fn start_frontier(distance: &mut [Cost], source: NodeId) -> RadixHeap {
    distance[source] = 0;
    let mut queue = RadixHeap::new(distance.len());
    queue.decrease(source, INF, 0);
    queue
}

/// O(V + E) shortest paths when every graph edge has weight zero or one.
pub fn zero_one_bfs(graph: &Graph<u64>, source: NodeId) -> Result<ShortestPaths, GraphError> {
    graph.validate_node(source)?;
    for (edge, item) in graph.edges().iter().enumerate() {
        if item.weight > 1 {
            return Err(GraphError::NonBinaryWeight {
                edge,
                weight: item.weight,
            });
        }
    }
    let mut paths = ShortestPaths::new(graph.node_count());
    paths.distance[source] = 0;
    let mut deque = VecDeque::from([(0, source)]);
    while let Some((cost, node)) = deque.pop_front() {
        if cost != paths.distance[node] {
            continue;
        }
        for arc in graph.out(node) {
            let weight = graph.edges()[arc.edge].weight;
            let candidate = cost + Cost::from(weight);
            if candidate >= paths.distance[arc.node] {
                continue;
            }
            paths.distance[arc.node] = candidate;
            paths.parent[arc.node] = (node, arc.edge);
            if weight == 0 {
                deque.push_front((candidate, arc.node));
            } else {
                deque.push_back((candidate, arc.node));
            }
        }
    }
    Ok(paths)
}

/// A* with an indexed 4-ary heap. The heuristic must be a stable, admissible
/// lower bound to the target and must return zero at the target.
/// Inconsistent admissible heuristics are supported by reopening vertices.
pub fn astar(
    graph: &Graph<u64>,
    source: NodeId,
    target: NodeId,
    mut heuristic: impl FnMut(NodeId) -> u64,
) -> Result<Option<Path>, GraphError> {
    graph_astar(graph, source, target, |node| Cost::from(heuristic(node)))
}

fn graph_astar(
    graph: &Graph<u64>,
    source: NodeId,
    target: NodeId,
    heuristic: impl FnMut(NodeId) -> Cost,
) -> Result<Option<Path>, GraphError> {
    graph.validate_node(source)?;
    graph.validate_node(target)?;
    Ok(astar_core(
        graph.node_count(),
        source,
        target,
        |node| {
            graph
                .out(node)
                .iter()
                .map(|arc| (arc.node, arc.edge, graph.edges()[arc.edge].weight))
        },
        heuristic,
    ))
}

pub(crate) fn astar_core<I: IntoIterator<Item = (NodeId, EdgeId, u64)>>(
    n: usize,
    source: NodeId,
    target: NodeId,
    neighbors: impl Fn(NodeId) -> I,
    mut heuristic: impl FnMut(NodeId) -> Cost,
) -> Option<Path> {
    if source == target {
        return Some(singleton_path(source));
    }
    let mut paths = ShortestPaths::new(n);
    let mut heap = IndexedHeap::new(n);
    paths.distance[source] = 0;
    heap.decrease(source, (heuristic(source), Reverse(0)));
    while let Some((_, node)) = heap.pop() {
        if node == target {
            return paths.path_to(target);
        }
        for (next, edge, weight) in neighbors(node) {
            let cost = paths.distance[node] + Cost::from(weight);
            if cost < paths.distance[next] {
                paths.distance[next] = cost;
                paths.parent[next] = (node, edge);
                heap.decrease(next, (cost + heuristic(next), Reverse(cost)));
            }
        }
    }
    None
}

/// Precomputes landmark distances for repeated exact A* queries (ALT).
/// Borrowing the graph keeps the distance tables tied to the graph that produced them.
#[derive(Debug)]
pub struct LandmarkRouter<'a> {
    graph: &'a Graph<u64>,
    landmarks: Vec<LandmarkDistances>,
}

#[derive(Debug)]
struct LandmarkDistances {
    from: Vec<Cost>,
    to: Option<Vec<Cost>>,
}

impl<'a> LandmarkRouter<'a> {
    /// Runs one single-source search per distinct landmark on undirected graphs,
    /// or one in each direction on directed graphs. Space is O(V × landmarks).
    /// An empty landmark list uses a zero heuristic.
    pub fn new(graph: &'a Graph<u64>, landmarks: &[NodeId]) -> Result<Self, GraphError> {
        for &node in landmarks {
            graph.validate_node(node)?;
        }
        let mut nodes = landmarks.to_vec();
        nodes.sort_unstable();
        nodes.dedup();
        let landmarks = nodes
            .into_iter()
            .map(|node| LandmarkDistances {
                from: dijkstra_search(graph, node, false, |_, _, _| {}),
                to: (graph.kind() == GraphKind::Directed)
                    .then(|| dijkstra_search(graph, node, true, |_, _, _| {})),
            })
            .collect();
        Ok(Self { graph, landmarks })
    }

    /// Finds an exact shortest path using triangle-inequality lower bounds.
    pub fn shortest_path(
        &self,
        source: NodeId,
        target: NodeId,
    ) -> Result<Option<Path>, GraphError> {
        graph_astar(self.graph, source, target, |node| {
            self.landmarks.iter().fold(0, |best, landmark| {
                let to = landmark.to.as_ref().unwrap_or(&landmark.from);
                best.max(landmark_bound(landmark.from[target], landmark.from[node]))
                    .max(landmark_bound(to[node], to[target]))
            })
        })
    }
}

fn landmark_bound(a: Cost, b: Cost) -> Cost {
    if a == INF || b == INF {
        0
    } else {
        a.saturating_sub(b)
    }
}

fn singleton_path(node: NodeId) -> Path {
    Path {
        cost: 0,
        nodes: vec![node],
        edges: Vec::new(),
    }
}

/// Two-frontier Dijkstra, with the reverse search following incoming arcs.
/// Stops when the sum of the frontier minima bounds the best complete route.
pub fn bidirectional_dijkstra(
    graph: &Graph<u64>,
    source: NodeId,
    target: NodeId,
) -> Result<Option<Path>, GraphError> {
    graph.validate_node(source)?;
    graph.validate_node(target)?;
    if source == target {
        return Ok(Some(singleton_path(source)));
    }

    let n = graph.node_count();
    let mut paths = [ShortestPaths::new(n), ShortestPaths::new(n)];
    let mut queues = [
        start_frontier(&mut paths[0].distance, source),
        start_frontier(&mut paths[1].distance, target),
    ];
    let mut best = (INF, NONE);
    let mut tie = 0;
    while let (Some(left), Some(right)) = (
        queues[0].peek_key(&paths[0].distance),
        queues[1].peek_key(&paths[1].distance),
    ) {
        if left.saturating_add(right) >= best.0 {
            break;
        }
        let side = match left.cmp(&right) {
            std::cmp::Ordering::Less => 0,
            std::cmp::Ordering::Greater => 1,
            std::cmp::Ordering::Equal => {
                let side = tie;
                tie ^= 1;
                side
            }
        };
        let (cost, node) = queues[side].pop(&paths[side].distance).unwrap();
        let arcs = if side == 0 {
            graph.out(node)
        } else {
            graph.inc(node)
        };
        for arc in arcs {
            let candidate = cost + Cost::from(graph.edges()[arc.edge].weight);
            if candidate < paths[side].distance[arc.node] {
                queues[side].decrease(arc.node, paths[side].distance[arc.node], candidate);
                paths[side].distance[arc.node] = candidate;
                paths[side].parent[arc.node] = (node, arc.edge);
                let combined = candidate.saturating_add(paths[side ^ 1].distance[arc.node]);
                if combined < best.0 {
                    best = (combined, arc.node);
                }
            }
        }
    }
    if best.1 == NONE {
        return Ok(None);
    }
    let mut node = best.1;
    let mut path = paths[0].path_to(node).unwrap();
    path.cost = best.0;
    while paths[1].parent[node] != NO_LINK {
        let (next, edge) = paths[1].parent[node];
        path.nodes.push(next);
        path.edges.push(edge);
        node = next;
    }
    Ok(Some(path))
}
