//! Stack-safe, insertion-ordered traversal and unweighted shortest paths.

use crate::{Graph, GraphError, NONE, NodeId};

/// The BFS order and forest, including the winning source of each reachable vertex.
#[derive(Clone, Debug)]
pub struct BreadthFirst {
    order: Vec<NodeId>,
    distance: Vec<usize>,
    parent: Vec<NodeId>,
    source: Vec<NodeId>,
}

impl BreadthFirst {
    pub fn order(&self) -> &[NodeId] {
        &self.order
    }

    pub fn distance(&self, node: NodeId) -> Option<usize> {
        present(&self.distance, node)
    }

    pub fn source(&self, node: NodeId) -> Option<NodeId> {
        present(&self.source, node)
    }

    pub fn path_to(&self, mut node: NodeId) -> Option<Vec<NodeId>> {
        self.distance(node)?;
        let mut path = vec![node];
        while self.parent[node] != NONE {
            node = self.parent[node];
            path.push(node);
        }
        path.reverse();
        Some(path)
    }
}

fn present(entries: &[usize], node: NodeId) -> Option<usize> {
    entries.get(node).copied().filter(|&value| value != NONE)
}

/// Returns the BFS discovery order from one vertex in O(V + E) time.
pub fn bfs<W>(graph: &Graph<W>, start: NodeId) -> Result<Vec<NodeId>, GraphError> {
    graph.validate_node(start)?;
    let mut visited = vec![false; graph.node_count()];
    Ok(bfs_order(
        &[start],
        |node| graph.out(node).iter().map(|arc| arc.node),
        |node, _| !std::mem::replace(&mut visited[node], true),
    ))
}

/// Finds the distance and path to a nearest source. Source order breaks ties.
/// Duplicate sources are accepted and an empty source list yields an empty forest.
pub fn multi_source_bfs<W>(
    graph: &Graph<W>,
    sources: &[NodeId],
) -> Result<BreadthFirst, GraphError> {
    for &source in sources {
        graph.validate_node(source)?;
    }
    Ok(bfs_core(graph.node_count(), sources, |node| {
        graph.out(node).iter().map(|arc| arc.node)
    }))
}

pub(crate) fn bfs_core<I: IntoIterator<Item = NodeId>>(
    n: usize,
    sources: &[NodeId],
    neighbors: impl Fn(NodeId) -> I,
) -> BreadthFirst {
    let mut distance = vec![NONE; n];
    let mut parent = vec![NONE; n];
    let mut source = vec![NONE; n];
    let order = bfs_order(sources, neighbors, |node, previous| {
        if distance[node] != NONE {
            return false;
        }
        if let Some(previous) = previous {
            distance[node] = distance[previous] + 1;
            parent[node] = previous;
            source[node] = source[previous];
        } else {
            distance[node] = 0;
            source[node] = node;
        }
        true
    });
    BreadthFirst {
        order,
        distance,
        parent,
        source,
    }
}

fn bfs_order<I: IntoIterator<Item = NodeId>>(
    sources: &[NodeId],
    neighbors: impl Fn(NodeId) -> I,
    mut discover: impl FnMut(NodeId, Option<NodeId>) -> bool,
) -> Vec<NodeId> {
    let mut order: Vec<_> = sources
        .iter()
        .copied()
        .filter(|&node| discover(node, None))
        .collect();
    let mut head = 0;
    while let Some(&node) = order.get(head) {
        head += 1;
        order.extend(
            neighbors(node)
                .into_iter()
                .filter(|&next| discover(next, Some(node))),
        );
    }
    order
}

/// Returns DFS preorder from one vertex, using O(V) auxiliary space.
pub fn dfs<W>(graph: &Graph<W>, start: NodeId) -> Result<Vec<NodeId>, GraphError> {
    graph.validate_node(start)?;
    let mut visited = vec![false; graph.node_count()];
    let mut stack = Vec::new();
    let mut order = Vec::new();
    depth_first(
        graph,
        start,
        false,
        &mut visited,
        &mut stack,
        |node, exiting| {
            if !exiting {
                order.push(node);
            }
        },
    );
    Ok(order)
}

pub(crate) fn depth_first<W>(
    graph: &Graph<W>,
    start: NodeId,
    reverse: bool,
    visited: &mut [bool],
    stack: &mut Vec<(NodeId, usize)>,
    mut visit: impl FnMut(NodeId, bool),
) {
    visited[start] = true;
    visit(start, false);
    stack.push((start, 0));
    while let Some((node, cursor)) = stack.last_mut() {
        let arcs = if reverse {
            graph.inc(*node)
        } else {
            graph.out(*node)
        };
        if let Some(arc) = arcs.get(*cursor) {
            *cursor += 1;
            if !visited[arc.node] {
                visited[arc.node] = true;
                visit(arc.node, false);
                stack.push((arc.node, 0));
            }
        } else {
            visit(*node, true);
            stack.pop();
        }
    }
}
