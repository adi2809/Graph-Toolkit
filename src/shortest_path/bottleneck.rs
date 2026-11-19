//! Minimax paths that minimize the largest edge weight along a route.

use crate::queue::IndexedHeap;
use crate::{EdgeId, Graph, GraphError, NodeId};

use super::{NO_LINK, trace_path};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BottleneckPath<W> {
    /// The largest edge weight on the path. A zero-edge path has no bottleneck.
    pub bottleneck: Option<W>,
    pub nodes: Vec<NodeId>,
    pub edges: Vec<EdgeId>,
}

/// Finds a path minimizing the largest edge weight in O((V + E) log V) time.
/// Ordered signed weights work as well as nonnegative weights.
pub fn minimum_bottleneck_path<W: Copy + Ord>(
    graph: &Graph<W>,
    source: NodeId,
    target: NodeId,
) -> Result<Option<BottleneckPath<W>>, GraphError> {
    graph.validate_node(source)?;
    graph.validate_node(target)?;
    if source == target {
        return Ok(Some(BottleneckPath {
            bottleneck: None,
            nodes: vec![source],
            edges: Vec::new(),
        }));
    }
    let (cost, predecessor) = bottleneck_search(graph, source, Some(target));
    Ok(cost[target].map(|bottleneck| {
        let (nodes, edges) = trace_path(&predecessor, target);
        BottleneckPath {
            bottleneck: Some(bottleneck),
            nodes,
            edges,
        }
    }))
}

pub(crate) fn bottleneck_search<W: Copy + Ord>(
    graph: &Graph<W>,
    source: NodeId,
    target: Option<NodeId>,
) -> (Vec<Option<W>>, Vec<(NodeId, EdgeId)>) {
    let n = graph.node_count();
    let mut labels = vec![None; n];
    let mut predecessor = vec![NO_LINK; n];
    let mut queue = IndexedHeap::new(n);
    queue.decrease(source, None);
    while let Some((cost, node)) = queue.pop() {
        if Some(node) == target {
            break;
        }
        for arc in graph.out(node) {
            if arc.node == source {
                continue;
            }
            let weight = graph.edges()[arc.edge].weight;
            let candidate = cost.map_or(weight, |cost: W| cost.max(weight));
            if labels[arc.node].is_none_or(|old| candidate < old) {
                labels[arc.node] = Some(candidate);
                predecessor[arc.node] = (node, arc.edge);
                queue.decrease(arc.node, Some(candidate));
            }
        }
    }
    (labels, predecessor)
}
