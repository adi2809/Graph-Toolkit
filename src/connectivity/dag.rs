use std::cmp::Reverse;
use std::collections::BinaryHeap;

use super::topological_sort;
use crate::{Graph, GraphBuilder, GraphError, GraphKind, NONE, NodeId};

/// Reverse DFS finish order as in Boost.Graph's output iterator convention
/// (sinks first). Uses an explicit stack and detects back edges in O(V + E).
pub fn reverse_topological_sort<W>(graph: &Graph<W>) -> Result<Vec<NodeId>, GraphError> {
    graph.require(GraphKind::Directed)?;
    let mut color = vec![0; graph.node_count()];
    let mut order = Vec::with_capacity(color.len());
    let mut stack = Vec::new();
    for root in 0..color.len() {
        if color[root] != 0 {
            continue;
        }
        color[root] = 1;
        stack.push((root, 0));
        while let Some((node, cursor)) = stack.last_mut() {
            if let Some(arc) = graph.out(*node).get(*cursor) {
                *cursor += 1;
                match color[arc.node] {
                    1 => return Err(GraphError::Cycle),
                    0 => {
                        color[arc.node] = 1;
                        stack.push((arc.node, 0));
                    }
                    _ => {}
                }
            } else {
                color[*node] = 2;
                order.push(*node);
                stack.pop();
            }
        }
    }
    Ok(order)
}

pub fn is_directed_acyclic_graph<W>(graph: &Graph<W>) -> bool {
    topological_sort(graph).is_ok()
}

/// Kahn traversal with the supplied ordering key and vertex-ID tie breaking.
pub fn lexicographical_topological_sort<W, K: Ord>(
    graph: &Graph<W>,
    mut key: impl FnMut(NodeId) -> K,
) -> Result<Vec<NodeId>, GraphError> {
    graph.require(GraphKind::Directed)?;
    let n = graph.node_count();
    let mut incoming: Vec<_> = (0..n).map(|node| graph.inc(node).len()).collect();
    let mut ready: BinaryHeap<_> = incoming
        .iter()
        .enumerate()
        .filter(|&(_, &degree)| degree == 0)
        .map(|(node, _)| Reverse((key(node), node)))
        .collect();
    let mut order = Vec::with_capacity(n);
    while let Some(Reverse((_, node))) = ready.pop() {
        order.push(node);
        for arc in graph.out(node) {
            incoming[arc.node] -= 1;
            if incoming[arc.node] == 0 {
                ready.push(Reverse((key(arc.node), arc.node)));
            }
        }
    }
    if order.len() != n {
        return Err(GraphError::Cycle);
    }
    Ok(order)
}

/// Earliest possible topological layer per vertex, grouped in increasing order.
pub fn topological_generations<W>(graph: &Graph<W>) -> Result<Vec<Vec<NodeId>>, GraphError> {
    let order = topological_sort(graph)?;
    let mut rank = vec![0; order.len()];
    let mut layers: Vec<Vec<_>> = Vec::new();
    for node in order {
        if rank[node] >= layers.len() {
            layers.resize_with(rank[node] + 1, Vec::new);
        }
        layers[rank[node]].push(node);
        for arc in graph.out(node) {
            rank[arc.node] = rank[arc.node].max(rank[node] + 1);
        }
    }
    Ok(layers)
}

/// Reachability closure of a DAG. Uses word-parallel propagation with
/// O(V²/word) working bits, followed by O(output edges) edge emission.
pub fn transitive_closure_dag<W>(graph: &Graph<W>) -> Result<Graph<()>, GraphError> {
    let order = topological_sort(graph)?;
    let n = order.len();
    let words = n.div_ceil(64);
    let cells = n.checked_mul(words).ok_or(GraphError::GraphSizeOverflow)?;
    let mut reachable = Vec::new();
    reachable
        .try_reserve_exact(cells)
        .map_err(|_| GraphError::GraphSizeOverflow)?;
    reachable.resize(cells, 0_u64);
    for &from in order.iter().rev() {
        for arc in graph.out(from) {
            let to = arc.node;
            if reachable[from * words + to / 64] & (1_u64 << (to % 64)) != 0 {
                continue;
            }
            reachable[from * words + to / 64] |= 1_u64 << (to % 64);
            for word in 0..words {
                reachable[from * words + word] |= reachable[to * words + word];
            }
        }
    }
    let mut builder = GraphBuilder::new(n, GraphKind::Directed);
    for (from, row) in reachable.chunks(words.max(1)).enumerate() {
        for (word, mut targets) in row.iter().copied().enumerate() {
            while targets != 0 {
                let bit = targets.trailing_zeros() as usize;
                targets &= targets - 1;
                builder.add_edge(from, word * 64 + bit, ()).unwrap();
            }
        }
    }
    Ok(builder.build())
}

/// DAG transitive reduction. Returns original edge IDs, retaining the first
/// inserted edge when parallel arcs have the same endpoints. Uses O(V + E)
/// space and lazily marks descendants of retained neighbors.
pub fn transitive_reduction<W>(graph: &Graph<W>) -> Result<Vec<usize>, GraphError> {
    let order = topological_sort(graph)?;
    let n = order.len();
    let mut rank = vec![0; n];
    for (index, node) in order.into_iter().enumerate() {
        rank[node] = index;
    }
    let mut stamp = vec![NONE; n];
    let mut stack = Vec::new();
    let mut targets = Vec::new();
    let mut edges = Vec::new();
    for from in 0..n {
        targets.clear();
        targets.extend_from_slice(graph.out(from));
        targets.sort_unstable_by_key(|arc| (rank[arc.node], arc.edge));
        let Some(last) = targets.last() else { continue };
        let last = rank[last.node];
        for arc in &targets {
            if stamp[arc.node] == from {
                continue;
            }
            edges.push(arc.edge);
            stamp[arc.node] = from;
            if rank[arc.node] < last {
                stack.push(arc.node);
            }
            while let Some(node) = stack.pop() {
                for next in graph.out(node) {
                    let to = next.node;
                    if rank[to] <= last && stamp[to] != from {
                        stamp[to] = from;
                        if rank[to] < last {
                            stack.push(to);
                        }
                    }
                }
            }
        }
    }
    edges.sort_unstable();
    Ok(edges)
}
