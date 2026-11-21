//! Maximum cardinality search, perfect elimination, and chordal certificates.
use crate::{Graph, GraphError, GraphKind, NONE};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChordalDecomposition {
    /// Each vertex's later neighbors form a clique.
    pub elimination_order: Vec<usize>,
    /// Minimum proper coloring of the simple projection, using IDs from zero.
    pub colors: Vec<usize>,
    /// A clique certifying the number of colors is optimal.
    pub maximum_clique: Vec<usize>,
}

/// Recognizes chordal graphs in O(V+E) time and space. Returns None for an
/// induced cycle of length at least four. Parallel edges act as one edge,
/// and self-loops are ignored. Works on disconnected undirected graphs.
pub fn chordal_decomposition<W>(
    graph: &Graph<W>,
) -> Result<Option<ChordalDecomposition>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    let mut buckets = vec![Vec::new(); n + 1];
    buckets[0].extend((0..n).rev());
    let mut label = vec![0; n];
    let mut rank = vec![NONE; n];
    let mut seen = vec![NONE; n];
    let mut highest = 0;
    let mut order = Vec::with_capacity(n);
    while order.len() < n {
        let Some(u) = buckets[highest].pop() else {
            highest -= 1;
            continue;
        };
        if rank[u] != NONE || label[u] != highest {
            continue;
        }
        rank[u] = n - 1 - order.len();
        order.push(u);
        for a in graph.out(u) {
            let v = a.node;
            if rank[v] == NONE && seen[v] != u {
                seen[v] = u;
                label[v] += 1;
                highest = highest.max(label[v]);
                buckets[label[v]].push(v);
            }
        }
    }
    order.reverse();
    // Verify later-neighbor containment grouped by the earliest later
    // neighbor. Mark each parent's neighborhood only once for all its children.
    let mut children = vec![Vec::new(); n];
    for u in 0..n {
        if let Some(v) = graph
            .out(u)
            .iter()
            .map(|a| a.node)
            .filter(|&v| rank[v] > rank[u])
            .min_by_key(|&v| rank[v])
        {
            children[v].push(u);
        }
    }
    seen.fill(NONE);
    for (parent, group) in children.into_iter().enumerate() {
        for a in graph.out(parent) {
            seen[a.node] = parent;
        }
        for u in group {
            if graph
                .out(u)
                .iter()
                .any(|a| rank[a.node] > rank[u] && a.node != parent && seen[a.node] != parent)
            {
                return Ok(None);
            }
        }
    }
    let mut colors = vec![NONE; n];
    let mut used = vec![NONE; n];
    seen.fill(NONE);
    let mut maximum_clique = Vec::new();
    let mut neighbors = Vec::new();
    for &u in order.iter().rev() {
        neighbors.clear();
        for a in graph.out(u) {
            let v = a.node;
            if colors[v] != NONE && seen[v] != u {
                seen[v] = u;
                used[colors[v]] = u;
                neighbors.push(v);
            }
        }
        let mut color = 0;
        while used[color] == u {
            color += 1;
        }
        colors[u] = color;
        if neighbors.len() + 1 > maximum_clique.len() {
            maximum_clique.clone_from(&neighbors);
            maximum_clique.push(u);
        }
    }
    Ok(Some(ChordalDecomposition {
        elimination_order: order,
        colors,
        maximum_clique,
    }))
}
