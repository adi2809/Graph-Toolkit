//! General bottleneck matching via monotone threshold queries.

use crate::{EdgeId, Graph, GraphError, GraphKind, NONE};

use super::{GeneralMatching, VACANT, maximum_general_matching_gabow};

#[derive(Clone, Debug)]
pub struct BottleneckMatching<W> {
    pub matching: GeneralMatching,
    /// The largest edge weight, or `None` for an empty matching.
    pub bottleneck: Option<W>,
}

/// Among maximum-cardinality matchings, minimizes the largest selected edge weight.
pub fn minimum_bottleneck_matching<W: Copy + Ord>(
    graph: &Graph<W>,
) -> Result<BottleneckMatching<W>, GraphError> {
    let maximum = maximum_general_matching_gabow(graph)?;
    let target = maximum.len();
    if target == 0 {
        return Ok(BottleneckMatching {
            matching: maximum,
            bottleneck: None,
        });
    }
    let upper = maximum
        .edges()
        .map(|id| graph.edges()[id].weight)
        .max()
        .unwrap();
    threshold_matching(graph, target, Some(upper))
        .map(|result| result.expect("the upper threshold contains a maximum matching"))
}

/// Minimizes the largest edge weight among matchings of exactly `size` edges.
/// Binary search uses the Gabow matcher as its threshold oracle. Returns
/// `None` if the requested matching size is infeasible.
pub fn minimum_bottleneck_matching_of_size<W: Copy + Ord>(
    graph: &Graph<W>,
    size: usize,
) -> Result<Option<BottleneckMatching<W>>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    threshold_matching(graph, size, None)
}

/// Minimum-bottleneck perfect matching, or `None` if no perfect matching exists.
pub fn minimum_bottleneck_perfect_matching<W: Copy + Ord>(
    graph: &Graph<W>,
) -> Result<Option<BottleneckMatching<W>>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    if graph.node_count() % 2 != 0 {
        return Ok(None);
    }
    threshold_matching(graph, graph.node_count() / 2, None)
}

fn threshold_matching<W: Copy + Ord>(
    graph: &Graph<W>,
    size: usize,
    upper: Option<W>,
) -> Result<Option<BottleneckMatching<W>>, GraphError> {
    let n = graph.node_count();
    if size > n / 2 {
        return Ok(None);
    }
    if size == 0 {
        return Ok(Some(BottleneckMatching {
            matching: GeneralMatching {
                mates: vec![VACANT; n],
                len: 0,
            },
            bottleneck: None,
        }));
    }

    let mut compact = vec![NONE; n];
    let mut count = 0;
    let mut candidates: Vec<_> = graph
        .edges()
        .iter()
        .enumerate()
        .filter_map(|(id, edge)| {
            if edge.from == edge.to || upper.is_some_and(|bound| edge.weight > bound) {
                return None;
            }
            for node in [edge.from, edge.to] {
                if compact[node] == NONE {
                    compact[node] = count;
                    count += 1;
                }
            }
            Some(id)
        })
        .collect();
    if size > count / 2 {
        return Ok(None);
    }
    candidates.sort_unstable_by_key(|&id| (graph.edges()[id].weight, id));
    let ends: Vec<_> = candidates
        .chunk_by(|&a, &b| graph.edges()[a].weight == graph.edges()[b].weight)
        .scan(0, |end, group| {
            *end += group.len();
            Some(*end)
        })
        .collect();

    let mut low = 0;
    let mut high = ends.len();
    let mut best = None;
    while low < high {
        let middle = low + (high - low) / 2;
        let prefix = &candidates[..ends[middle]];
        let edges = prefix_matching(graph, &compact, count, prefix)?;
        if edges.len() >= size {
            best = Some(edges);
            high = middle;
        } else {
            low = middle + 1;
        }
    }
    Ok(best.map(|edges| {
        let mut matching = GeneralMatching {
            mates: vec![VACANT; n],
            len: size,
        };
        let mut bottleneck = None;
        for id in edges.into_iter().take(size) {
            let edge = graph.edges()[id];
            matching.pair(edge.from, edge.to, id);
            bottleneck = Some(bottleneck.map_or(edge.weight, |old: W| old.max(edge.weight)));
        }
        BottleneckMatching {
            matching,
            bottleneck,
        }
    }))
}

fn prefix_matching<W>(
    graph: &Graph<W>,
    compact: &[usize],
    count: usize,
    edges: &[EdgeId],
) -> Result<Vec<EdgeId>, GraphError> {
    let mut builder = Graph::builder(count, GraphKind::Undirected);
    for &id in edges {
        let edge = &graph.edges()[id];
        builder.add_edge(compact[edge.from], compact[edge.to], id)?;
    }
    let threshold = builder.build();
    let matching = maximum_general_matching_gabow(&threshold)?;
    Ok(matching
        .edges()
        .map(|id| threshold.edges()[id].weight)
        .collect())
}
