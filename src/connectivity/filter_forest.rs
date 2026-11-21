//! Karger–Klein–Tarjan sampling and heavy-edge filtering. The path-maximum
//! oracle uses a sorted reconstruction tree, costing O(n log n) for a forest.
use super::{DisjointSet, SpanningForest, ThresholdConnectivityIndex};
use crate::{EdgeId, Graph, GraphError, GraphKind, NONE, random::Random};

#[derive(Copy, Clone)]
struct Candidate {
    endpoints: [usize; 2],
    key: (u64, EdgeId),
}

/// Randomized exact minimum spanning forest. Each round contracts two
/// Borůvka phases, recursively solves a half sample, and discards edges
/// heavier than their sampled tree path. Breaking ties by original edge ID
/// gives the same forest for every seed. Returns edges in Kruskal order.
pub fn minimum_spanning_forest_filtered(
    graph: &Graph<u64>,
    seed: u64,
) -> Result<SpanningForest, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let candidates: Vec<_> = graph
        .edges()
        .iter()
        .enumerate()
        .filter(|(_, e)| e.from != e.to)
        .map(|(id, e)| Candidate {
            endpoints: [e.from, e.to],
            key: (e.weight, id),
        })
        .collect();
    let mut edges: Vec<_> = select_forest(graph.node_count(), &candidates, &mut Random(seed))
        .into_iter()
        .map(|i| candidates[i].key.1)
        .collect();
    edges.sort_unstable_by_key(|&id| (graph.edges()[id].weight, id));
    Ok(SpanningForest {
        total_weight: edges
            .iter()
            .map(|&id| u128::from(graph.edges()[id].weight))
            .sum(),
        component_count: graph.node_count() - edges.len(),
        edges,
    })
}

fn select_forest(n: usize, edges: &[Candidate], rng: &mut Random) -> Vec<usize> {
    let mut selected = Vec::new();
    if edges.is_empty() {
        return selected;
    }
    let mut components = DisjointSet::new(n);
    for _ in 0..2 {
        let mut cheapest = vec![NONE; n];
        for (i, e) in edges.iter().enumerate() {
            let [a, b] = e.endpoints.map(|v| components.find(v));
            if a == b {
                continue;
            }
            for c in [a, b] {
                let old = cheapest[c];
                if old == NONE || edges[old].key > e.key {
                    cheapest[c] = i;
                }
            }
        }
        let before = selected.len();
        for i in cheapest.into_iter().filter(|&i| i != NONE) {
            let [a, b] = edges[i].endpoints;
            if components.union(a, b).is_some() {
                selected.push(i);
            }
        }
        if before == selected.len() {
            return selected;
        }
    }

    let mut local = vec![NONE; n];
    let mut count = 0;
    let mut sample = Vec::new();
    let mut contracted = Vec::new();
    for (i, e) in edges.iter().enumerate() {
        let [a, b] = e.endpoints.map(|v| components.find(v));
        if a == b {
            continue;
        }
        let endpoints = [a, b].map(|c| {
            if local[c] == NONE {
                local[c] = count;
                count += 1;
            }
            local[c]
        });
        let arc = Candidate {
            endpoints,
            key: e.key,
        };
        if rng.next() & 1 == 0 {
            sample.push(arc);
        }
        contracted.push((i, arc));
    }
    if contracted.is_empty() {
        return selected;
    }

    let mut forest = Graph::builder(count, GraphKind::Undirected);
    for i in select_forest(count, &sample, rng) {
        let e = sample[i];
        let [a, b] = e.endpoints;
        forest.add_edge(a, b, e.key).unwrap();
    }
    let index = ThresholdConnectivityIndex::new(&forest.build()).unwrap();
    let (original, light): (Vec<_>, Vec<_>) = contracted
        .into_iter()
        .filter(|(_, e)| {
            let [a, b] = e.endpoints;
            index
                .merge_weight(a, b)
                .is_none_or(|largest| largest >= e.key)
        })
        .unzip();
    selected.extend(
        select_forest(count, &light, rng)
            .into_iter()
            .map(|i| original[i]),
    );
    selected
}
