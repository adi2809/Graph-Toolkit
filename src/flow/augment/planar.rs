//! Hassin's cofacial directed flow reduction to a single dual shortest path computation.
use super::{MaxFlow, Network, validate_terminals};
use crate::{
    Graph, GraphError, GraphKind, planarity::check_planarity, queue::RadixHeap,
    shortest_path::dijkstra_core,
};
use std::collections::BTreeMap;

/// Exact maximum flow when the underlying simple graph has a planar embedding
/// with both terminals on one face. Returns None if no such embedding exists.
/// Parallel and antiparallel capacities are aggregated before embedding and
/// flows are distributed back to original edge IDs. O((V+E) log(V+2)) time
/// and O(V+E) space, including the embedding and cut certificate.
pub fn max_flow_st_planar(
    graph: &Graph<u64>,
    source: usize,
    sink: usize,
) -> Result<Option<MaxFlow>, GraphError> {
    validate_terminals(graph, source, sink)?;
    let n = graph.node_count();
    let mut primal = Graph::builder(
        n.checked_add(1).ok_or(GraphError::GraphSizeOverflow)?,
        GraphKind::Undirected,
    );
    let mut ids = BTreeMap::new();
    let mut capacities = Vec::<[u128; 2]>::new();
    let mut group = Vec::with_capacity(graph.edge_count());
    for e in graph.edges() {
        if e.from == e.to {
            group.push(None);
            continue;
        }
        let key = (e.from.min(e.to), e.from.max(e.to));
        let id = *ids.entry(key).or_insert_with(|| {
            capacities.push([0; 2]);
            primal.add_edge(key.0, key.1, ()).unwrap()
        });
        let direction = usize::from(e.from > e.to);
        capacities[id][direction] += u128::from(e.weight);
        group.push(Some((id, direction)));
    }
    // Subdivision preserves an independent return edge even for adjacent terminals.
    let closing = primal.add_edge(sink, n, ())?;
    primal.add_edge(n, source, ())?;
    let primal = primal.build();
    let Some(embedding) = check_planarity(&primal) else {
        return Ok(None);
    };
    let (face_count, faces) = embedding.dart_face_ids();
    let mut sides = vec![[0; 2]; primal.edge_count()];
    let mut darts = faces.into_iter();
    for u in 0..primal.node_count() {
        for a in embedding.clockwise_neighbors(u)? {
            sides[a.edge][usize::from(primal.edges()[a.edge].from != u)] = darts.next().unwrap();
        }
    }
    // Flow is the potential difference of incident faces. Their differences
    // telescope around each primal vertex and enforce conservation.
    let mut dual = Graph::builder(face_count, GraphKind::Directed);
    for (&[forward, backward], &[a, b]) in capacities.iter().zip(&sides) {
        dual.add_edge(b, a, forward)?;
        dual.add_edge(a, b, backward)?;
    }
    let dual = dual.build();
    let [target, root] = sides[closing];
    let mut distance = vec![u128::MAX; face_count];
    dijkstra_core(
        &dual,
        root,
        false,
        |arc| dual.edges()[arc.edge].weight,
        &mut distance,
        &mut RadixHeap::new(face_count),
        |_, _, _| {},
    );
    let mut state = Network::new(graph);
    state.value = distance[target];
    let mut remaining: Vec<_> = sides[..closing]
        .iter()
        .map(|&[a, b]| {
            let (x, y) = (distance[a], distance[b]);
            (usize::from(x < y), x.abs_diff(y))
        })
        .collect();
    for (i, item) in group.into_iter().enumerate() {
        if let Some((id, direction)) = item {
            let (wanted, amount) = &mut remaining[id];
            if direction == *wanted {
                let used = (*amount).min(u128::from(state.capacity[i][0])) as u64;
                state.capacity[i][0] -= used;
                state.capacity[i][1] = used;
                *amount -= u128::from(used);
            }
        }
    }
    Ok(Some(state.finish(source)))
}
