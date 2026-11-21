//! Replace a core edge by a short two-edge walk with reserved spokes.
//! Protected forest edges survive every level. Each reserved
//! edge serves at most MAX_LOAD omitted edges at this level.
use super::{Builder, CoreEdge, Forest, Walk};
use crate::{Graph, GraphKind};

const MAX_LOAD: u8 = 8;
const NEIGHBOR_PROBES: usize = 8;

fn oriented(edge: &CoreEdge, from: usize) -> Walk {
    if edge.u == from {
        edge.walk
    } else {
        edge.walk.reverse()
    }
}

fn covers(b: &Builder, shortcut: Walk, path: [Walk; 2]) -> bool {
    for backwards in [false, true] {
        let orient = |w: Walk| if backwards { w.reverse() } else { w };
        let direct = b.arena.metric(orient(shortcut));
        if direct.open() {
            let [a, z] = path.map(|w| b.arena.metric(orient(w)));
            if !a.open() || !z.open() || a.length + z.length > direct.length.saturating_mul(4) {
                return false;
            }
        }
    }
    true
}

pub(super) fn two_hop(b: &mut Builder, edges: Vec<CoreEdge>, forest: &Forest) -> Vec<CoreEdge> {
    if edges.len() < 3 {
        return edges;
    }
    let mut graph = Graph::builder(forest.parent.len(), GraphKind::Undirected);
    for e in &edges {
        graph.add_edge(e.u, e.v, ()).unwrap();
    }
    let graph = graph.build();
    // Edges arrive in lexicographic order. Each adjacency list therefore
    // contains first its lower neighbors and then its upper neighbors, sorted.
    let mut load = vec![0u8; edges.len()];
    let mut removed = vec![false; edges.len()];
    for (i, e) in edges.iter().enumerate() {
        if load[i] > 0 || forest.parent[e.u] == e.v || forest.parent[e.v] == e.u {
            continue;
        }
        let (start, end) = if graph.out(e.u).len() <= graph.out(e.v).len() {
            (e.u, e.v)
        } else {
            (e.v, e.u)
        };
        let outgoing = graph.out(start);
        let incoming = graph.out(end);
        for a in outgoing.iter().take(NEIGHBOR_PROBES) {
            if a.node == end || removed[a.edge] || load[a.edge] == MAX_LOAD {
                continue;
            }
            let Ok(pos) = incoming.binary_search_by_key(&a.node, |z| z.node) else {
                continue;
            };
            let second = incoming[pos].edge;
            if removed[second] || load[second] == MAX_LOAD {
                continue;
            }
            let path = [
                oriented(&edges[a.edge], start),
                oriented(&edges[second], a.node),
            ];
            let direct = oriented(e, start);
            if !covers(b, direct, path) {
                continue;
            }
            if b.record_cycles {
                let through = b.arena.join(path[0], path[1]);
                let cycle = b.arena.join(direct, through.reverse());
                b.close(cycle);
            }
            load[a.edge] += 1;
            load[second] += 1;
            removed[i] = true;
            b.shortcuts += 1;
            break;
        }
    }
    edges
        .into_iter()
        .zip(removed)
        .filter_map(|(e, dropped)| (!dropped).then_some(e))
        .collect()
}
