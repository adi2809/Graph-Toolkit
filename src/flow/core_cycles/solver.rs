use super::{Chain, Forest, Metric};
use crate::{
    Graph, GraphError, NONE,
    flow::{MaxFlow, augment::recovery::round_and_certify, residual_neighbors, validate_terminals},
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoreCycleFlow {
    pub flow: MaxFlow,
    /// Number of distinct terminal forests that routed positive flow.
    pub epochs: usize,
    /// Number of positive-throughput lifted cycles.
    pub augmentations: usize,
    /// Original arcs carrying nonintegral flow before exact recovery.
    pub fractional_edges: usize,
}

pub(super) fn residual_metric(cap: [u128; 2], threshold: u128) -> [Metric; 2] {
    let total = cap[0] + cap[1];
    cap.map(|room| {
        if room < threshold {
            super::embedding::CLOSED
        } else {
            Metric {
                gradient: 0,
                length: u128::from(1 + total.ilog2() - room.ilog2()),
            }
        }
    })
}

struct Routes {
    chain: Chain,
    paths: Option<super::paths::PathFlow>,
}

enum Search {
    Routes(Box<Routes>),
    Cut(u128),
}

fn epoch(
    graph: &Graph<u64>,
    cap: &[[u128; 2]],
    s: usize,
    t: usize,
    threshold: u128,
) -> Result<Search, GraphError> {
    let n = graph.node_count();
    // Unit-hop distances provide a dual lower bound for every admissible
    // s-t route. Pushing only on tight arcs leaves this labeling feasible
    // when reverse residual arcs open. Every route lifted in this epoch
    // is therefore globally shortest by hop count at this threshold.
    let mut level = vec![NONE; n];
    level[s] = 0;
    let mut queue = vec![s];
    let mut head = 0;
    while let Some(&u) = queue.get(head) {
        head += 1;
        if level[u] >= level[t] {
            break;
        }
        for a in residual_neighbors(graph, u) {
            if level[a.to] == NONE && cap[a.edge][a.reversed] >= threshold {
                level[a.to] = level[u] + 1;
                queue.push(a.to);
            }
        }
    }
    if level[t] == NONE {
        let cut = graph
            .edges()
            .iter()
            .filter(|e| level[e.from] != NONE && level[e.to] == NONE)
            .map(|e| u128::from(e.weight))
            .sum();
        return Ok(Search::Cut(cut));
    }
    // Split at a distance layer. On the source side choose any tight
    // predecessor. On the sink side use reverse BFS to keep only vertices
    // with a tight suffix to t. Thus any surviving bridge can carry flow.
    let split = (level[t] - 1) / 2;
    let mut parent = vec![NONE; n];
    let mut depth = vec![0; n];
    let mut incoming = vec![NONE; n];
    parent[s] = s;
    parent[t] = t;
    for &u in &queue {
        if level[u] >= split {
            break;
        }
        for a in residual_neighbors(graph, u) {
            if level[a.to] == level[u] + 1
                && parent[a.to] == NONE
                && cap[a.edge][a.reversed] >= threshold
            {
                parent[a.to] = u;
                incoming[a.to] = a.edge;
                depth[a.to] = level[a.to];
            }
        }
    }
    queue.clear();
    queue.push(t);
    head = 0;
    while let Some(&u) = queue.get(head) {
        head += 1;
        for a in residual_neighbors(graph, u) {
            let v = a.to;
            if level[v] != level[u] - 1 || cap[a.edge][a.reversed ^ 1] < threshold {
                continue;
            }
            if level[v] != split && parent[v] == NONE {
                parent[v] = u;
                incoming[v] = a.edge;
                depth[v] = depth[u] + 1;
                queue.push(v);
            }
        }
    }
    // Short routes cost less to enumerate than to split into heavy paths.
    // For long routes keep capacities and notification thresholds implicit.
    let use_paths = level[t] > 2 * (1 + n.ilog2() as usize);
    let path_parent = use_paths.then(|| parent.clone());
    let reached = |u: usize| incoming[u] != NONE || u == s || u == t;
    let input = graph.edges().iter().enumerate().filter_map(|(id, e)| {
        let (u, v) = (e.from, e.to);
        let bridge = reached(u)
            && reached(v)
            && ((level[u] == split && level[v] == split + 1 && cap[id][0] >= threshold)
                || (level[v] == split && level[u] == split + 1 && cap[id][1] >= threshold));
        // A throughput-increasing fundamental cycle consists of two tree
        // paths, a single bridge and the closing arc. Other chords can
        // enter a future basis after the next forest rebuild.
        (incoming[u] == id || incoming[v] == id || bridge).then(|| {
            let mut metric = residual_metric(cap[id], threshold);
            if use_paths {
                metric[usize::from(level[u] <= level[v])] = super::embedding::CLOSED;
            }
            (id, u, v, metric)
        })
    });
    let chain = Chain::build(
        Forest { parent, depth },
        graph.edge_count(),
        input,
        Some((s, t)),
    )?;
    let paths = path_parent
        .map(|p| super::paths::PathFlow::new(graph, p, &incoming, &level, cap, threshold, &chain));
    Ok(Search::Routes(Box::new(Routes { chain, paths })))
}

/// Exact deterministic s–t flow by dynamically maintained core walks and
/// dyadic Euler rounding. A virtual t→s edge has gradient -1. Among the
/// feasible terminal cycles encoded by a forest, minimizing the ratio
/// chooses the least-length retained route. All have minimum hop count in
/// the threshold residual graph, so their length is within the maximum
/// residual length bucket of the globally shortest weighted route (all
/// lengths are at least one). Repair updates both lengths and directed
/// feasibility. A new forest exposes a route whenever one exists at the
/// current scaling threshold. Cut bounds permit recovery as soon as rounded
/// throughput equals the capacity of a known original cut.
///
/// The hierarchy can incur linear recourse on one change, and the number of
/// forest rebuilds is input dependent. Capacity scaling bounds the number
/// of successful augmentations by O(m log(2+mU)) for maximum capacity U.
/// This implementation has no general
/// almost-linear running-time guarantee.
pub fn max_flow_core_cycles(
    graph: &Graph<u64>,
    source: usize,
    sink: usize,
) -> Result<CoreCycleFlow, GraphError> {
    validate_terminals(graph, source, sink)?;
    finish_from_residual(
        graph,
        source,
        sink,
        graph.edges().iter().map(|e| [e.weight, 0]),
        0,
    )
}

pub(in crate::flow) fn finish_from_residual(
    graph: &Graph<u64>,
    source: usize,
    sink: usize,
    residual: impl Iterator<Item = [u64; 2]>,
    initial_value: u128,
) -> Result<CoreCycleFlow, GraphError> {
    let m = graph.edge_count();
    let bits = usize::BITS - m.leading_zeros();
    let unit = 1u128 << bits;
    // Every cut loses at most m/unit < 1. The perturbed optimum therefore
    // rounds upward to the original integral optimum.
    // A warm start stays exactly feasible. Shrink only edges that still
    // have positive forward room in the initial residual network.
    let mut cap: Vec<_> = residual
        .map(|[forward, back]| {
            [
                (u128::from(forward) * unit).saturating_sub(1),
                u128::from(back) * unit,
            ]
        })
        .collect();
    let mut path = Vec::new();
    let mut epochs = 0;
    let mut augmentations = 0;
    // The two singleton terminal cuts give an initial upper bound on the
    // original integral optimum. Keep fractional throughput in two limbs.
    let terminal_cut = |vertex, outgoing| {
        let arcs = if outgoing {
            graph.out(vertex)
        } else {
            graph.inc(vertex)
        };
        arcs.iter()
            .filter(|a| a.node != vertex)
            .map(|a| u128::from(graph.edges()[a.edge].weight))
            .sum::<u128>()
    };
    let mut upper = terminal_cut(source, true).min(terminal_cut(sink, false));
    let mut whole = initial_value;
    let mut remainder = 0u128;
    // At the first threshold at most 2m augmentations fit under the bound.
    // After an unsuccessful search at threshold d, its cut has less than
    // m*d residual room in scaled units. Halving d admits at most 2m more.
    let count = m.max(1) as u128;
    let gap = upper - initial_value;
    let per_edge = (gap / count) * unit + (gap % count) * unit / count;
    let mut threshold = 1u128 << per_edge.max(1).ilog2();
    while whole + u128::from(remainder > 0) < upper {
        let mut routes = match epoch(graph, &cap, source, sink, threshold)? {
            Search::Routes(chain) => chain,
            Search::Cut(cut) => {
                upper = upper.min(cut);
                if whole + u128::from(remainder > 0) == upper {
                    break;
                }
                // At threshold 1 this cut is saturated in the perturbed
                // network and has lost strictly less than one original unit.
                assert!(threshold > 1, "saturated cut must certify upward rounding");
                threshold /= 2;
                continue;
            }
        };
        let chain = &mut routes.chain;
        let best = chain.winner.expect("a meeting yields a terminal core edge");
        debug_assert!(chain.arena.metric(best).open());
        epochs += 1;
        while chain.arena.metric(best).open() {
            let amount = if let Some(paths) = &mut routes.paths {
                paths.push(graph, chain, best, &mut cap)
            } else {
                chain.arena.lift(best, &mut path);
                let amount = path
                    .iter()
                    .map(|&(edge, back)| cap[edge][usize::from(back)])
                    .min()
                    .unwrap();
                debug_assert!(amount >= threshold && path.len() < graph.node_count());
                for &(id, back) in &path {
                    let dir = usize::from(back);
                    cap[id][dir] -= amount;
                    cap[id][dir ^ 1] += amount;
                    chain.arena.update(
                        chain.leaves[id].unwrap(),
                        residual_metric(cap[id], threshold),
                    );
                }
                amount
            };
            augmentations += 1;
            remainder += amount % unit;
            whole += amount / unit + remainder / unit;
            remainder %= unit;
            if whole + u128::from(remainder > 0) == upper {
                break;
            }
            chain.arena.repair();
        }
        if let Some(paths) = routes.paths {
            paths.flush(&mut cap);
        }
    }
    let amounts: Vec<_> = cap.into_iter().map(|pair| pair[1]).collect();
    let fractional_edges = amounts.iter().filter(|x| **x % unit != 0).count();
    let flow = round_and_certify(graph, source, sink, amounts, bits)
        .expect("rounded throughput attains the original cut");
    Ok(CoreCycleFlow {
        flow,
        epochs,
        augmentations,
        fractional_edges,
    })
}
