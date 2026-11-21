#[path = "support/flow_certificate.rs"]
mod certificate;
use graph_library::{Graph, GraphError, GraphKind, flow::*};

#[test]
fn exhaust_small_terminal_cuts_and_recover_original_flows() {
    let mut state = 67u64;
    let mut rng = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for i in 0..700 {
        let n = 2 + i % 7;
        let mut g = Graph::builder(n, GraphKind::Directed);
        for _ in 0..rng() % 30 {
            g.add_edge(rng() as usize % n, rng() as usize % n, rng() % 23)
                .unwrap();
        }
        let g = g.build();
        let s = i % n;
        let t = (s + 1 + i % (n - 1)) % n;
        let min = (0..1usize << n)
            .filter(|mask| mask >> s & 1 != 0 && mask >> t & 1 == 0)
            .map(|mask| {
                g.edges()
                    .iter()
                    .filter(|e| mask >> e.from & 1 != 0 && mask >> e.to & 1 == 0)
                    .map(|e| u128::from(e.weight))
                    .sum::<u128>()
            })
            .min()
            .unwrap();
        let f = max_flow_tree_cycle(&g, s, t, 12, TreeCycleBasis::ResidualAndLowStretch).unwrap();
        assert!(f.value_before_completion <= min);
        assert_eq!(f.flow.value, min);
        certificate::verify(&g, s, t, &f.flow);
        let dynamic = max_flow_dynamic_trees(&g, s, t).unwrap();
        certificate::verify(&g, s, t, &dynamic);
        assert_eq!(dynamic.value, min);
        let planar = max_flow_st_planar(&g, s, t).unwrap();
        if n <= 4 {
            assert!(planar.is_some());
        }
        if let Some(f) = planar {
            assert_eq!(f.value, min);
            certificate::verify(&g, s, t, &f);
        }
    }
}
#[test]
fn planar_obstructions_bridges_wide_flows_and_descent() {
    // Opposite poles of an octahedron cannot share a face.
    let mut g = Graph::builder(6, GraphKind::Directed);
    for v in 0..4 {
        g.add_edge(v, (v + 1) % 4, 1).unwrap();
        g.add_edge(4, v, 2).unwrap();
        g.add_edge(v, 5, 3).unwrap();
    }
    let g = g.build();
    assert!(max_flow_st_planar(&g, 4, 5).unwrap().is_none());
    let f = max_flow_tree_cycle(&g, 4, 5, 1, TreeCycleBasis::ResidualAndLowStretch).unwrap();
    assert_eq!(
        (f.descent_steps, f.value_before_completion, f.flow.value),
        (1, 2, 8)
    );
    assert_eq!(
        max_flow_tree_cycle(&g, 4, 5, 0, TreeCycleBasis::ResidualAndLowStretch)
            .unwrap()
            .value_before_completion,
        0
    );
    let mut g = Graph::builder(5, GraphKind::Directed);
    for (u, v) in [(0, 1), (0, 1), (1, 0), (2, 2), (2, 3), (2, 3)] {
        g.add_edge(u, v, u64::MAX).unwrap();
    }
    let g = g.build();
    for (s, t, want) in [(0, 1, 2 * u128::from(u64::MAX)), (0, 3, 0), (4, 1, 0)] {
        for f in [
            max_flow_st_planar(&g, s, t).unwrap().unwrap(),
            max_flow_tree_cycle(&g, s, t, 20, TreeCycleBasis::ResidualAndLowStretch)
                .unwrap()
                .flow,
        ] {
            assert_eq!(f.value, want);
            certificate::verify(&g, s, t, &f);
        }
    }
    assert_eq!(
        max_flow_st_planar(&g, 1, 1),
        Err(GraphError::IdenticalTerminals)
    );
}
#[test]
fn ratio_cycles_match_exhaustive_simple_cycle_search() {
    // Enumerate all traversals without repeated vertices or logical edges,
    // including self-loops, reversed arcs, and length-two parallel cycles.
    fn visit(
        g: &Graph<CycleRatioWeight>,
        start: usize,
        u: usize,
        seen: usize,
        used: usize,
        total: (i128, i128),
        best: &mut Option<(i128, i128)>,
    ) {
        for a in g.outgoing(u).unwrap() {
            if used >> a.edge & 1 != 0 {
                continue;
            }
            let e = &g.edges()[a.edge];
            let signed = i128::from(e.weight.cost);
            let cost = if e.from == e.to && !g.is_directed() {
                -signed.abs()
            } else if e.from == u {
                signed
            } else {
                -signed
            };
            let val = (total.0 + cost, total.1 + i128::from(e.weight.length));
            if a.node == start {
                if best.is_none_or(|(p, q)| val.0 * q < p * val.1) {
                    *best = Some(val);
                }
            } else if seen >> a.node & 1 == 0 {
                visit(
                    g,
                    start,
                    a.node,
                    seen | 1 << a.node,
                    used | 1 << a.edge,
                    val,
                    best,
                );
            }
        }
    }
    let mut seed = 97531u64;
    let mut rng = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for kind in [GraphKind::Undirected, GraphKind::Directed] {
        for _ in 0..200 {
            let mut b = Graph::builder(4, kind);
            for _ in 0..7 {
                b.add_edge(
                    rng() as usize % 4,
                    rng() as usize % 4,
                    CycleRatioWeight {
                        cost: (rng() % 16) as i64 - 5,
                        length: 1 + rng() % 7,
                    },
                )
                .unwrap();
            }
            let g = b.build();
            let mut best = None;
            for s in 0..4 {
                visit(&g, s, s, 1 << s, 0, (0, 0), &mut best);
            }
            let answer = if g.is_directed() {
                minimum_ratio_cycle(&g)
            } else {
                best = best.filter(|x| x.0 < 0);
                minimum_undirected_ratio_cycle(&g)
            }
            .unwrap();
            assert_eq!(answer.is_some(), best.is_some());
            if let (Some(c), Some((p, q))) = (answer, best) {
                assert_eq!(c.total_cost * q, c.total_length as i128 * p);
                let oriented: Vec<_> = c
                    .steps
                    .iter()
                    .map(|s| {
                        let e = &g.edges()[s.edge];
                        if s.reversed {
                            (
                                e.to,
                                e.from,
                                -i128::from(e.weight.cost),
                                u128::from(e.weight.length),
                            )
                        } else {
                            (
                                e.from,
                                e.to,
                                i128::from(e.weight.cost),
                                u128::from(e.weight.length),
                            )
                        }
                    })
                    .collect();
                for i in 0..oriented.len() {
                    assert_eq!(oriented[i].1, oriented[(i + 1) % oriented.len()].0);
                }
                assert_eq!(c.total_cost, oriented.iter().map(|e| e.2).sum::<i128>());
                assert_eq!(c.total_length, oriented.iter().map(|e| e.3).sum::<u128>());
            }
        }
    }
}

#[test]
fn karp_matches_general_ratio_oracle_on_signed_multigraphs() {
    let mut r = 307u64;
    let mut rng = || {
        r ^= r << 13;
        r ^= r >> 7;
        r ^= r << 17;
        r
    };
    for _ in 0..600 {
        let n = 1 + rng() as usize % 12;
        let mut graph = Graph::builder(n, GraphKind::Directed);
        let mut ratio = Graph::builder(n, GraphKind::Directed);
        for _ in 0..rng() % 40 {
            let u = rng() as usize % n;
            let v = rng() as usize % n;
            let cost = (rng() % 57) as i64 - 18;
            graph.add_edge(u, v, cost).unwrap();
            ratio
                .add_edge(u, v, CycleRatioWeight { cost, length: 1 })
                .unwrap();
        }
        let graph = graph.build();
        let a = minimum_mean_cycle(&graph).unwrap();
        let b = minimum_ratio_cycle(&ratio.build()).unwrap();
        assert_eq!(a.is_some(), b.is_some());
        if let (Some(a), Some(b)) = (a, b) {
            assert_eq!(
                a.total_cost * b.total_length as i128,
                b.total_cost * a.total_length as i128
            );
            for (i, step) in a.steps.iter().enumerate() {
                assert_eq!(
                    graph.edges()[step.edge].to,
                    graph.edges()[a.steps[(i + 1) % a.steps.len()].edge].from
                );
            }
        }
    }
}
