#[path = "support/flow_certificate.rs"]
mod certificate;
use graph_library::{
    Graph, GraphError, GraphKind,
    flow::{
        CycleRatioWeight, EmbeddedCycleHierarchy, max_flow_core_cycles, recover_dyadic_max_flow,
    },
};

#[test]
fn compare_to_every_cut_on_weighted_directed_multigraphs() {
    let mut seed = 932852321u64;
    let mut rng = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for case in 0..1200 {
        let n = 2 + case % 7;
        let s = case % n;
        let t = (s + 1 + case % (n - 1)) % n;
        let mut b = Graph::builder(n, GraphKind::Directed);
        for _ in 0..rng() % 45 {
            b.add_edge(rng() as usize % n, rng() as usize % n, rng() % 43)
                .unwrap();
        }
        let g = b.build();
        let opt = (0..1usize << n)
            .filter(|a| a >> s & 1 == 1 && a >> t & 1 == 0)
            .map(|a| {
                g.edges()
                    .iter()
                    .filter(|e| a >> e.from & 1 == 1 && a >> e.to & 1 == 0)
                    .map(|e| u128::from(e.weight))
                    .sum::<u128>()
            })
            .min()
            .unwrap();
        let answer = max_flow_core_cycles(&g, s, t).unwrap();
        assert_eq!(answer.flow.value, opt);
        certificate::verify(&g, s, t, &answer.flow);
        if opt != 0 {
            assert!(answer.fractional_edges > 0 && answer.epochs > 0);
        }
    }
}

#[test]
fn recovery_rounds_up_and_checks_the_residual_certificate() {
    let mut b = Graph::builder(5, GraphKind::Directed);
    // Two 3/8 branches merge into a 3/4 unit bottleneck. A separate
    // fractional circulation includes antiparallel edges and a loop.
    for (a, z) in [
        (0, 1),
        (1, 3),
        (0, 2),
        (2, 3),
        (3, 4),
        (1, 2),
        (2, 1),
        (1, 1),
    ] {
        b.add_edge(a, z, 1).unwrap();
    }
    let g = b.build();
    let f = recover_dyadic_max_flow(&g, 0, 4, &[3, 3, 3, 3, 6, 7, 7, 5], 3)
        .unwrap()
        .unwrap();
    certificate::verify(&g, 0, 4, &f);
    assert_eq!(f.value, 1);
    assert!(
        recover_dyadic_max_flow(&g, 0, 4, &[0; 8], 3)
            .unwrap()
            .is_none()
    );
    assert_eq!(
        recover_dyadic_max_flow(&g, 0, 4, &[0, 1, 0, 0, 0, 0, 0, 0], 3),
        Err(GraphError::InvalidModelParameters)
    );
    assert!(recover_dyadic_max_flow(&g, 0, 4, &[9; 8], 3).is_err());
    assert!(recover_dyadic_max_flow(&g, 0, 4, &[], 3).is_err());
    assert!(recover_dyadic_max_flow(&g, 0, 4, &[0; 8], 65).is_err());
    let mut b = Graph::builder(2, GraphKind::Directed);
    for _ in 0..4 {
        b.add_edge(0, 1, u64::MAX).unwrap();
    }
    let g = b.build();
    let x = (u128::from(u64::MAX) << 64) - 1;
    let f = recover_dyadic_max_flow(&g, 0, 1, &[x; 4], 64)
        .unwrap()
        .unwrap();
    assert_eq!(f.value, 4 * u128::from(u64::MAX));
    certificate::verify(&g, 0, 1, &f);
    let independent = max_flow_core_cycles(&g, 0, 1).unwrap();
    assert_eq!(independent.flow.value, f.value);
    certificate::verify(&g, 0, 1, &independent.flow);
}

#[test]
fn dynamic_walks_have_current_gradients_lengths_and_original_edges() {
    let n = 37;
    let mut r = 701u64;
    let mut next = || {
        r ^= r << 13;
        r ^= r >> 7;
        r ^= r << 17;
        r
    };
    let mut b = Graph::builder(n, GraphKind::Undirected);
    for i in 1..n {
        b.add_edge(i - 1, i, CycleRatioWeight { cost: 0, length: 1 })
            .unwrap();
    }
    for _ in 0..120 {
        b.add_edge(
            next() as usize % n,
            next() as usize % n,
            CycleRatioWeight {
                cost: (next() % 21) as i64 - 10,
                length: 1 + next() % 15,
            },
        )
        .unwrap();
    }
    let g = b.build();
    let mut weights: Vec<_> = g.edges().iter().map(|e| e.weight).collect();
    let mut hierarchy = EmbeddedCycleHierarchy::new(&g).unwrap();
    for _ in 0..300 {
        let changes: Vec<_> = (0..4)
            .map(|_| {
                (
                    next() as usize % weights.len(),
                    CycleRatioWeight {
                        cost: (next() % 101) as i64 - 50,
                        length: 1 + next() % 50,
                    },
                )
            })
            .collect();
        for &(id, w) in &changes {
            weights[id] = w;
        }
        hierarchy.update(&changes).unwrap();
        let c = hierarchy.improving_walk().unwrap();
        let mut cost = 0;
        let mut len = 0;
        let mut ends = Vec::new();
        for a in c.steps {
            let e = &g.edges()[a.edge];
            let w = weights[a.edge];
            len += u128::from(w.length);
            if a.reversed {
                cost -= i128::from(w.cost);
                ends.push((e.to, e.from));
            } else {
                cost += i128::from(w.cost);
                ends.push((e.from, e.to));
            }
        }
        for i in 0..ends.len() {
            assert_eq!(ends[i].1, ends[(i + 1) % ends.len()].0);
        }
        assert_eq!((cost, len), (c.total_cost, c.total_length));
        assert!(cost < 0);
    }
    let zero: Vec<_> = weights
        .iter()
        .enumerate()
        .map(|(i, w)| (i, CycleRatioWeight { cost: 0, ..*w }))
        .collect();
    hierarchy.update(&zero).unwrap();
    assert!(hierarchy.improving_walk().is_none());
}

fn check_simple_cycle(g: &Graph<CycleRatioWeight>, cycle: &graph_library::flow::RatioCycle) {
    let mut vertices = vec![false; g.node_count()];
    let mut cost = 0;
    let mut length = 0;
    let mut endpoints = Vec::new();
    for step in &cycle.steps {
        let e = &g.edges()[step.edge];
        let (u, v, sign) = if step.reversed {
            (e.to, e.from, -1)
        } else {
            (e.from, e.to, 1)
        };
        assert!(!std::mem::replace(&mut vertices[u], true));
        endpoints.push((u, v));
        cost += sign * i128::from(e.weight.cost);
        length += u128::from(e.weight.length);
    }
    assert!(!endpoints.is_empty());
    for (i, &(_, v)) in endpoints.iter().enumerate() {
        assert_eq!(v, endpoints[(i + 1) % endpoints.len()].0);
    }
    assert_eq!((cost, length), (cycle.total_cost, cycle.total_length));
    assert!(cost < 0);
}

// Enumerate all vertex-simple closed walks, independent of the production
// contraction hierarchy, its dual bounds, and its threshold decision routine.
fn brute_ratio(g: &Graph<CycleRatioWeight>) -> (i128, i128) {
    let mut best = (0, 1);
    for root in 0..g.node_count() {
        let mut stack = vec![(root, 1u64 << root, 0i128, 0i128)];
        while let Some((u, used, cost, length)) = stack.pop() {
            for a in g.outgoing(u).unwrap() {
                let e = &g.edges()[a.edge];
                let c = i128::from(e.weight.cost);
                let mut sum = cost + if e.from == u { c } else { -c };
                let len = length + i128::from(e.weight.length);
                if a.node == u {
                    sum = cost - c.abs();
                }
                if a.node == root {
                    if sum * best.1 < best.0 * len {
                        best = (sum, len);
                    }
                } else if a.node > root && used >> a.node & 1 == 0 {
                    stack.push((a.node, used | (1 << a.node), sum, len));
                }
            }
        }
    }
    best
}

#[test]
fn signed_hierarchy_certifies_global_ratios_after_edits() {
    let mut edges = vec![];
    for u in 0..6 {
        for v in u..6 {
            edges.push((
                u,
                v,
                CycleRatioWeight {
                    cost: 3 * (v as i64 - u as i64),
                    length: 4,
                },
            ));
        }
    }
    edges.push((2, 3, CycleRatioWeight { cost: 3, length: 4 }));
    let build = |edges: &[(usize, usize, CycleRatioWeight)]| {
        let mut b = Graph::builder(6, GraphKind::Undirected);
        for &(u, v, w) in edges {
            b.add_edge(u, v, w).unwrap();
        }
        b.build()
    };
    let mut hierarchy = EmbeddedCycleHierarchy::new(&build(&edges)).unwrap();
    assert!(hierarchy.spanner_replacements() > 0);
    assert!(hierarchy.certified_cycle(1).unwrap().is_none());
    assert!(hierarchy.certified_cycle(0).is_err());
    let mut seed = 134764u64;
    for round in 0..200 {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        let i = seed as usize % edges.len();
        let w = CycleRatioWeight {
            cost: ((seed >> 9) % 77) as i64 - 38,
            length: if round % 13 == 0 { 47 } else { 3 + seed % 4 },
        };
        hierarchy.update(&[(i, edges[i].2), (i, w)]).unwrap();
        edges[i].2 = w;
        let g = build(&edges);
        let (num, denom) = brute_ratio(&g);
        for factor in [1, 2, 7] {
            let result = hierarchy.certified_cycle(factor).unwrap();
            if num == 0 {
                assert!(result.is_none());
                continue;
            }
            let c = result.unwrap();
            check_simple_cycle(&g, &c);
            let length = c.total_length as i128;
            assert!(num * length <= c.total_cost * denom);
            assert!(i128::from(factor) * c.total_cost * denom <= num * length);
        }
        let previous = hierarchy.improving_walk();
        let invalid = CycleRatioWeight {
            cost: i64::MIN,
            length: 0,
        };
        assert!(hierarchy.update(&[(i, invalid), (edges.len(), w)]).is_err());
        assert_eq!(previous, hierarchy.improving_walk());
    }
}

#[test]
fn flow_handles_a_long_route_and_an_interior_cut_above_u64() {
    let n = 10_003;
    let mut b = Graph::builder(n, GraphKind::Directed);
    for u in 0..n - 1 {
        let cap = if u == n / 2 {
            u64::MAX / 2 + 2
        } else {
            u64::MAX
        };
        b.add_edge(u, u + 1, cap).unwrap();
        b.add_edge(u, u + 1, cap).unwrap();
        b.add_edge(u + 1, u, 9).unwrap();
    }
    let g = b.build();
    let answer = max_flow_core_cycles(&g, 0, n - 1).unwrap();
    assert_eq!(answer.flow.value, 2 * u128::from(u64::MAX / 2 + 2));
    certificate::verify(&g, 0, n - 1, &answer.flow);
}

#[test]
fn shared_long_paths_and_dynamic_barrier_keep_original_edge_flows_feasible() {
    use graph_library::flow::{TreeCycleBasis, max_flow_tree_cycle};
    let layers = 66;
    let width = 9;
    let source = layers * width;
    let sink = source + 1;
    let mut b = Graph::builder(sink + 1, GraphKind::Directed);
    for i in 0..width {
        b.add_edge(source, i, 77 + i as u64).unwrap();
        b.add_edge((layers - 1) * width + i, sink, 191).unwrap();
    }
    for k in 1..layers {
        for i in 0..width {
            let u = (k - 1) * width + i;
            let v = k * width + i;
            let capacity = if k == layers / 2 { 7 + i as u64 } else { 280 };
            b.add_edge(u, v, capacity).unwrap();
            b.add_edge(u, k * width + (i + 1) % width, capacity)
                .unwrap();
            b.add_edge(v, u, 11).unwrap();
        }
    }
    let g = b.build();
    let expected = 2 * (0..width).map(|i| 7 + i as u128).sum::<u128>();
    let f = max_flow_core_cycles(&g, source, sink).unwrap();
    assert!(f.augmentations >= width);
    assert_eq!(f.flow.value, expected);
    certificate::verify(&g, source, sink, &f.flow);
    let other = max_flow_tree_cycle(&g, source, sink, 5, TreeCycleBasis::DynamicHierarchy).unwrap();
    assert!(other.descent_steps > 0);
    assert_eq!(other.flow.value, expected);
    certificate::verify(&g, source, sink, &other.flow);
}
