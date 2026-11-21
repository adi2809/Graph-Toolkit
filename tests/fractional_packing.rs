use graph_library::{Graph, GraphKind, multicommodity::*};

#[test]
fn primal_dual_bracket_and_round_limit() {
    let mut g = Graph::builder(5, GraphKind::Directed);
    for (u, v, c) in [
        (0, 1, 20),
        (0, 2, 20),
        (1, 3, 4),
        (2, 3, 4),
        (3, 4, 100),
        (0, 4, 0),
        (0, 0, 2),
    ] {
        g.add_edge(u, v, c).unwrap();
    }
    let g = g.build();
    let demands = [
        Commodity {
            source: 0,
            sink: 3,
            demand: 3,
        },
        Commodity {
            source: 0,
            sink: 4,
            demand: 1,
        },
    ];
    let opt = 2.;
    let cut_short = approximate_concurrent_flow(&g, &demands, 0.05, 1).unwrap();
    assert!(!cut_short.reached_tolerance);
    assert!(cut_short.flow.fraction <= opt && cut_short.upper_bound >= opt);
    let r = approximate_concurrent_flow(&g, &demands, 0.05, 10_000).unwrap();
    assert!(r.reached_tolerance);
    assert!(
        r.flow.fraction >= 0.95 * r.upper_bound
            && r.flow.fraction <= opt + 1e-10
            && r.upper_bound >= opt - 1e-10
    );
    let zero = approximate_concurrent_flow(
        &g,
        &[Commodity {
            source: 3,
            sink: 0,
            demand: 2,
        }],
        0.1,
        10,
    )
    .unwrap();
    assert_eq!((zero.flow.fraction, zero.upper_bound), (0., 0.));
    assert!(approximate_concurrent_flow(&g, &demands, 0., 10).is_err());
}
#[cfg(feature = "multicommodity")]
#[test]
fn generated_instances_enclose_lp_optimum() {
    let mut rng = 4301u64;
    let mut next = || {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        rng
    };
    for _ in 0..25 {
        let mut g = Graph::builder(7, GraphKind::Directed);
        for u in 0..7 {
            g.add_edge(u, (u + 1) % 7, 2 + next() % 6).unwrap();
        }
        for _ in 0..20 {
            g.add_edge(next() as usize % 7, next() as usize % 7, next() % 6)
                .unwrap();
        }
        let g = g.build();
        let demands = [
            Commodity {
                source: 0,
                sink: 4,
                demand: 2,
            },
            Commodity {
                source: 2,
                sink: 5,
                demand: 1,
            },
            Commodity {
                source: 6,
                sink: 1,
                demand: 3,
            },
        ];
        let optimal = maximum_concurrent_flow(&g, &demands).unwrap().fraction;
        let bounds = approximate_concurrent_flow(&g, &demands, 0.1, 3000).unwrap();
        assert!(bounds.reached_tolerance);
        assert!(bounds.flow.fraction <= optimal + 1e-8 && bounds.upper_bound >= optimal - 1e-8);
    }
}
