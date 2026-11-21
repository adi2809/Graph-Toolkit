use graph_library::{Graph, GraphKind, community::*};
#[test]
fn weighted_connected_partitions_and_quality() {
    let mut b = Graph::builder(13, GraphKind::Undirected);
    for offset in [0, 4, 8] {
        for u in offset..offset + 4 {
            for v in u + 1..offset + 4 {
                b.add_edge(u, v, 2.).unwrap();
            }
        }
    }
    for (u, v, w) in [
        (3, 4, 0.01),
        (7, 8, 0.01),
        (7, 8, 0.01),
        (0, 0, 2.),
        (12, 3, 0.),
    ] {
        b.add_edge(u, v, w).unwrap();
    }
    let g = b.build();
    for seed in 0..80 {
        let opt = LeidenOptions {
            seed,
            ..Default::default()
        };
        let f = leiden_communities(&g, opt).unwrap();
        assert_eq!(f, leiden_communities(&g, opt).unwrap());
        assert!(f.contractions > 0);
        assert_eq!(f.community_of, [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3]);
        assert!((modularity(&g, &f.community_of, 1.).unwrap() - f.modularity).abs() < 1e-12);
    }
}
#[test]
fn random_graphs_preserve_connectivity_and_modularity_after_contraction() {
    let mut r = 321u64;
    let mut next = || {
        r ^= r << 13;
        r ^= r >> 7;
        r ^= r << 17;
        r
    };
    for seed in 0..100 {
        let mut b = Graph::builder(24, GraphKind::Undirected);
        for _ in 0..64 {
            b.add_edge(
                next() as usize % 24,
                next() as usize % 24,
                (next() % 8) as f64,
            )
            .unwrap();
        }
        let g = b.build();
        let f = leiden_communities(
            &g,
            LeidenOptions {
                seed,
                ..Default::default()
            },
        )
        .unwrap();
        assert!((modularity(&g, &f.community_of, 1.).unwrap() - f.modularity).abs() < 1e-12);
        let singletons: Vec<_> = (0..24).collect();
        assert!(f.modularity >= modularity(&g, &singletons, 1.).unwrap() - 1e-12);
        for u in 0..24 {
            let mut seen = [false; 24];
            seen[u] = true;
            let mut stack = vec![u];
            while let Some(v) = stack.pop() {
                for a in g.outgoing(v).unwrap() {
                    if !seen[a.node]
                        && g.edges()[a.edge].weight > 0.
                        && f.community_of[a.node] == f.community_of[u]
                    {
                        seen[a.node] = true;
                        stack.push(a.node);
                    }
                }
            }
            assert!((0..24).all(|v| seen[v] == (f.community_of[v] == f.community_of[u])));
        }
    }
}
