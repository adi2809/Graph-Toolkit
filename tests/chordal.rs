use graph_library::{Graph, GraphKind, connectivity::chordal_decomposition};
#[test]
fn chordality_clique_and_coloring_match_induced_subgraphs() {
    let n = 5;
    for mask in 0u32..1 << 10 {
        let mut g = Graph::builder(n, GraphKind::Undirected);
        let mut adj = [0usize; 5];
        let mut bit = 0;
        for u in 0..n {
            g.add_edge(u, u, ()).unwrap();
            for v in u + 1..n {
                if mask >> bit & 1 != 0 {
                    g.add_edge(u, v, ()).unwrap();
                    g.add_edge(v, u, ()).unwrap();
                    adj[u] |= 1 << v;
                    adj[v] |= 1 << u;
                }
                bit += 1;
            }
        }
        let g = g.build();
        let mut chordal = true;
        let mut omega = 0;
        for set in 1usize..1 << n {
            let k = set.count_ones();
            let degrees = (0..n)
                .filter(|&u| set >> u & 1 != 0)
                .map(|u| (adj[u] & set).count_ones());
            // With 4 or 5 vertices a simple 2-regular graph is a single cycle.
            if k >= 4 && degrees.clone().all(|d| d == 2) {
                chordal = false;
            }
            if degrees.clone().all(|d| d == k - 1) {
                omega = omega.max(k);
            }
        }
        let result = chordal_decomposition(&g).unwrap();
        assert_eq!(result.is_some(), chordal);
        if let Some(r) = result {
            assert_eq!(r.maximum_clique.len(), omega as usize);
            assert_eq!(*r.colors.iter().max().unwrap() + 1, omega as usize);
            assert!(
                g.edges()
                    .iter()
                    .all(|e| e.from == e.to || r.colors[e.from] != r.colors[e.to])
            );
            for &a in &r.maximum_clique {
                for &b in &r.maximum_clique {
                    assert!(a == b || adj[a] >> b & 1 != 0);
                }
            }
        }
    }
}
