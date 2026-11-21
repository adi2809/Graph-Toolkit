use graph_library::{
    Graph, GraphError, GraphKind,
    connectivity::{DisjointSet, low_stretch_spanning_forest},
};
#[test]
fn multiscale_contractions_return_original_edges_and_span_each_component() {
    let mut r = 331u64;
    let mut next = || {
        r ^= r << 13;
        r ^= r >> 7;
        r ^= r << 17;
        r
    };
    for _ in 0..100 {
        let n = 61;
        let mut g = Graph::builder(n, GraphKind::Undirected);
        let mut forest = DisjointSet::new(n);
        for _ in 0..100 {
            let u = next() as usize % n;
            let v = next() as usize % n;
            g.add_edge(u, v, 1u64 << (next() % 64)).unwrap();
        }
        let g = g.build();
        let ids = low_stretch_spanning_forest(&g).unwrap();
        assert_eq!(ids, low_stretch_spanning_forest(&g).unwrap());
        for id in ids {
            let e = &g.edges()[id];
            assert_ne!(forest.find(e.from), forest.find(e.to));
            forest.union(e.from, e.to);
        }
        for e in g.edges() {
            assert_eq!(forest.find(e.from), forest.find(e.to));
        }
    }
    let mut b = Graph::builder(3, GraphKind::Undirected);
    b.add_edge(0, 2, u64::MAX).unwrap();
    b.add_edge(0, 1, 1).unwrap();
    b.add_edge(1, 2, 1).unwrap();
    assert_eq!(low_stretch_spanning_forest(&b.build()).unwrap(), [1, 2]);
    let mut b = Graph::builder(1, GraphKind::Undirected);
    b.add_edge(0, 0, 0).unwrap();
    assert_eq!(
        low_stretch_spanning_forest(&b.build()),
        Err(GraphError::InvalidModelParameters)
    );
}
