//! Compare preprocessing and a fixed mix of positive and negative queries.
use graph_library::connectivity::{ChainReachabilityIndex, ReachabilityIndex, st_connected};
use graph_library::{Graph, GraphKind};
use std::{hint::black_box, time::Instant};

fn measure(f: impl FnOnce() -> Vec<bool>) -> (Vec<bool>, f64) {
    let start = Instant::now();
    let results = black_box(f());
    (results, start.elapsed().as_secs_f64() * 1000.)
}
fn main() {
    let n = 100_000;
    let mut state = 42_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state as usize % n
    };
    let queries: Vec<_> = (0..2000).map(|_| (next(), next())).collect();
    println!(
        "family,vertices,arcs,queries,positive,bfs_ms,interval_build_ms,interval_query_ms,chain_build_ms,chain_query_ms,chains"
    );
    for (family, reverse) in [("lanes", false), ("lanes with SCCs", true)] {
        let mut b = Graph::builder(n, GraphKind::Directed);
        for u in 0..n - 4 {
            b.add_edge(u, u + 4, ()).unwrap();
            if reverse && (u / 4) % 16 != 15 {
                b.add_edge(u + 4, u, ()).unwrap();
            }
            if u % 113 == 0 && u + 5 < n {
                b.add_edge(u, u + 5, ()).unwrap();
            }
        }
        let g = b.build();
        let (expected, bfs) = measure(|| {
            queries
                .iter()
                .map(|&(u, v)| st_connected(&g, u, v).unwrap())
                .collect()
        });
        let start = Instant::now();
        let intervals = ReachabilityIndex::new(&g, 3, 17);
        let interval_build = start.elapsed().as_secs_f64() * 1000.;
        let (answers, interval_query) = measure(|| {
            queries
                .iter()
                .map(|&(u, v)| intervals.reachable(u, v).unwrap())
                .collect()
        });
        assert_eq!(expected, answers);
        let start = Instant::now();
        let chains = ChainReachabilityIndex::with_budget(&g, 128 * 1024 * 1024).unwrap();
        let chain_build = start.elapsed().as_secs_f64() * 1000.;
        let (answers, chain_query) = measure(|| {
            queries
                .iter()
                .map(|&(u, v)| chains.reachable(u, v).unwrap())
                .collect()
        });
        assert_eq!(expected, answers);
        let positive = expected.iter().filter(|&&b| b).count();
        println!(
            "{family},{n},{},{},{positive},{bfs:.3},{interval_build:.3},{interval_query:.3},{chain_build:.3},{chain_query:.3},{}",
            g.edge_count(),
            queries.len(),
            chains.chain_count()
        );
    }
}
