#[path = "../tests/support/weighted_oracle.rs"]
mod weighted_oracle;

use std::collections::HashSet;
use std::hint::black_box;
use std::time::{Duration, Instant};

use graph_library::matching::{WeightedMatching, maximum_weight_matching, minimum_cost_matching};
use graph_library::{Graph, GraphBuilder, GraphKind};

const SAMPLES: usize = 7;

struct Random(u64);

impl Random {
    fn below(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % n as u64) as usize
    }

    fn weight(&mut self) -> i64 {
        self.below(2_001) as i64 - 1_000
    }
}

fn sparse(n: usize) -> Graph<i64> {
    let mut rng = Random(6_906_193);
    let m = n * 3;
    let mut seen = HashSet::with_capacity(m);
    let mut graph = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, m);
    for a in 0..n {
        let b = (a + 1) % n;
        seen.insert((a.min(b), a.max(b)));
        graph.add_edge(a, b, rng.weight()).unwrap();
    }
    while seen.len() < m {
        let a = rng.below(n);
        let b = rng.below(n);
        if a != b && seen.insert((a.min(b), a.max(b))) {
            graph.add_edge(a, b, rng.weight()).unwrap();
        }
    }
    graph.build()
}

fn dense(n: usize) -> Graph<i64> {
    let mut rng = Random(56_009_347);
    let mut graph = Graph::builder(n, GraphKind::Undirected);
    for a in 0..n {
        for b in a + 1..n {
            if rng.below(3) != 0 {
                graph.add_edge(a, b, rng.weight()).unwrap();
            }
        }
    }
    graph.build()
}

fn worst_case(n: usize) -> Graph<i64> {
    // A complete graph with disjoint signed ladders of distinct heavy edges.
    // Each profitable rung forces a separate augmentation. Background edges
    // stay slack during those stages, so the free-vertex scans are cubic.
    let mut graph = Graph::builder(n, GraphKind::Undirected);
    for a in 0..n {
        for b in a + 1..n {
            let weight = if a % 2 == 0 && b == a + 1 {
                let rung = 10_000 + a as i64;
                if a < n / 2 { rung } else { -rung }
            } else {
                1
            };
            graph.add_edge(a, b, weight).unwrap();
        }
    }
    graph.build()
}

const OBJECTIVES: [&str; 3] = ["maximum weight", "cardinality first", "minimum cost"];

fn solve(graph: &Graph<i64>, mode: usize) -> WeightedMatching {
    match mode {
        0 => maximum_weight_matching(graph, false),
        1 => maximum_weight_matching(graph, true),
        _ => minimum_cost_matching(graph),
    }
    .unwrap()
}

fn verify(graph: &Graph<i64>, result: &WeightedMatching) -> (usize, i128) {
    let mut used = vec![false; graph.node_count()];
    let mut count = 0;
    let mut weight = 0;
    for id in result.matching.edges() {
        let edge = &graph.edges()[id];
        assert!(edge.from != edge.to && !used[edge.from] && !used[edge.to]);
        used[edge.from] = true;
        used[edge.to] = true;
        count += 1;
        weight += i128::from(edge.weight);
    }
    assert_eq!(count, result.matching.len());
    assert_eq!(weight, result.total_weight);
    (count, weight)
}

fn main() {
    let run_networkx = std::env::args().any(|argument| argument == "--networkx");
    let mut labels = Vec::new();
    let mut graphs = Vec::new();
    for (label, sizes, generator) in [
        ("sparse", [128, 256, 512], sparse as fn(usize) -> Graph<i64>),
        ("dense", [64, 128, 256], dense),
        ("worst-case ladders", [64, 128, 256], worst_case),
    ] {
        for n in sizes {
            labels.push(label);
            graphs.push(generator(n));
        }
    }

    println!(
        "graph,vertices,edges,objective,matching_size,value,rust_median_ms,networkx_median_ms"
    );
    let mut timed = Vec::new();
    for graph in &graphs {
        let baseline: [(usize, i128); 3] =
            std::array::from_fn(|mode| verify(graph, &solve(graph, mode)));
        let mut elapsed = [[Duration::ZERO; 3]; SAMPLES];
        for (sample, row) in elapsed.iter_mut().enumerate() {
            for offset in 0..3 {
                let mode = (sample + offset) % 3;
                let start = Instant::now();
                let answer = black_box(solve(black_box(graph), mode));
                row[mode] = start.elapsed();
                assert_eq!(verify(graph, &answer), baseline[mode]);
            }
        }
        let median: [f64; 3] = std::array::from_fn(|mode| {
            let mut data = elapsed.map(|row| row[mode]);
            data.sort_unstable();
            data[SAMPLES / 2].as_secs_f64() * 1e3
        });
        timed.push((baseline, median));
    }
    let oracle = run_networkx.then(|| {
        let (version, results) = weighted_oracle::networkx(&graphs, 3);
        eprintln!("NetworkX {version}, three samples per workload");
        results
    });
    for (case, (graph, (baseline, median))) in graphs.iter().zip(timed).enumerate() {
        for (mode, &name) in OBJECTIVES.iter().enumerate() {
            let (size, value) = baseline[mode];
            let nx = oracle.as_ref().map_or(String::new(), |all| {
                let expected = &all[case][mode];
                assert_eq!(value, expected.weight);
                if mode != 0 {
                    assert_eq!(size, expected.size);
                }
                format!("{:.3}", expected.elapsed.as_secs_f64() * 1e3)
            });
            println!(
                "{},{},{},{name},{size},{value},{:.3},{nx}",
                labels[case],
                graph.node_count(),
                graph.edge_count(),
                median[mode]
            );
        }
    }
}
