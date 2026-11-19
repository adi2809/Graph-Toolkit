#[path = "../tests/support/shortest_oracle.rs"]
mod shortest_oracle;

use std::collections::HashSet;
use std::hint::black_box;
use std::time::{Duration, Instant};

use graph_library::shortest_path::{
    AllPairsShortestPaths, SignedShortestPaths, bellman_ford, duong_shortest_paths, floyd_warshall,
    johnson,
};
use graph_library::{Graph, GraphBuilder, GraphKind};

const SAMPLES: usize = 7;

struct Random(u64);

impl Random {
    fn below(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 as usize % n
    }

    fn potentials(&mut self, n: usize) -> Vec<i64> {
        (0..n).map(|_| self.below(3_000) as i64 - 1_500).collect()
    }
}

fn sparse(n: usize) -> Graph<i64> {
    let mut rng = Random(6_532_170);
    let potential = rng.potentials(n);
    let m = 4 * n;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Directed, m);
    let mut seen = HashSet::with_capacity(m);
    for a in 0..n {
        let b = (a + 1) % n;
        seen.insert((a, b));
        let cost = rng.below(41) as i64 + potential[b] - potential[a];
        builder.add_edge(a, b, cost).unwrap();
    }
    while seen.len() < m {
        let a = rng.below(n);
        let b = rng.below(n);
        if a != b && seen.insert((a, b)) {
            let cost = rng.below(41) as i64 + potential[b] - potential[a];
            builder.add_edge(a, b, cost).unwrap();
        }
    }
    builder.build()
}

fn dense(n: usize) -> Graph<i64> {
    let mut rng = Random(7_482_307);
    let potential = rng.potentials(n);
    let mut builder = Graph::builder(n, GraphKind::Directed);
    for a in 0..n {
        for b in 0..n {
            if a != b && rng.below(3) != 0 {
                let cost = rng.below(41) as i64 + potential[b] - potential[a];
                builder.add_edge(a, b, cost).unwrap();
            }
        }
    }
    builder.build()
}

fn worst_case(n: usize) -> Graph<i64> {
    let hub = n / 2;
    let mut builder = Graph::builder(n, GraphKind::Directed);
    // Every step down the negative chain sends a better label into the
    // zero-cost clique. Its dense adjacency is scanned in Θ(V) rounds.
    for a in 0..hub {
        if a + 1 < hub {
            builder.add_edge(a, a + 1, -1).unwrap();
        }
        builder.add_edge(a, hub, 0).unwrap();
    }
    builder.add_edge(hub, 0, hub as i64).unwrap();
    for a in hub..n {
        for b in hub..n {
            if a != b {
                builder.add_edge(a, b, 0).unwrap();
            }
        }
    }
    builder.build()
}

enum ResultPaths {
    Source(SignedShortestPaths),
    Matrix(AllPairsShortestPaths),
}

impl ResultPaths {
    fn values(&self) -> Vec<Option<i128>> {
        match self {
            Self::Source(paths) => paths.distances().collect(),
            Self::Matrix(paths) => (0..paths.node_count())
                .flat_map(|source| paths.distances_from(source).unwrap())
                .collect(),
        }
    }
}

fn solve(graph: &Graph<i64>, source: usize, algorithm: usize) -> ResultPaths {
    match algorithm {
        0 => ResultPaths::Source(bellman_ford(graph, source).unwrap()),
        1 => ResultPaths::Matrix(johnson(graph).unwrap()),
        _ => ResultPaths::Matrix(floyd_warshall(graph).unwrap()),
    }
}

fn compare_duong() {
    type Generator = fn(usize) -> Graph<i64>;
    type Solver = fn(
        &Graph<i64>,
        usize,
    )
        -> Result<SignedShortestPaths, graph_library::shortest_path::SignedPathError>;
    let solvers: [Solver; 2] = [bellman_ford, duong_shortest_paths];
    println!("family,vertices,edges,bellman_ford_ms,duong_ms");
    let instances: [(&str, [usize; 3], Generator); 3] = [
        ("sparse", [1024, 8192, 32768], sparse),
        ("dense", [128, 256, 512], dense),
        ("worst-case chain and clique", [256, 512, 1024], worst_case),
    ];
    for (family, sizes, generate) in instances {
        for n in sizes {
            let graph = generate(n);
            let expected = bellman_ford(&graph, 0)
                .unwrap()
                .distances()
                .collect::<Vec<_>>();
            let mut samples = [[Duration::ZERO; SAMPLES]; 2];
            for round in 0..=SAMPLES {
                for algorithm in [round % 2, (round + 1) % 2] {
                    let start = Instant::now();
                    let result = black_box(solvers[algorithm](black_box(&graph), 0).unwrap());
                    let elapsed = start.elapsed();
                    assert!(result.distances().eq(expected.iter().copied()));
                    if round > 0 {
                        samples[algorithm][round - 1] = elapsed;
                    }
                }
            }
            let [bellman, duong] = samples.map(|mut values| {
                values.sort_unstable();
                values[SAMPLES / 2].as_secs_f64() * 1000.0
            });
            println!(
                "{family},{n},{},{bellman:.3},{duong:.3}",
                graph.edge_count()
            );
        }
    }
}

fn main() {
    if std::env::args().any(|argument| argument == "--duong") {
        compare_duong();
        return;
    }

    let compare = std::env::args().any(|argument| argument == "--networkx");
    let names = ["Bellman-Ford", "Johnson", "Floyd-Warshall"];
    let mut families = Vec::new();
    let mut cases = Vec::new();
    for (family, sizes, generate) in [
        ("sparse", [96, 192, 384], sparse as fn(usize) -> Graph<i64>),
        ("dense", [64, 128, 224], dense),
        ("worst-case chain and clique", [64, 128, 192], worst_case),
    ] {
        for n in sizes {
            families.push(family);
            cases.push((generate(n), 0));
        }
    }

    println!("family,vertices,edges,algorithm,rust_median_ms,networkx_median_ms");
    let mut measurements = Vec::new();
    for (graph, source) in &cases {
        let expected: [Vec<Option<i128>>; 3] =
            std::array::from_fn(|mode| solve(graph, *source, mode).values());
        assert_eq!(expected[1], expected[2]);
        assert_eq!(expected[0], expected[1][..graph.node_count()]);
        let mut samples = [[Duration::ZERO; 3]; SAMPLES];
        for (round, row) in samples.iter_mut().enumerate() {
            for offset in 0..3 {
                let algorithm = (round + offset) % 3;
                let start = Instant::now();
                let result = black_box(solve(black_box(graph), *source, algorithm));
                row[algorithm] = start.elapsed();
                assert_eq!(result.values(), expected[algorithm]);
            }
        }
        let median = [0, 1, 2].map(|mode| {
            let mut durations = samples.map(|row| row[mode]);
            durations.sort_unstable();
            durations[SAMPLES / 2].as_secs_f64() * 1e3
        });
        measurements.push((expected, median));
    }

    let reference = compare.then(|| {
        let (version, results) = shortest_oracle::networkx(&cases, 3);
        eprintln!("NetworkX {version}, three samples per workload");
        results
    });
    for (case, ((graph, _), (expected, rust))) in cases.iter().zip(measurements).enumerate() {
        for (mode, name) in names.iter().enumerate() {
            let nx = reference.as_ref().map_or(String::new(), |all| {
                let row = &all[case][mode];
                assert_eq!(row.distances.as_ref(), Some(&expected[mode]));
                format!("{:.3}", row.elapsed.as_secs_f64() * 1e3)
            });
            println!(
                "{},{},{},{name},{:.3},{nx}",
                families[case],
                graph.node_count(),
                graph.edge_count(),
                rust[mode]
            );
        }
    }
}
