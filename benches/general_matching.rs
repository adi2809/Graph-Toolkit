use std::collections::HashSet;
use std::hint::black_box;
use std::time::{Duration, Instant};

use graph_library::matching::{
    GeneralMatching, maximum_general_matching, maximum_general_matching_gabow,
};
use graph_library::{Graph, GraphBuilder, GraphError, GraphKind};

const SAMPLES: usize = 7;

struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn index(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn sparse(n: usize, seed: u64) -> Graph<()> {
    let mut random = Random(seed);
    let m = n * 3;
    let mut used = HashSet::with_capacity(m);
    let mut graph = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, m);
    while used.len() < m {
        let a = random.index(n);
        let b = random.index(n);
        if a != b && used.insert((a.min(b), a.max(b))) {
            graph.add_edge(a, b, ()).unwrap();
        }
    }
    graph.build()
}

fn dense(n: usize, seed: u64) -> Graph<()> {
    let mut random = Random(seed);
    let mut graph = Graph::builder(n, GraphKind::Undirected);
    for a in 0..n {
        for b in a + 1..n {
            if random.index(3) == 0 {
                graph.add_edge(a, b, ()).unwrap();
            }
        }
    }
    graph.build()
}

fn adversarial(scale: usize) -> Graph<()> {
    // The short-and-long chain family from Mehlhorn and Nobahari (2026).
    // Give internal vertices lower IDs so both greedy initializers leave
    // exactly two free vertices per chain.
    let clique = 2 * scale.isqrt();
    let sizes: Vec<_> = (0..scale / 8)
        .map(|_| 8)
        .chain((5..scale.isqrt()).map(|i| 2 * i))
        .collect();
    let n = clique + sizes.iter().sum::<usize>();
    let mut graph = Graph::builder(n, GraphKind::Undirected);
    for a in 0..clique {
        for b in a + 1..clique {
            graph.add_edge(a, b, ()).unwrap();
        }
    }
    let mut first = clique;
    for size in sizes {
        let last = first + size - 1;
        let left = last - 1;
        for a in first..left - 1 {
            graph.add_edge(a, a + 1, ()).unwrap();
        }
        graph.add_edge(left, first, ()).unwrap();
        graph.add_edge(left - 1, last, ()).unwrap();
        graph.add_edge(0, left, ()).unwrap();
        for a in first..left {
            graph.add_edge(0, a, ()).unwrap();
        }
        first += size;
    }
    graph.build()
}

type Matcher = fn(&Graph<()>) -> Result<GeneralMatching, GraphError>;

fn sample(matcher: Matcher, graph: &Graph<()>, expected: Option<usize>) -> (usize, Duration) {
    let start = Instant::now();
    let matching = black_box(matcher(black_box(graph)).unwrap());
    let elapsed = start.elapsed();

    let mut seen = vec![false; graph.node_count()];
    let mut count = 0;
    for id in matching.edges() {
        let edge = graph.edge(id).unwrap();
        assert!(edge.from != edge.to && !seen[edge.from] && !seen[edge.to]);
        seen[edge.from] = true;
        seen[edge.to] = true;
        let forward = matching.mate(edge.from).unwrap();
        let backward = matching.mate(edge.to).unwrap();
        assert_eq!((forward.node, forward.edge), (edge.to, id));
        assert_eq!((backward.node, backward.edge), (edge.from, id));
        count += 1;
    }
    assert_eq!(count, matching.len());
    if let Some(expected) = expected {
        assert_eq!(count, expected);
    }
    (count, elapsed)
}

fn bench(name: &str, graph: Graph<()>, known_size: Option<usize>) {
    let matchers = [
        maximum_general_matching as Matcher,
        maximum_general_matching_gabow as Matcher,
    ];
    let expected = sample(matchers[0], &graph, known_size).0;
    sample(matchers[1], &graph, Some(expected));

    let mut times = [[Duration::ZERO; 2]; SAMPLES];
    for (sample_id, row) in times.iter_mut().enumerate() {
        for offset in 0..2 {
            // Alternate who runs first to lessen cache and ordering effects.
            let algorithm = (sample_id + offset) % 2;
            row[algorithm] = sample(matchers[algorithm], &graph, Some(expected)).1;
        }
    }
    let [edmonds, gabow] = [0, 1].map(|algorithm| {
        let mut data = times.map(|row| row[algorithm]);
        data.sort_unstable();
        data[SAMPLES / 2].as_secs_f64() * 1e3
    });
    println!(
        "{name},{},{},{expected},{edmonds:.3},{gabow:.3},{:.2}",
        graph.node_count(),
        graph.edge_count(),
        edmonds / gabow
    );
}

fn main() {
    println!("family,vertices,edges,matched_edges,edmonds_median_ms,gabow_median_ms,speedup");
    for n in [6_000, 24_000, 60_000] {
        bench("sparse", sparse(n, 746_219), None);
    }
    for n in [400, 800, 1_600] {
        bench("dense", dense(n, 4_752_091), Some(n / 2));
    }
    for scale in [1_200, 3_600, 9_600, 20_000] {
        let graph = adversarial(scale);
        let known_size = graph.node_count() / 2;
        bench("short-and-long chains", graph, Some(known_size));
    }
}
