#[path = "../tests/support/structure_oracle.rs"]
mod oracle;

use std::fmt::Write as _;
use std::hint::black_box;
use std::time::{Duration, Instant};

use graph_library::functions::induced_subgraph;
use graph_library::generators::{
    complete_bipartite_graph, complete_graph, gnm_random_graph, grid_2d_graph, shrikhande_graph,
    tutte_graph,
};
use graph_library::isomorphism::vf2pp_subgraph_isomorphism;
use graph_library::planarity::check_planarity_with_witness;
use graph_library::{Graph, GraphBuilder, GraphKind};

const SAMPLES: usize = 7;
const UNDIRECTED: GraphKind = GraphKind::Undirected;

struct Measurement {
    algorithm: &'static str,
    family: &'static str,
    vertices: usize,
    edges: usize,
    pattern_vertices: usize,
    pattern_edges: usize,
    elapsed: Duration,
    outcome: bool,
    witness_edges: usize,
}

fn timing<T>(mut operation: impl FnMut() -> T, check: impl Fn(&T)) -> Duration {
    check(&operation());
    let mut samples = [Duration::ZERO; SAMPLES];
    for sample in &mut samples {
        let start = Instant::now();
        let result = black_box(operation());
        *sample = start.elapsed();
        check(&result);
    }
    samples.sort_unstable();
    samples[SAMPLES / 2]
}

fn sparse_witness_host(copies: usize) -> Graph<()> {
    let planar = tutte_graph();
    let mut builder = GraphBuilder::new(planar.node_count() * copies + 6, UNDIRECTED);
    let n = planar.node_count();
    for layer in 0..copies {
        for edge in planar.edges() {
            builder
                .add_edge(n * layer + edge.from, n * layer + edge.to, ())
                .unwrap();
        }
        builder.add_edge(layer * n, (layer + 1) * n, ()).unwrap();
    }
    for a in 0..3 {
        for b in 3..6 {
            builder
                .add_edge(n * copies + a, n * copies + b, ())
                .unwrap();
        }
    }
    builder.build()
}

fn diamond_subdivision(cells: usize) -> Graph<()> {
    let base = complete_bipartite_graph(3, 3).unwrap();
    let mut builder = GraphBuilder::new(6 + base.edge_count() * (3 * cells - 1), UNDIRECTED);
    let mut next = 6;
    for edge in base.edges() {
        let mut u = edge.from;
        for i in 0..cells {
            let (a, b) = (next, next + 1);
            next += 2;
            let v = if i + 1 == cells {
                edge.to
            } else {
                next += 1;
                next - 1
            };
            for (x, y) in [(u, a), (u, b), (a, b), (a, v), (b, v)] {
                builder.add_edge(x, y, ()).unwrap();
            }
            u = v;
        }
    }
    builder.build()
}

fn sparse_induced_host(side: usize) -> (Graph<()>, Graph<()>) {
    let grid = grid_2d_graph(side, side, [false; 2]).unwrap();
    let mut motif = GraphBuilder::new(8, UNDIRECTED);
    for u in 0..4 {
        for v in u + 1..4 {
            motif.add_edge(u, v, ()).unwrap();
        }
    }
    for u in 3..7 {
        motif.add_edge(u, u + 1, ()).unwrap();
    }
    let motif = motif.build();
    let n = grid.node_count();
    let mut host = GraphBuilder::new(n + motif.node_count(), UNDIRECTED);
    for edge in grid.edges() {
        host.add_edge(edge.from, edge.to, ()).unwrap();
    }
    for edge in motif.edges() {
        host.add_edge(n + edge.from, n + edge.to, ()).unwrap();
    }
    host.add_edge(n - 1, n, ()).unwrap();
    (host.build(), motif)
}

fn dense_induced_host(n: usize, k: usize) -> (Graph<()>, Graph<()>) {
    let host = gnm_random_graph(n, n * (n - 1) / 3, UNDIRECTED, 5061).unwrap();
    let selected = (0..k).map(|i| i * (n / k));
    let pattern = induced_subgraph(&host, selected)
        .unwrap()
        .materialize()
        .graph;
    (host, pattern)
}

fn rook() -> Graph<()> {
    let mut builder = GraphBuilder::new(16, UNDIRECTED);
    for u in 0..16 {
        for v in u + 1..16 {
            if u / 4 == v / 4 || u % 4 == v % 4 {
                builder.add_edge(u, v, ()).unwrap();
            }
        }
    }
    builder.build()
}

fn with_clique(base: &Graph<()>, layers: usize, isolates: usize) -> Graph<()> {
    let n = base.node_count();
    let mut builder = GraphBuilder::new(n * layers + isolates, UNDIRECTED);
    for i in 0..layers {
        for edge in base.edges() {
            builder
                .add_edge(n * i + edge.from, n * i + edge.to, ())
                .unwrap();
        }
        for j in 0..i {
            for u in 0..n {
                builder.add_edge(n * i + u, n * j + u, ()).unwrap();
            }
        }
    }
    builder.build()
}

fn main() {
    let compare = std::env::args().any(|arg| arg == "--networkx");
    let mut input = String::new();
    let mut samples = Vec::new();
    for (family, graph) in [
        ("sparse connected blocks", sparse_witness_host(4)),
        ("sparse connected blocks", sparse_witness_host(8)),
        ("dense complete", complete_graph(32, UNDIRECTED).unwrap()),
        ("dense complete", complete_graph(64, UNDIRECTED).unwrap()),
        ("worst-case diamond subdivision", diamond_subdivision(4)),
        ("worst-case diamond subdivision", diamond_subdivision(8)),
    ] {
        let witness = check_planarity_with_witness(&graph).unwrap_err();
        let elapsed = timing(
            || check_planarity_with_witness(black_box(&graph)).unwrap_err(),
            |result| assert_eq!(result, &witness),
        );
        if compare {
            input.push_str("C\n");
            oracle::graph(&mut input, &graph, |_| 1.0);
            writeln!(input, "{}", witness.edge_ids().count()).unwrap();
            for id in witness.edge_ids() {
                let edge = graph.edges()[id];
                writeln!(input, "{} {}", edge.from, edge.to).unwrap();
            }
        }
        samples.push(Measurement {
            algorithm: "kuratowski",
            family,
            vertices: graph.node_count(),
            edges: graph.edge_count(),
            pattern_vertices: 0,
            pattern_edges: 0,
            elapsed,
            outcome: false,
            witness_edges: witness.edge_ids().count(),
        });
    }

    let small = sparse_induced_host(24);
    let large = sparse_induced_host(40);
    let dense_small = dense_induced_host(96, 7);
    let dense_large = dense_induced_host(192, 8);
    for (family, (host, pattern), expected) in [
        ("sparse planted motif", small, true),
        ("sparse planted motif", large, true),
        ("dense planted pattern", dense_small, true),
        ("dense planted pattern", dense_large, true),
        (
            "adversarial strongly regular",
            (with_clique(&rook(), 1, 12), shrikhande_graph()),
            false,
        ),
        (
            "adversarial strongly regular product",
            (
                with_clique(&rook(), 3, 16),
                with_clique(&shrikhande_graph(), 3, 0),
            ),
            false,
        ),
    ] {
        let elapsed = timing(
            || vf2pp_subgraph_isomorphism(black_box(&host), black_box(&pattern)),
            |result| assert_eq!(result.is_some(), expected),
        );
        if compare {
            input.push_str("S\n");
            oracle::graph(&mut input, &host, |_| 1.0);
            oracle::graph(&mut input, &pattern, |_| 1.0);
            input.push_str("0\n");
        }
        samples.push(Measurement {
            algorithm: "induced vf2pp",
            family,
            vertices: host.node_count(),
            edges: host.edge_count(),
            pattern_vertices: pattern.node_count(),
            pattern_edges: pattern.edge_count(),
            elapsed,
            outcome: expected,
            witness_edges: 0,
        });
    }
    let oracle = compare.then(|| {
        let (version, rows) = oracle::invoke("3", &format!("{}\n{input}", samples.len()));
        eprintln!("NetworkX {version}, three samples plus a warm-up");
        assert_eq!(rows.len(), samples.len());
        rows
    });
    println!(
        "algorithm,family,vertices,edges,pattern_vertices,pattern_edges,rust_ms,networkx_ms,outcome,rust_witness_edges,networkx_witness_edges"
    );
    for (index, sample) in samples.iter().enumerate() {
        let nx = oracle.as_ref().map(|rows| &rows[index]);
        let baseline = nx.map_or(String::new(), |row| {
            assert_eq!(sample.outcome, row[1] == "1");
            format!("{:.3}", row[0].parse::<f64>().unwrap() / 1e6)
        });
        let nx_witness = nx.and_then(|row| row.get(2)).map_or("", String::as_str);
        println!(
            "{},{},{},{},{},{},{:.3},{},{},{},{}",
            sample.algorithm,
            sample.family,
            sample.vertices,
            sample.edges,
            sample.pattern_vertices,
            sample.pattern_edges,
            sample.elapsed.as_secs_f64() * 1e3,
            baseline,
            sample.outcome,
            sample.witness_edges,
            nx_witness
        );
    }
}
