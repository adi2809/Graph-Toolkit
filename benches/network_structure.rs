#[path = "../tests/support/structure_oracle.rs"]
mod oracle;

use std::fmt::Write as _;
use std::hint::black_box;
use std::time::{Duration, Instant};

use graph_library::community::{LouvainOptions, louvain_communities_by};
use graph_library::generators::{
    complete_graph, gnm_random_graph, grid_2d_graph, ring_of_cliques, shrikhande_graph,
    stochastic_block_model,
};
use graph_library::isomorphism::vf2pp_isomorphism;
use graph_library::planarity::check_planarity;
use graph_library::{Graph, GraphBuilder, GraphKind};

const SAMPLES: usize = 7;
const UNDIRECTED: GraphKind = GraphKind::Undirected;

struct Measurement {
    algorithm: &'static str,
    family: &'static str,
    vertices: usize,
    edges: usize,
    rust_elapsed: Duration,
    result: Expected,
}

enum Expected {
    Boolean(bool),
    Modularity {
        score: f64,
        count: usize,
        levels: usize,
        sweeps: usize,
    },
}

fn time<T>(mut run: impl FnMut() -> T, validate: impl Fn(&T)) -> Duration {
    validate(&run());
    let mut samples = [Duration::ZERO; SAMPLES];
    for sample in &mut samples {
        let start = Instant::now();
        let actual = black_box(run());
        *sample = start.elapsed();
        validate(&actual);
    }
    samples.sort_unstable();
    samples[SAMPLES / 2]
}

fn stacked_triangulation(n: usize) -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(n, UNDIRECTED, 3 * n - 6);
    for u in 1..n {
        for v in u.saturating_sub(3)..u {
            builder.add_edge(u, v, ()).unwrap();
        }
    }
    builder.build()
}

fn relabel(graph: &Graph<()>) -> Graph<()> {
    let n = graph.node_count();
    let mut builder = GraphBuilder::with_edge_capacity(n, graph.kind(), graph.edge_count());
    let permute = |node| n - 1 - ((node + n / 3) % n);
    for edge in graph.edges().iter().rev() {
        builder
            .add_edge(permute(edge.from), permute(edge.to), ())
            .unwrap();
    }
    builder.build()
}

fn rook() -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(16, UNDIRECTED, 48);
    for u in 0..16 {
        for v in u + 1..16 {
            if u / 4 == v / 4 || u % 4 == v % 4 {
                builder.add_edge(u, v, ()).unwrap();
            }
        }
    }
    builder.build()
}

fn with_clique(graph: &Graph<()>, layers: usize) -> Graph<()> {
    let n = graph.node_count();
    let mut builder = GraphBuilder::new(n * layers, UNDIRECTED);
    for i in 0..layers {
        for edge in graph.edges() {
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

fn blocks(n: usize, block_size: usize, inside: f64, outside: f64, kind: GraphKind) -> Graph<()> {
    let k = n / block_size;
    let sizes = vec![block_size; k];
    let probabilities: Vec<_> = (0..k)
        .map(|i| {
            (0..k)
                .map(|j| if i == j { inside } else { outside })
                .collect()
        })
        .collect();
    stochastic_block_model(&sizes, &probabilities, kind, 884921).unwrap()
}

fn main() {
    let compare = std::env::args().any(|arg| arg == "--networkx");
    let mut oracle_input = String::new();
    let mut measurements = Vec::new();

    for (case, graph, expected) in [
        (
            "sparse grid",
            grid_2d_graph(64, 64, [false; 2]).unwrap(),
            true,
        ),
        (
            "sparse grid",
            grid_2d_graph(128, 128, [false; 2]).unwrap(),
            true,
        ),
        (
            "dense complete",
            complete_graph(128, UNDIRECTED).unwrap(),
            false,
        ),
        (
            "dense complete",
            complete_graph(512, UNDIRECTED).unwrap(),
            false,
        ),
        (
            "worst-case stacked triangulation",
            stacked_triangulation(4096),
            true,
        ),
        (
            "worst-case stacked triangulation",
            stacked_triangulation(16384),
            true,
        ),
    ] {
        let elapsed = time(
            || check_planarity(black_box(&graph)),
            |result| assert_eq!(result.is_some(), expected),
        );
        if compare {
            oracle_input.push_str("P\n");
            oracle::graph(&mut oracle_input, &graph, |_| 1.0);
        }
        measurements.push(Measurement {
            algorithm: "planarity",
            family: case,
            vertices: graph.node_count(),
            edges: graph.edge_count(),
            rust_elapsed: elapsed,
            result: Expected::Boolean(expected),
        });
    }

    let sparse = |n, kind| gnm_random_graph(n, 4 * n, kind, 304911).unwrap();
    let dense = |n, kind| gnm_random_graph(n, n * (n - 1) / 6, kind, 906511).unwrap();
    for (case, left, right, expected) in [
        ("sparse undirected", sparse(1024, UNDIRECTED), None, true),
        (
            "sparse directed",
            sparse(8192, GraphKind::Directed),
            None,
            true,
        ),
        ("dense undirected", dense(128, UNDIRECTED), None, true),
        (
            "dense directed",
            dense(512, GraphKind::Directed),
            None,
            true,
        ),
        (
            "adversarial strongly regular",
            rook(),
            Some(shrikhande_graph()),
            false,
        ),
        (
            "adversarial strongly regular product",
            with_clique(&rook(), 3),
            Some(with_clique(&shrikhande_graph(), 3)),
            false,
        ),
    ] {
        let right = right.unwrap_or_else(|| relabel(&left));
        let elapsed = time(
            || vf2pp_isomorphism(black_box(&left), black_box(&right)),
            |mapping| assert_eq!(mapping.is_some(), expected),
        );
        if compare {
            oracle_input.push_str("I\n");
            oracle::graph(&mut oracle_input, &left, |_| 1.0);
            oracle::graph(&mut oracle_input, &right, |_| 1.0);
            oracle_input.push_str("0\n");
        }
        measurements.push(Measurement {
            algorithm: "vf2pp",
            family: case,
            vertices: left.node_count(),
            edges: left.edge_count(),
            rust_elapsed: elapsed,
            result: Expected::Boolean(expected),
        });
    }

    let options = LouvainOptions {
        seed: 613,
        ..Default::default()
    };
    for (case, graph) in [
        ("sparse blocks", blocks(2048, 128, 0.04, 0.0005, UNDIRECTED)),
        (
            "sparse directed blocks",
            blocks(8192, 256, 0.05, 0.00025, GraphKind::Directed),
        ),
        ("dense blocks", blocks(256, 64, 0.8, 0.3, UNDIRECTED)),
        (
            "dense directed blocks",
            blocks(512, 128, 0.8, 0.3, GraphKind::Directed),
        ),
        (
            "adversarial ring of cliques",
            ring_of_cliques(192, 5).unwrap(),
        ),
        (
            "adversarial ring of cliques",
            ring_of_cliques(1024, 5).unwrap(),
        ),
    ] {
        let run = || louvain_communities_by(black_box(&graph), |_| 1.0, options).unwrap();
        let expected = run();
        let elapsed = time(run, |result| assert_eq!(result, &expected));
        if compare {
            oracle_input.push_str("L\n");
            oracle::graph(&mut oracle_input, &graph, |_| 1.0);
            writeln!(
                oracle_input,
                "{} {} {} 0",
                options.resolution, options.threshold, options.seed
            )
            .unwrap();
            for id in &expected.community_of {
                write!(oracle_input, "{id} ").unwrap();
            }
            oracle_input.push('\n');
        }
        measurements.push(Measurement {
            algorithm: "louvain",
            family: case,
            vertices: graph.node_count(),
            edges: graph.edge_count(),
            rust_elapsed: elapsed,
            result: Expected::Modularity {
                score: expected.modularity,
                count: expected.community_count(),
                levels: expected.level_modularities.len(),
                sweeps: expected.sweeps,
            },
        });
    }

    let reference = compare.then(|| {
        let request = format!("{}\n{oracle_input}", measurements.len());
        let (version, rows) = oracle::invoke("3", &request);
        eprintln!("NetworkX {version}, three samples per case");
        assert_eq!(measurements.len(), rows.len());
        rows
    });

    println!(
        "algorithm,family,vertices,edges,rust_median_ms,networkx_median_ms,result,rust_modularity,networkx_modularity,rust_communities,networkx_communities,rust_levels,rust_sweeps"
    );
    for (i, row) in measurements.into_iter().enumerate() {
        let nx = reference.as_ref().map(|rows| &rows[i]);
        let nx_millis = nx.map_or(String::new(), |row| {
            format!("{:.3}", row[0].parse::<f64>().unwrap() / 1e6)
        });
        let millis = row.rust_elapsed.as_secs_f64() * 1e3;
        let mut details = String::new();
        match row.result {
            Expected::Boolean(value) => {
                if let Some(nx) = nx {
                    assert_eq!(value, nx[1] == "1");
                }
                write!(details, "{value},,,,,,").unwrap();
            }
            Expected::Modularity {
                score,
                count,
                levels,
                sweeps,
            } => {
                let (mut nx_score, mut nx_count) = ("", "");
                if let Some(nx) = nx {
                    assert!((nx[3].parse::<f64>().unwrap() - score).abs() < 1e-10);
                    (nx_score, nx_count) = (&nx[1], &nx[2]);
                }
                write!(
                    details,
                    ",{score:.9},{nx_score},{count},{nx_count},{levels},{sweeps}"
                )
                .unwrap();
            }
        }
        println!(
            "{},{},{},{},{millis:.3},{nx_millis},{details}",
            row.algorithm, row.family, row.vertices, row.edges,
        );
    }
}
