use std::fmt::Write as _;
use std::io::Write as _;
use std::process::{Command, Stdio};
use std::time::Duration;

use graph_library::Graph;

pub struct Reference {
    pub size: usize,
    pub weight: i128,
    pub elapsed: Duration,
}

pub fn networkx(graphs: &[Graph<i64>], repetitions: usize) -> (String, Vec<[Reference; 3]>) {
    let mut input = String::new();
    writeln!(input, "{}", graphs.len()).unwrap();
    for graph in graphs {
        writeln!(input, "{} {}", graph.node_count(), graph.edge_count()).unwrap();
        for edge in graph.edges() {
            writeln!(input, "{} {} {}", edge.from, edge.to, edge.weight).unwrap();
        }
    }
    let python = std::env::var_os("NETWORKX_PYTHON").unwrap_or_else(|| "python3".into());
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/support/networkx_matching.py"
    );
    let mut process = Command::new(python)
        .arg(path)
        .arg(repetitions.to_string())
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("start Python with NetworkX installed");
    process
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    let output = process.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "NetworkX oracle exited with {}",
        output.status
    );
    let output = String::from_utf8(output.stdout).unwrap();
    let mut lines = output.lines();
    let version = lines.next().unwrap().to_owned();
    let data = (0..graphs.len())
        .map(|_| {
            std::array::from_fn(|_| {
                let mut row = lines.next().unwrap().split_whitespace();
                Reference {
                    size: row.next().unwrap().parse().unwrap(),
                    weight: row.next().unwrap().parse().unwrap(),
                    elapsed: Duration::from_nanos(row.next().unwrap().parse().unwrap()),
                }
            })
        })
        .collect();
    assert!(lines.next().is_none());
    (version, data)
}
