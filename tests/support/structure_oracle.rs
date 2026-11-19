use std::fmt::Write as _;
use std::io::Write as _;
use std::process::{Command, Stdio};

use graph_library::Graph;

pub fn graph<W>(buffer: &mut String, graph: &Graph<W>, weight: impl Fn(&W) -> f64) {
    writeln!(
        buffer,
        "{} {} {}",
        graph.node_count(),
        graph.edge_count(),
        u8::from(graph.is_directed())
    )
    .unwrap();
    for edge in graph.edges() {
        writeln!(
            buffer,
            "{} {} {:.17e}",
            edge.from,
            edge.to,
            weight(&edge.weight)
        )
        .unwrap();
    }
}

pub fn invoke(mode: &str, input: &str) -> (String, Vec<Vec<String>>) {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/support/networkx_structure.py"
    );
    let executable = std::env::var_os("NETWORKX_PYTHON").unwrap_or_else(|| "python3".into());
    let mut process = Command::new(executable)
        .arg(path)
        .arg(mode)
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
    assert!(output.status.success(), "NetworkX oracle failed");
    let text = String::from_utf8(output.stdout).unwrap();
    let mut lines = text.lines();
    let version = lines.next().unwrap().to_owned();
    (
        version,
        lines
            .map(|line| line.split_whitespace().map(str::to_owned).collect())
            .collect(),
    )
}
