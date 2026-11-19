use std::fmt::Write as _;
use std::io::Write as _;
use std::process::{Command, Stdio};
use std::time::Duration;

use graph_library::{Graph, GraphKind};

pub struct Reference {
    pub distances: Option<Vec<Option<i128>>>,
    pub elapsed: Duration,
}

pub fn networkx(cases: &[(Graph<i64>, usize)], repeats: usize) -> (String, Vec<[Reference; 3]>) {
    let mut input = String::new();
    writeln!(input, "{}", cases.len()).unwrap();
    for (graph, source) in cases {
        let directed = u8::from(graph.kind() == GraphKind::Directed);
        writeln!(
            input,
            "{} {} {source} {directed}",
            graph.node_count(),
            graph.edge_count()
        )
        .unwrap();
        for edge in graph.edges() {
            writeln!(input, "{} {} {}", edge.from, edge.to, edge.weight).unwrap();
        }
    }
    let python = std::env::var_os("NETWORKX_PYTHON").unwrap_or_else(|| "python3".into());
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/support/networkx_shortest.py"
    );
    let mut process = Command::new(python)
        .arg(path)
        .arg(repeats.to_string())
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
    let results = cases
        .iter()
        .map(|(graph, _)| {
            let n = graph.node_count();
            let mut mode = 0;
            std::array::from_fn(|_| {
                let mut fields = lines.next().unwrap().split_whitespace();
                let elapsed = Duration::from_nanos(fields.next().unwrap().parse().unwrap());
                let count = if mode == 0 { n } else { n * n };
                mode += 1;
                let first = fields.next().unwrap();
                let distances = (first != "cycle").then(|| {
                    let distances: Vec<_> = std::iter::once(first)
                        .chain(fields)
                        .map(|value| (value != "x").then(|| value.parse::<i128>().unwrap()))
                        .collect();
                    assert_eq!(distances.len(), count);
                    distances
                });
                Reference { distances, elapsed }
            })
        })
        .collect();
    assert!(lines.next().is_none());
    (version, results)
}
