#[path = "support/shortest_oracle.rs"]
mod shortest_oracle;

use graph_library::shortest_path::{
    AllPairsShortestPaths, NegativeCycle, SignedPath, SignedPathError, bellman_ford,
    dag_longest_path, dag_shortest_paths, duong_shortest_paths, find_negative_cycle,
    floyd_warshall, johnson,
};
use graph_library::{Graph, GraphError, GraphKind};

fn graph(
    n: usize,
    kind: GraphKind,
    edges: impl IntoIterator<Item = (usize, usize, i64)>,
) -> Graph<i64> {
    let mut builder = Graph::builder(n, kind);
    for (a, b, cost) in edges {
        builder.add_edge(a, b, cost).unwrap();
    }
    builder.build()
}

struct Samples(u64);

impl Samples {
    fn index(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 as usize % n
    }
}

fn allowed_cases() -> Vec<(Graph<i64>, usize)> {
    let mut rng = Samples(675_903);
    (0..175)
        .map(|case| {
            let n = if case % 28 == 0 {
                33 + rng.index(41)
            } else {
                1 + rng.index(15)
            };
            let kind = if case % 4 == 0 {
                GraphKind::Undirected
            } else {
                GraphKind::Directed
            };
            let potential: Vec<_> = (0..n).map(|_| rng.index(100) as i64 - 50).collect();
            let source = rng.index(n);
            let count = rng.index(2 * n * n);
            let edges = (0..count).map(|_| {
                let a = rng.index(n);
                let b = rng.index(n);
                let base = rng.index(11) as i64;
                let cost = if kind == GraphKind::Directed {
                    base + potential[b] - potential[a]
                } else {
                    base
                };
                (a, b, cost)
            });
            (graph(n, kind, edges), source)
        })
        .collect()
}

fn cycle_cases() -> Vec<(Graph<i64>, usize)> {
    let directed = GraphKind::Directed;
    vec![
        (graph(1, directed, [(0, 0, i64::MIN)]), 0),
        (
            graph(
                7,
                directed,
                [(0, 1, 5), (2, 3, -7), (3, 4, 1), (4, 2, 3), (3, 2, 123)],
            ),
            0,
        ),
        (
            graph(
                7,
                directed,
                [(0, 1, 5), (1, 2, 2), (2, 3, -7), (3, 4, 1), (4, 2, 3)],
            ),
            0,
        ),
        (
            graph(4, GraphKind::Undirected, [(0, 1, 4), (2, 3, i64::MIN)]),
            2,
        ),
        (
            graph(
                3,
                directed,
                [(0, 1, i64::MAX), (1, 2, i64::MIN), (2, 0, -1)],
            ),
            0,
        ),
        (
            graph(
                70,
                directed,
                (0..69).map(|a| (a, a + 1, -1)).chain([(69, 0, 68)]),
            ),
            0,
        ),
    ]
}

fn check_walk(graph: &Graph<i64>, nodes: &[usize], edges: &[usize], cost: i128) {
    assert_eq!(nodes.len(), edges.len() + 1);
    let mut sum = 0;
    for (pair, &id) in nodes.windows(2).zip(edges) {
        let edge = graph.edge(id).unwrap();
        assert!(
            (edge.from, edge.to) == (pair[0], pair[1])
                || (graph.kind() == GraphKind::Undirected
                    && (edge.to, edge.from) == (pair[0], pair[1]))
        );
        sum += i128::from(edge.weight);
    }
    assert_eq!(cost, sum);
}

fn check_path(graph: &Graph<i64>, path: &SignedPath, from: usize, to: usize, cost: i128) {
    assert_eq!(path.nodes.first(), Some(&from));
    assert_eq!(path.nodes.last(), Some(&to));
    assert_eq!(path.cost, cost);
    let mut used = vec![false; graph.node_count()];
    for &node in &path.nodes {
        assert!(!std::mem::replace(&mut used[node], true));
    }
    check_walk(graph, &path.nodes, &path.edges, cost);
}

fn check_cycle(graph: &Graph<i64>, cycle: &NegativeCycle, source: Option<usize>) {
    assert!(!cycle.edges.is_empty());
    assert_eq!(cycle.nodes.first(), cycle.nodes.last());
    assert!(cycle.cost < 0);
    check_walk(graph, &cycle.nodes, &cycle.edges, cycle.cost);
    if let Some(source) = source {
        let reached = graph_library::traversal::bfs(graph, source).unwrap();
        assert!(reached.contains(cycle.nodes.first().unwrap()));
    }
}

fn extract_cycle(graph: &Graph<i64>, error: SignedPathError, source: Option<usize>) {
    let SignedPathError::NegativeCycle(cycle) = error else {
        panic!("unexpected error {error}")
    };
    check_cycle(graph, &cycle, source);
}

fn dp_distances(graph: &Graph<i64>, from: usize) -> Vec<Option<i128>> {
    let n = graph.node_count();
    let mut current = vec![None; n];
    current[from] = Some(0_i128);
    for _ in 1..n {
        let mut next = current.clone();
        for (node, &value) in current.iter().enumerate() {
            if let Some(cost) = value {
                for arc in graph.outgoing(node).unwrap() {
                    let candidate = cost + i128::from(graph.edges()[arc.edge].weight);
                    if next[arc.node].is_none_or(|old| candidate < old) {
                        next[arc.node] = Some(candidate);
                    }
                }
            }
        }
        if next == current {
            break;
        }
        current = next;
    }
    current
}

fn flatten(paths: &AllPairsShortestPaths) -> Vec<Option<i128>> {
    (0..paths.node_count())
        .flat_map(|from| paths.distances_from(from).unwrap())
        .collect()
}

#[test]
fn signed_algorithms_match_edge_budget_dp_and_return_valid_paths() {
    for (graph, source) in allowed_cases() {
        assert!(find_negative_cycle(&graph).is_none());
        let johnson = johnson(&graph).unwrap();
        let floyd = floyd_warshall(&graph).unwrap();
        assert_eq!(flatten(&johnson), flatten(&floyd));
        let duong = duong_shortest_paths(&graph, source).unwrap();
        let single = bellman_ford(&graph, source).unwrap();
        assert_eq!(
            duong.distances().collect::<Vec<_>>(),
            single.distances().collect::<Vec<_>>()
        );
        for target in 0..graph.node_count() {
            if let Some(path) = duong.path_to(target) {
                check_path(
                    &graph,
                    &path,
                    source,
                    target,
                    single.distance(target).unwrap(),
                );
            }
        }
        for from in 0..graph.node_count() {
            for (to, expected) in dp_distances(&graph, from).into_iter().enumerate() {
                assert_eq!(johnson.distance(from, to), expected);
                if from == source {
                    assert_eq!(single.distance(to), expected);
                }
                if let Some(cost) = expected {
                    check_path(&graph, &johnson.path(from, to).unwrap(), from, to, cost);
                    check_path(&graph, &floyd.path(from, to).unwrap(), from, to, cost);
                    if from == source {
                        check_path(&graph, &single.path_to(to).unwrap(), from, to, cost);
                    }
                }
            }
        }
    }
}

#[test]
fn negative_cycles_are_global_for_apsp_and_reachable_for_bellman_ford() {
    for (case, (graph, source)) in cycle_cases().into_iter().enumerate() {
        if case == 1 {
            assert_eq!(bellman_ford(&graph, source).unwrap().distance(1), Some(5));
            assert_eq!(
                duong_shortest_paths(&graph, source).unwrap().distance(1),
                Some(5)
            );
        } else {
            extract_cycle(
                &graph,
                duong_shortest_paths(&graph, source).unwrap_err(),
                Some(source),
            );
            extract_cycle(
                &graph,
                bellman_ford(&graph, source).unwrap_err(),
                Some(source),
            );
        }
        extract_cycle(&graph, johnson(&graph).unwrap_err(), None);
        extract_cycle(&graph, floyd_warshall(&graph).unwrap_err(), None);
        check_cycle(&graph, &find_negative_cycle(&graph).unwrap(), None);
    }
}

#[test]
fn wide_signed_distances_and_dag_extrema() {
    let wide = graph(
        5,
        GraphKind::Directed,
        [
            (0, 1, i64::MIN),
            (1, 2, i64::MIN),
            (0, 2, i64::MAX),
            (3, 4, i64::MAX),
        ],
    );
    let expected = 2 * i128::from(i64::MIN);
    assert_eq!(bellman_ford(&wide, 0).unwrap().distance(2), Some(expected));
    assert_eq!(johnson(&wide).unwrap().distance(0, 2), Some(expected));
    assert_eq!(
        floyd_warshall(&wide).unwrap().distance(0, 2),
        Some(expected)
    );
    assert_eq!(
        dag_shortest_paths(&wide, 0).unwrap().distance(2),
        Some(expected)
    );
    assert_eq!(
        dag_longest_path(&wide).unwrap().unwrap().cost,
        i128::from(i64::MAX)
    );
    assert!(matches!(
        bellman_ford(&wide, 7),
        Err(SignedPathError::Graph(GraphError::InvalidNode { .. }))
    ));
    let empty = graph(0, GraphKind::Directed, []);
    assert_eq!(johnson(&empty).unwrap().node_count(), 0);
    assert_eq!(floyd_warshall(&empty).unwrap().node_count(), 0);
    assert!(dag_longest_path(&empty).unwrap().is_none());
}

#[test]
fn bellman_ford_handles_a_long_negative_dag_without_repeated_sweeps() {
    let n = 50_000;
    let chain = graph(
        n,
        GraphKind::Directed,
        (1..n).rev().map(|node| (node - 1, node, -1)),
    );
    let paths = bellman_ford(&chain, 0).unwrap();
    let last = paths.path_to(n - 1).unwrap();
    assert_eq!(last.cost, -((n - 1) as i128));
    assert_eq!(last.nodes.len(), n);
    assert_eq!(last.edges.first(), Some(&(n - 2)));
    assert_eq!(last.edges.last(), Some(&0));
}

#[test]
#[ignore = "requires Python and NetworkX"]
fn signed_algorithms_agree_with_networkx() {
    let mut rng = Samples(286_151);
    let random_signed = (0..90).map(|case| {
        let n = 1 + rng.index(19);
        let kind = if case % 5 == 0 {
            GraphKind::Undirected
        } else {
            GraphKind::Directed
        };
        let m = rng.index(n * 5);
        let source = rng.index(n);
        let edges = (0..m).map(|_| (rng.index(n), rng.index(n), rng.index(21) as i64 - 10));
        (graph(n, kind, edges), source)
    });
    let long_scc = [false, true].map(|negative| {
        let n = 81;
        let arcs = (0..n - 1).map(|u| (u, u + 1, -2)).chain([(
            n - 1,
            0,
            2 * (n - 1) as i64 - i64::from(negative),
        )]);
        (graph(n, GraphKind::Directed, arcs), 0)
    });
    let cases: Vec<_> = allowed_cases()
        .into_iter()
        .chain(cycle_cases())
        .chain(random_signed)
        .chain(long_scc)
        .collect();
    let (version, reference) = shortest_oracle::networkx(&cases, 1);
    let oracle_time: std::time::Duration =
        reference.iter().flatten().map(|value| value.elapsed).sum();
    println!(
        "NetworkX {version}, {} graphs and four algorithms in {oracle_time:.2?}",
        cases.len()
    );
    for (case, ((graph, source), reference)) in cases.iter().zip(reference).enumerate() {
        check_duong(graph, *source, &reference[0].distances);
        let results = [
            bellman_ford(graph, *source).map(|paths| paths.distances().collect()),
            johnson(graph).map(|paths| flatten(&paths)),
            floyd_warshall(graph).map(|paths| flatten(&paths)),
        ];
        for (mode, (result, expected)) in results.into_iter().zip(reference).enumerate() {
            match result {
                Ok(distances) => assert_eq!(
                    Some(distances),
                    expected.distances,
                    "case {case} mode {mode}"
                ),
                Err(error) => {
                    assert!(expected.distances.is_none(), "case {case} mode {mode}");
                    extract_cycle(graph, error, (mode == 0).then_some(*source));
                }
            }
        }
    }
}

fn check_duong(graph: &Graph<i64>, source: usize, expected: &Option<Vec<Option<i128>>>) {
    match duong_shortest_paths(graph, source) {
        Ok(paths) => {
            let actual = paths.distances().collect::<Vec<_>>();
            assert_eq!(Some(&actual), expected.as_ref());
            for (target, distance) in actual.into_iter().enumerate() {
                if let Some(distance) = distance {
                    check_path(
                        graph,
                        &paths.path_to(target).unwrap(),
                        source,
                        target,
                        distance,
                    );
                }
            }
        }
        Err(error) => {
            assert!(expected.is_none());
            extract_cycle(graph, error, Some(source));
        }
    }
}
