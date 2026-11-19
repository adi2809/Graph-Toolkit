#[path = "support/weighted_oracle.rs"]
mod weighted_oracle;

use graph_library::matching::{
    WeightedMatching, bipartite_coloring, is_edge_cover, is_matching, is_maximal_matching,
    is_perfect_matching, maximal_matching, maximum_weight_matching,
    minimum_cost_bipartite_matching, minimum_cost_full_bipartite_matching, minimum_cost_matching,
    minimum_cost_perfect_matching, minimum_edge_cover,
};
use graph_library::{Graph, GraphError, GraphKind};

type InputEdge = (usize, usize, i64);

fn graph(n: usize, edges: impl IntoIterator<Item = InputEdge>) -> Graph<i64> {
    let mut builder = Graph::builder(n, GraphKind::Undirected);
    for (a, b, weight) in edges {
        builder.add_edge(a, b, weight).unwrap();
    }
    builder.build()
}

struct Random(u64);

impl Random {
    fn below(&mut self, count: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % count as u64) as usize
    }
}

fn cases(small: usize, large: usize) -> Vec<Graph<i64>> {
    let mut rng = Random(134_737_155);
    let mut cases = vec![
        graph(0, []),
        graph(9, [(0, 0, i64::MAX)]),
        graph(
            6,
            [
                (0, 1, i64::MAX),
                (2, 3, i64::MAX),
                (4, 5, i64::MIN),
                (3, 4, -1),
                (1, 2, i64::MIN + 1),
                (0, 1, i64::MIN),
            ],
        ),
        graph(
            4,
            [(0, 1, -5), (0, 1, -7), (1, 2, 99), (2, 3, -4), (3, 3, 500)],
        ),
    ];
    for index in 0..small + large {
        let n = if index < small {
            1 + rng.below(11)
        } else {
            12 + rng.below(64)
        };
        let count = rng.below(if index % 2 == 0 { n * n } else { n * 3 });
        let weight_table = [i64::MIN, i64::MAX, -7, -1, 0, 1, 5, 8, 13];
        let edges = (0..count).map(|_| {
            let weight = if index % 7 == 0 {
                weight_table[rng.below(weight_table.len())]
            } else {
                rng.below(101) as i64 - 50
            };
            (rng.below(n), rng.below(n), weight)
        });
        cases.push(graph(n, edges));
    }
    cases
}

fn results(graph: &Graph<i64>) -> [(usize, i128); 3] {
    [
        maximum_weight_matching(graph, false).unwrap(),
        maximum_weight_matching(graph, true).unwrap(),
        minimum_cost_matching(graph).unwrap(),
    ]
    .map(|result| validate(graph, &result))
}

fn validate(graph: &Graph<i64>, result: &WeightedMatching) -> (usize, i128) {
    let mut count = 0;
    let mut weight = 0;
    let mut seen = vec![false; graph.node_count()];
    for id in result.matching.edges() {
        let edge = graph.edge(id).unwrap();
        assert_ne!(edge.from, edge.to);
        assert!(!seen[edge.from] && !seen[edge.to]);
        seen[edge.from] = true;
        seen[edge.to] = true;
        assert_eq!(result.matching.mate(edge.from).unwrap().edge, id);
        assert_eq!(result.matching.mate(edge.to).unwrap().edge, id);
        assert_eq!(result.matching.mate(edge.from).unwrap().node, edge.to);
        assert_eq!(result.matching.mate(edge.to).unwrap().node, edge.from);
        weight += i128::from(edge.weight);
        count += 1;
    }
    for (vertex, seen) in seen.into_iter().enumerate() {
        assert_eq!(seen, result.matching.mate(vertex).is_some());
    }
    assert_eq!(count, result.matching.len());
    assert_eq!(weight, result.total_weight);
    (count, weight)
}

fn exhaustive(graph: &Graph<i64>, mode: usize) -> (usize, i128) {
    fn search(
        graph: &Graph<i64>,
        vertices: usize,
        mode: usize,
        memo: &mut [Option<(usize, i128)>],
    ) -> (usize, i128) {
        if vertices == 0 {
            return (0, 0);
        }
        if let Some(answer) = memo[vertices] {
            return answer;
        }
        let vertex = vertices.trailing_zeros() as usize;
        let rest = vertices & !(1 << vertex);
        let mut answer = search(graph, rest, mode, memo);
        for arc in graph.outgoing(vertex).unwrap() {
            if rest & (1 << arc.node) == 0 {
                continue;
            }
            let (count, weight) = search(graph, rest & !(1 << arc.node), mode, memo);
            let next = (
                count + 1,
                weight + i128::from(graph.edges()[arc.edge].weight),
            );
            let improvement = match mode {
                0 => next.1 > answer.1,
                1 => next > answer,
                _ => (next.0, -next.1) > (answer.0, -answer.1),
            };
            if improvement {
                answer = next;
            }
        }
        memo[vertices] = Some(answer);
        answer
    }
    let limit = 1 << graph.node_count();
    search(graph, limit - 1, mode, &mut vec![None; limit])
}

#[test]
fn weighted_matching_agrees_with_exhaustive_optima() {
    for (case, graph) in cases(500, 0).iter().enumerate() {
        let answers = results(graph);
        for (mode, &actual) in answers.iter().enumerate() {
            let expected = exhaustive(graph, mode);
            assert_eq!(actual.1, expected.1, "case {case}, mode {mode}");
            if mode != 0 {
                assert_eq!(actual.0, expected.0, "case {case}, mode {mode}");
            }
        }
        let perfect = minimum_cost_perfect_matching(graph).unwrap();
        let expected = answers[2];
        assert_eq!(perfect.is_some(), expected.0 * 2 == graph.node_count());
        if let Some(perfect) = perfect {
            assert_eq!(validate(graph, &perfect), expected);
        }
    }
}

#[test]
#[ignore = "requires Python and NetworkX"]
fn weighted_matching_agrees_with_networkx() {
    let graphs = cases(450, 160);
    let (version, expected) = weighted_oracle::networkx(&graphs, 1);
    let oracle_time: std::time::Duration =
        expected.iter().flatten().map(|entry| entry.elapsed).sum();
    println!(
        "NetworkX {version}, {} graphs and three objectives in {oracle_time:.2?}",
        graphs.len()
    );
    for (case, (graph, expected)) in graphs.iter().zip(expected).enumerate() {
        for (mode, (actual, expected)) in results(graph).into_iter().zip(expected).enumerate() {
            assert_eq!(actual.1, expected.weight, "case {case}, mode {mode}");
            if mode != 0 {
                assert_eq!(actual.0, expected.size, "case {case}, mode {mode}");
            }
        }
    }
}

#[test]
fn sparse_rectangular_assignment_agrees_with_exhaustive_optima() {
    let mut rng = Random(221_163);
    for case in 0..350 {
        let (a, b) = (1 + rng.below(5), 1 + rng.below(6));
        let size = rng.below(a * b * 2);
        let graph = graph(
            a + b,
            (0..size).map(|_| {
                let u = rng.below(a);
                let v = a + rng.below(b);
                let weight = [i64::MIN, -5, -1, 0, 1, 7, i64::MAX][rng.below(7)];
                if rng.below(2) == 0 {
                    (u, v, weight)
                } else {
                    (v, u, weight)
                }
            }),
        );
        let left: Vec<_> = (0..a).collect();
        let right: Vec<_> = (a..a + b).collect();
        let expected = exhaustive(&graph, 2);
        for side in [&left, &right] {
            let actual = minimum_cost_bipartite_matching(&graph, side).unwrap();
            assert_eq!(validate(&graph, &actual), expected, "case {case}");
            let full = minimum_cost_full_bipartite_matching(&graph, side).unwrap();
            assert_eq!(full.is_some(), expected.0 == a.min(b));
            if let Some(full) = full {
                assert_eq!(validate(&graph, &full), expected);
            }
        }
        let colors = bipartite_coloring(&graph).unwrap();
        assert!(
            graph
                .edges()
                .iter()
                .all(|edge| colors[edge.from] != colors[edge.to])
        );
    }
    let triangle = graph(3, [(0, 1, 0), (1, 2, 3), (2, 0, 4)]);
    assert!(matches!(
        bipartite_coloring(&triangle),
        Err(GraphError::InvalidBipartition { .. })
    ));
    assert_eq!(
        minimum_cost_bipartite_matching(&triangle, &[0]).unwrap_err(),
        GraphError::InvalidBipartition { edge: 1 }
    );
}

#[test]
fn matching_checks_and_minimum_edge_cover_use_original_ids() {
    let path = graph(4, [(1, 2, 10), (0, 1, 3), (2, 3, 1), (3, 3, 7)]);
    let greedy = maximal_matching(&path).unwrap();
    let chosen: Vec<_> = greedy.edges().collect();
    assert!(is_matching(&path, &chosen).unwrap());
    assert!(is_maximal_matching(&path, &chosen).unwrap());
    assert!(is_perfect_matching(&path, &[1, 2]).unwrap());
    assert!(is_maximal_matching(&path, &[0]).unwrap());
    assert!(!is_perfect_matching(&path, &[0]).unwrap());
    assert!(!is_maximal_matching(&path, &[]).unwrap());
    assert!(!is_matching(&path, &[1, 1]).unwrap());
    assert!(!is_matching(&path, &[3]).unwrap());
    assert!(!is_matching(&path, &[4]).unwrap());
    assert!(is_edge_cover(&path, &[1, 2, 2]));
    assert!(!is_edge_cover(&path, &[1, 3]));

    let mut rng = Random(116_507);
    for _ in 0..120 {
        let n = 1 + rng.below(8);
        let size = rng.below(12);
        let graph = graph(n, (0..size).map(|_| (rng.below(n), rng.below(n), 1)));
        let actual = minimum_edge_cover(&graph).unwrap();
        let mut expected = None;
        for mask in 0..1_usize << size {
            let mut vertices = 0;
            for (id, edge) in graph.edges().iter().enumerate() {
                if mask & (1 << id) != 0 {
                    vertices |= (1 << edge.from) | (1 << edge.to);
                }
            }
            if vertices == (1 << n) - 1 {
                let count = mask.count_ones() as usize;
                expected = Some(expected.map_or(count, |old: usize| old.min(count)));
            }
        }
        assert_eq!(actual.as_ref().map(Vec::len), expected);
        if let Some(actual) = actual {
            assert!(is_edge_cover(&graph, &actual));
        }
    }
}
