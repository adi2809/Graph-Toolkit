#[path = "support/structure_oracle.rs"]
mod oracle;

use std::collections::{BTreeSet, HashSet};
use std::fmt::Write as _;

use graph_library::functions::induced_subgraph;
use graph_library::generators::{
    complete_bipartite_graph, complete_graph, cycle_graph, empty_graph, gnp_random_graph,
    path_graph, tutte_graph,
};
use graph_library::isomorphism::{
    vf2pp_all_subgraph_isomorphisms, vf2pp_all_subgraph_isomorphisms_with_labels,
    vf2pp_subgraph_is_isomorphic, vf2pp_subgraph_isomorphism,
    vf2pp_subgraph_isomorphism_with_labels,
};
use graph_library::planarity::{
    KuratowskiKind, KuratowskiWitness, check_planarity_with_witness, is_planar, kuratowski_witness,
};
use graph_library::{Graph, GraphKind};

const UNDIRECTED: GraphKind = GraphKind::Undirected;

fn graph<W>(
    n: usize,
    kind: GraphKind,
    edges: impl IntoIterator<Item = (usize, usize, W)>,
) -> Graph<W> {
    let mut builder = Graph::builder(n, kind);
    for (u, v, w) in edges {
        builder.add_edge(u, v, w).unwrap();
    }
    builder.build()
}

fn verify_witness<W>(graph: &Graph<W>, witness: &KuratowskiWitness) {
    let pairs: BTreeSet<_> = match witness.kind {
        KuratowskiKind::K5 => {
            assert_eq!(witness.branch_vertices.len(), 5);
            (0..5)
                .flat_map(|i| (i + 1..5).map(move |j| (i, j)))
                .collect()
        }
        KuratowskiKind::K33 => {
            assert_eq!(witness.branch_vertices.len(), 6);
            (0..3).flat_map(|i| (3..6).map(move |j| (i, j))).collect()
        }
    };
    let mut interiors = BTreeSet::new();
    let mut used_edges = BTreeSet::new();
    let mut represented = BTreeSet::new();
    let branch = &witness.branch_vertices;
    assert_eq!(branch.iter().collect::<BTreeSet<_>>().len(), branch.len());
    for path in &witness.paths {
        assert!(!path.edges.is_empty());
        assert_eq!(path.edges.len() + 1, path.nodes.len());
        let u = branch.iter().position(|v| v == &path.nodes[0]).unwrap();
        let v = branch
            .iter()
            .position(|v| v == path.nodes.last().unwrap())
            .unwrap();
        assert!(represented.insert((u.min(v), u.max(v))));
        for &interior in &path.nodes[1..path.nodes.len() - 1] {
            assert!(!branch.contains(&interior));
            assert!(interiors.insert(interior));
        }
        for (&edge, pair) in path.edges.iter().zip(path.nodes.windows(2)) {
            let endpoints = graph.edge(edge).unwrap();
            assert!(used_edges.insert(edge));
            assert!(
                (pair[0] == endpoints.from && pair[1] == endpoints.to)
                    || (pair[0] == endpoints.to && pair[1] == endpoints.from)
            );
        }
    }
    assert_eq!(pairs, represented);
    assert_eq!(used_edges, witness.edge_ids().collect());
}

fn subdivision(base: &Graph<()>, length: usize) -> Graph<()> {
    let n = base.node_count() + (length - 1) * base.edge_count();
    let mut builder = Graph::builder(n, UNDIRECTED);
    let mut next = base.node_count();
    for item in base.edges() {
        let mut previous = item.from;
        for _ in 1..length {
            builder.add_edge(previous, next, ()).unwrap();
            previous = next;
            next += 1;
        }
        builder.add_edge(previous, item.to, ()).unwrap();
    }
    builder.build()
}

#[test]
fn certificate_lifts_long_chains_to_original_edge_ids() {
    for (base, kind, length) in [
        (
            complete_graph(5, UNDIRECTED).unwrap(),
            KuratowskiKind::K5,
            1,
        ),
        (
            complete_bipartite_graph(3, 3).unwrap(),
            KuratowskiKind::K33,
            4096,
        ),
    ] {
        let obstruction = subdivision(&base, length);
        let graph = graph(
            obstruction.node_count() + 5,
            GraphKind::Directed,
            [(0, 0, ()), (1, 3, ()), (3, 2, ()), (2, 1, ()), (1, 5, ())]
                .into_iter()
                .chain(obstruction.edges().iter().flat_map(|edge| {
                    [
                        (edge.to + 5, edge.from + 5, ()),
                        (edge.from + 5, edge.to + 5, ()),
                    ]
                })),
        );
        let result = kuratowski_witness(&graph).unwrap();
        assert_eq!(result.kind, kind);
        assert_eq!(result.edge_ids().count(), base.edge_count() * length);
        assert!(result.edge_ids().all(|id| id % 2 == 1));
        verify_witness(&graph, &result);
        assert_eq!(check_planarity_with_witness(&graph).unwrap_err(), result);
    }
    let planar = tutte_graph();
    assert!(kuratowski_witness(&planar).is_none());
    assert_eq!(
        check_planarity_with_witness(&planar).unwrap().edge_count(),
        planar.edge_count()
    );
}

fn multiplicity<W>(graph: &Graph<W>) -> Vec<Vec<usize>> {
    let mut matrix = vec![vec![0; graph.node_count()]; graph.node_count()];
    for edge in graph.edges() {
        matrix[edge.from][edge.to] += 1;
        if !graph.is_directed() && edge.from != edge.to {
            matrix[edge.to][edge.from] += 1;
        }
    }
    matrix
}

fn verify_mapping<A, B>(host: &Graph<A>, pattern: &Graph<B>, mapping: &[usize]) {
    let a = multiplicity(host);
    let b = multiplicity(pattern);
    assert_eq!(mapping.len(), pattern.node_count());
    assert_eq!(mapping.iter().collect::<HashSet<_>>().len(), mapping.len());
    for u in pattern.nodes() {
        for v in pattern.nodes() {
            assert_eq!(a[mapping[u]][mapping[v]], b[u][v]);
        }
    }
}

fn exhaustive<A, B>(
    host: &Graph<A>,
    pattern: &Graph<B>,
    labels: Option<(&[usize], &[usize])>,
) -> BTreeSet<Vec<usize>> {
    let a = multiplicity(host);
    let b = multiplicity(pattern);
    let mut found = BTreeSet::new();
    let mut pending = vec![Vec::new()];
    while let Some(path) = pending.pop() {
        let u = path.len();
        if u == pattern.node_count() {
            found.insert(path);
            continue;
        }
        for v in host.nodes().filter(|v| !path.contains(v)) {
            if labels.is_some_and(|(x, y)| x[v] != y[u]) || a[v][v] != b[u][u] {
                continue;
            }
            if path
                .iter()
                .enumerate()
                .all(|(i, &j)| a[v][j] == b[u][i] && a[j][v] == b[i][u])
            {
                let mut next = path.clone();
                next.push(v);
                pending.push(next);
            }
        }
    }
    found
}

#[test]
fn induced_rules_preserve_nonedges_direction_loops_and_multiplicity() {
    let host = complete_graph(7, UNDIRECTED).unwrap();
    assert!(!vf2pp_subgraph_is_isomorphic(
        &host,
        &path_graph(3, UNDIRECTED)
    ));
    assert_eq!(
        vf2pp_all_subgraph_isomorphisms(&host, &complete_graph(3, UNDIRECTED).unwrap()).count(),
        210
    );
    let ring = cycle_graph(5, UNDIRECTED);
    assert_eq!(
        vf2pp_all_subgraph_isomorphisms(&ring, &path_graph(3, UNDIRECTED)).count(),
        10
    );
    let unmatched = empty_graph(3, UNDIRECTED);
    assert!(!vf2pp_subgraph_is_isomorphic(&ring, &unmatched));
    let empty = empty_graph(0, UNDIRECTED);
    assert_eq!(
        vf2pp_all_subgraph_isomorphisms(&host, &empty).collect::<Vec<_>>(),
        [Vec::<usize>::new()]
    );

    let host = graph(
        5,
        GraphKind::Directed,
        [
            (1, 1, ()),
            (1, 3, ()),
            (1, 3, ()),
            (1, 4, ()),
            (0, 1, ()),
            (4, 3, ()),
            (0, 0, ()),
            (3, 0, ()),
        ],
    );
    let pattern = graph(2, GraphKind::Directed, [(0, 0, ()), (0, 1, ()), (0, 1, ())]);
    let mapping =
        vf2pp_subgraph_isomorphism_with_labels(&host, &pattern, &[0, 5, 0, 9, 8], &[5, 9])
            .unwrap()
            .unwrap();
    assert_eq!(mapping, [1, 3]);
    let reverse = graph(2, GraphKind::Directed, [(0, 0, ()), (1, 0, ()), (1, 0, ())]);
    assert!(!vf2pp_subgraph_is_isomorphic(&host, &reverse));
    assert!(vf2pp_subgraph_isomorphism_with_labels(&host, &pattern, &[5], &[5, 9]).is_err());
}

fn induced_pair(trial: usize) -> (Graph<()>, Graph<()>, Vec<usize>, Vec<usize>) {
    let n = 3 + trial % 15;
    let k = (trial * 19 + trial / 5) % (n + 1);
    let kind = if trial % 2 == 0 {
        UNDIRECTED
    } else {
        GraphKind::Directed
    };
    let p = 0.08 * (trial % 11 + 1) as f64;
    let base = gnp_random_graph(n, p, kind, 10_008 + trial as u64).unwrap();
    let host = graph(
        n,
        kind,
        base.edges()
            .iter()
            .map(|e| (e.from, e.to, ()))
            .chain((trial % 3 == 0).then_some((0, 0, ())))
            .chain(
                (trial % 4 == 0)
                    .then(|| base.edges().first().map(|e| (e.from, e.to, ())))
                    .flatten(),
            ),
    );
    let positions: Vec<_> = (0..k).map(|u| (u + trial / 7) % n).collect();
    let mut local = vec![usize::MAX; n];
    for (i, &u) in positions.iter().enumerate() {
        local[u] = i;
    }
    let pattern = if trial % 3 == 1 {
        gnp_random_graph(k, p, kind, 24_056 + trial as u64).unwrap()
    } else {
        graph(
            k,
            kind,
            host.edges()
                .iter()
                .filter(|e| local[e.from] != usize::MAX && local[e.to] != usize::MAX)
                .map(|e| (local[e.from], local[e.to], ())),
        )
    };
    let host_labels: Vec<_> = host.nodes().map(|u| (u + trial / 3) % 3).collect();
    let mut pattern_labels: Vec<_> = positions.iter().map(|&u| host_labels[u]).collect();
    if k != 0 && trial % 13 == 0 {
        pattern_labels[0] = (pattern_labels[0] + 1) % 3;
    }
    (host, pattern, host_labels, pattern_labels)
}

#[test]
fn induced_iterator_agrees_with_exhaustive_injections() {
    for trial in (0..165).filter(|trial| trial % 15 < 5) {
        let (host, pattern, host_labels, pattern_labels) = induced_pair(trial);
        let expected = exhaustive(&host, &pattern, None);
        let actual: Vec<_> = vf2pp_all_subgraph_isomorphisms(&host, &pattern).collect();
        assert_eq!(actual.len(), expected.len(), "unlabelled case {trial}");
        assert_eq!(
            actual.into_iter().collect::<BTreeSet<_>>(),
            expected,
            "unlabelled case {trial}"
        );
        let expected = exhaustive(&host, &pattern, Some((&host_labels, &pattern_labels)));
        let actual: Vec<_> = vf2pp_all_subgraph_isomorphisms_with_labels(
            &host,
            &pattern,
            &host_labels,
            &pattern_labels,
        )
        .unwrap()
        .collect();
        assert_eq!(actual.len(), expected.len(), "labelled case {trial}");
        assert_eq!(
            actual.into_iter().collect::<BTreeSet<_>>(),
            expected,
            "labelled case {trial}"
        );
    }
    let host = subdivision(&path_graph(250, UNDIRECTED), 50);
    let pattern = induced_subgraph(&host, 0..500).unwrap().materialize().graph;
    let mapping = vf2pp_subgraph_isomorphism(&host, &pattern).unwrap();
    assert_eq!(mapping.len(), pattern.node_count());
}

#[test]
#[ignore = "requires Python and NetworkX"]
fn kuratowski_subdivisions_agree_with_networkx() {
    let (version, atlas) = oracle::invoke("atlas", "");
    let mut cases = Vec::new();
    for row in atlas {
        let n = row[0].parse().unwrap();
        let edges = row[2..]
            .chunks_exact(2)
            .map(|chunk| (chunk[0].parse().unwrap(), chunk[1].parse().unwrap(), ()));
        cases.push(graph(n, UNDIRECTED, edges));
    }
    for trial in 0..110 {
        let n = 5 + trial % 19;
        let p = 0.04 * (trial % 13 + 1) as f64;
        let kind = if trial % 2 == 0 {
            UNDIRECTED
        } else {
            GraphKind::Directed
        };
        let base = gnp_random_graph(n, p, kind, 10_791 + trial as u64).unwrap();
        cases.push(graph(
            n,
            kind,
            base.edges()
                .iter()
                .rev()
                .map(|e| (e.from, e.to, ()))
                .chain([(0, 0, ()), (1, 0, ()), (0, 1, ())]),
        ));
    }
    let mut query = format!("{}\n", cases.len());
    let mut expected = Vec::new();
    for graph in &cases {
        let witness = kuratowski_witness(graph);
        assert_eq!(witness.is_none(), is_planar(graph));
        if let Some(ref witness) = witness {
            verify_witness(graph, witness);
        }
        expected.push(witness.is_none());
        query.push_str("C\n");
        oracle::graph(&mut query, graph, |_| 1.0);
        let edges: Vec<_> = witness.iter().flat_map(|w| w.edge_ids()).collect();
        writeln!(query, "{}", edges.len()).unwrap();
        for id in edges {
            let edge = graph.edges()[id];
            writeln!(query, "{} {}", edge.from, edge.to).unwrap();
        }
    }
    let (_, output) = oracle::invoke("1", &query);
    assert_eq!(expected.len(), output.len());
    for (id, (planar, result)) in expected.into_iter().zip(output).enumerate() {
        assert_eq!(planar, result[1] == "1", "case {id}");
    }
    println!(
        "NetworkX {version}, {} planarity decisions and independent subdivision checks",
        cases.len()
    );
}

#[test]
#[ignore = "requires Python and NetworkX"]
fn induced_vf2pp_agrees_with_networkx() {
    let n = 300;
    let mut query = format!("{n}\n");
    let mut expected = Vec::new();
    for trial in 0..n {
        let (host, pattern, host_labels, pattern_labels) = induced_pair(trial);
        let labels = trial % 3 != 2;
        let result = if labels {
            vf2pp_subgraph_isomorphism_with_labels(&host, &pattern, &host_labels, &pattern_labels)
                .unwrap()
        } else {
            vf2pp_subgraph_isomorphism(&host, &pattern)
        };
        if let Some(ref mapping) = result {
            verify_mapping(&host, &pattern, mapping);
            if labels {
                assert!(
                    mapping
                        .iter()
                        .enumerate()
                        .all(|(p, &h)| pattern_labels[p] == host_labels[h])
                );
            }
        }
        expected.push(result.is_some());
        query.push_str("S\n");
        oracle::graph(&mut query, &host, |_| 1.0);
        oracle::graph(&mut query, &pattern, |_| 1.0);
        writeln!(query, "{}", u8::from(labels)).unwrap();
        if labels {
            for label in host_labels.into_iter().chain(pattern_labels) {
                write!(query, "{label} ").unwrap();
            }
            query.push('\n');
        }
    }
    let (version, output) = oracle::invoke("1", &query);
    assert_eq!(expected.len(), output.len());
    for (id, (answer, nx)) in expected.into_iter().zip(output).enumerate() {
        assert_eq!(answer, nx[1] == "1", "case {id}");
    }
    println!("NetworkX {version}, {n} induced subgraph comparisons");
}
