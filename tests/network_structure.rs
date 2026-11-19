#[path = "support/structure_oracle.rs"]
mod oracle;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use graph_library::community::{
    CommunityError, LouvainOptions, louvain_communities, louvain_communities_by, modularity,
};
use graph_library::connectivity::weakly_connected_components;
use graph_library::generators::*;
use graph_library::isomorphism::{
    vf2pp_all_isomorphisms, vf2pp_all_isomorphisms_with_labels, vf2pp_is_isomorphic,
    vf2pp_isomorphism, vf2pp_isomorphism_with_labels,
};
use graph_library::planarity::{check_planarity, is_planar};
use graph_library::{Graph, GraphKind};

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

fn project<W>(graph: &Graph<W>, map: &[usize]) -> BTreeMap<(usize, usize), usize> {
    let mut edges = BTreeMap::new();
    for edge in graph.edges() {
        let (u, v) = (map[edge.from], map[edge.to]);
        let key = if graph.is_directed() || u <= v {
            (u, v)
        } else {
            (v, u)
        };
        *edges.entry(key).or_default() += 1;
    }
    edges
}

fn verify_bijection<A, B>(a: &Graph<A>, b: &Graph<B>, mapping: &[usize]) {
    assert_eq!(
        mapping.iter().copied().collect::<BTreeSet<_>>(),
        b.nodes().collect()
    );
    assert_eq!(
        project(a, mapping),
        project(b, &b.nodes().collect::<Vec<_>>())
    );
}

fn verify_rotation<W>(graph: &Graph<W>, planar: bool) {
    assert_eq!(is_planar(graph), planar);
    let embedding = check_planarity(graph);
    assert_eq!(embedding.is_some(), planar);
    let Some(embedding) = embedding else { return };
    assert_eq!(embedding.node_count(), graph.node_count());
    let mut stable = BTreeMap::new();
    for (id, edge) in graph.edges().iter().enumerate() {
        if edge.from != edge.to {
            let u = edge.from.min(edge.to);
            let v = edge.from.max(edge.to);
            stable.entry((u, v)).or_insert(id);
        }
    }
    let mut seen = BTreeSet::new();
    let mut count = 0;
    for u in graph.nodes() {
        for arc in embedding.clockwise_neighbors(u).unwrap() {
            assert_eq!(stable[&(u.min(arc.node), u.max(arc.node))], arc.edge);
            assert!(seen.insert((u, arc.node)));
            count += 1;
        }
    }
    assert_eq!(count, 2 * stable.len());
    assert_eq!(embedding.edge_count(), stable.len());
    let components = weakly_connected_components(graph);
    let mut vertices = vec![0; components.count];
    let mut edges = vec![0; components.count];
    let mut faces = vec![0; components.count];
    for id in &components.component_of {
        vertices[*id] += 1;
    }
    for &(a, _) in stable.keys() {
        edges[components.component_of[a]] += 1;
    }
    let mut visited = BTreeSet::new();
    for face in embedding.face_boundaries() {
        faces[components.component_of[face[0]]] += 1;
        for (&a, &b) in face.iter().zip(face.iter().cycle().skip(1)) {
            assert!(seen.contains(&(a, b)));
            assert!(visited.insert((a, b)));
        }
    }
    assert_eq!(visited, seen);
    for component in 0..components.count {
        if edges[component] != 0 {
            assert_eq!(vertices[component] + faces[component], edges[component] + 2);
        }
    }
}

fn subdivided_biclique(length: usize) -> Graph<()> {
    let mut n = 6;
    let mut edges = Vec::new();
    for a in 0..3 {
        for b in 3..6 {
            let mut u = a;
            for _ in 1..length {
                edges.push((u, n, ()));
                u = n;
                n += 1;
            }
            edges.push((u, b, ()));
        }
    }
    graph(n, GraphKind::Undirected, edges)
}

#[test]
fn planar_embeddings_cover_disconnected_and_deep_inputs() {
    verify_rotation(&null_graph(), true);
    verify_rotation(&complete_graph(4, GraphKind::Undirected).unwrap(), true);
    verify_rotation(&complete_graph(5, GraphKind::Undirected).unwrap(), false);
    verify_rotation(&subdivided_biclique(100), false);
    verify_rotation(&tutte_graph(), true);
    verify_rotation(&icosahedral_graph(), true);
    verify_rotation(&heawood_graph(), false);
    verify_rotation(&hexagonal_lattice_graph(15, 17, false).unwrap(), true);
    verify_rotation(&hexagonal_lattice_graph(6, 6, true).unwrap(), false);
    verify_rotation(&path_graph(45_000, GraphKind::Undirected), true);
    let mixed = graph(
        8,
        GraphKind::Directed,
        [
            (0, 1, ()),
            (1, 0, ()),
            (0, 1, ()),
            (1, 2, ()),
            (2, 0, ()),
            (0, 3, ()),
            (3, 3, ()),
            (4, 5, ()),
        ],
    );
    verify_rotation(&mixed, true);
}

#[test]
fn vf2pp_handles_symmetry_labels_and_deep_search() {
    let cycle = cycle_graph(7, GraphKind::Undirected);
    let solutions: Vec<_> = vf2pp_all_isomorphisms(&cycle, &cycle).collect();
    assert_eq!(solutions.len(), 14);
    assert_eq!(solutions.iter().collect::<BTreeSet<_>>().len(), 14);
    for solution in solutions {
        verify_bijection(&cycle, &cycle, &solution);
    }
    let complete = complete_graph(5, GraphKind::Undirected).unwrap();
    let permutations: BTreeSet<_> = vf2pp_all_isomorphisms(&complete, &complete).collect();
    assert_eq!(permutations.len(), 120);
    for mapping in &permutations {
        verify_bijection(&complete, &complete, mapping);
    }
    let labels = [0, 1, 1, 1, 1, 1, 1];
    assert_eq!(
        vf2pp_all_isomorphisms_with_labels(&cycle, &cycle, &labels, &labels)
            .unwrap()
            .count(),
        2
    );
    let path = path_graph(25_000, GraphKind::Directed);
    assert_eq!(
        vf2pp_isomorphism(&path, &path).unwrap(),
        path.nodes().collect::<Vec<_>>()
    );
    let isolates = empty_graph(30_000, GraphKind::Undirected);
    assert!(vf2pp_is_isomorphic(&isolates, &isolates));
    assert_eq!(
        vf2pp_all_isomorphisms(&null_graph(), &null_graph()).collect::<Vec<_>>(),
        [Vec::<usize>::new()]
    );
    assert!(!vf2pp_is_isomorphic(
        &petersen_graph(),
        &circular_ladder_graph(5).unwrap()
    ));
}

#[test]
fn weighted_louvain_handles_loops_parallel_edges_and_resolution() {
    for kind in [GraphKind::Undirected, GraphKind::Directed] {
        let mut edges = vec![
            (0, 0, 4.0),
            (3, 3, 4.0),
            (2, 3, 0.1),
            (0, 2, 3.0),
            (0, 2, 4.0),
        ];
        for i in [0, 3] {
            for u in i..i + 3 {
                for v in u + 1..i + 3 {
                    edges.push((u, v, 10.0));
                    if kind == GraphKind::Directed {
                        edges.push((v, u, 10.0));
                    }
                }
            }
        }
        let graph = graph(7, kind, edges);
        let result = louvain_communities(
            &graph,
            LouvainOptions {
                seed: 533,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            result.communities(),
            [vec![0, 1, 2], vec![3, 4, 5], vec![6]]
        );
        assert!(
            (modularity(&graph, &result.community_of, 1.0).unwrap() - result.modularity).abs()
                < 1e-12
        );
        assert!(result.level_modularities.is_sorted());
        let high = LouvainOptions {
            resolution: 100.0,
            ..Default::default()
        };
        assert_eq!(
            louvain_communities(&graph, high).unwrap().community_count(),
            7
        );
    }
    let empty = empty_graph(4, GraphKind::Directed);
    assert_eq!(
        louvain_communities_by(&empty, |_| 1.0, LouvainOptions::default())
            .unwrap()
            .community_of,
        [0, 1, 2, 3]
    );
    let invalid = graph(2, GraphKind::Undirected, [(0, 1, -1.0)]);
    assert_eq!(
        louvain_communities(&invalid, LouvainOptions::default()).unwrap_err(),
        CommunityError::InvalidWeight { edge: 0 }
    );
    let enormous = graph(3, GraphKind::Directed, [(0, 1, f64::MAX), (1, 2, f64::MAX)]);
    assert!(
        louvain_communities(&enormous, LouvainOptions::default())
            .unwrap()
            .modularity
            .is_finite()
    );
}

#[test]
fn added_generators_preserve_degrees_and_copy_structure() {
    let base = graph(
        3,
        GraphKind::Directed,
        [(0, 1, 8), (0, 1, 12), (1, 2, 4), (2, 2, 3)],
    );
    let k = 9;
    let lifted = random_k_lift(&base, k, 222).unwrap();
    for (id, original) in base.edges().iter().enumerate() {
        let block = &lifted.edges()[id * k..(id + 1) * k];
        let mut destinations = BTreeSet::new();
        for (copy, edge) in block.iter().enumerate() {
            assert_eq!(edge.from, original.from * k + copy);
            assert_eq!(edge.to / k, original.to);
            assert_eq!(edge.weight, original.weight);
            assert!(destinations.insert(edge.to % k));
        }
    }
    assert_eq!(
        lifted.edges(),
        random_k_lift(&base, k, 222).unwrap().edges()
    );
    let config = directed_configuration_model(&[3, 0, 1, 4], &[0, 5, 1, 2], 123).unwrap();
    assert_eq!(
        config
            .nodes()
            .map(|u| config.incoming(u).unwrap().len())
            .collect::<Vec<_>>(),
        [3, 0, 1, 4]
    );
    assert_eq!(
        config
            .nodes()
            .map(|u| config.outgoing(u).unwrap().len())
            .collect::<Vec<_>>(),
        [0, 5, 1, 2]
    );
    assert_eq!(
        margulis_gabber_galil_graph(7)
            .unwrap()
            .degree_histogram()
            .last(),
        Some(&49)
    );
    let periodic = triangular_lattice_graph(12, 30, true).unwrap();
    assert_eq!(periodic.degree_histogram().last(), Some(&180));
    assert_eq!(periodic.degree(0).unwrap(), 6);
    assert!(triangular_lattice_graph(2, 5, true).is_err());
    assert!(hexagonal_lattice_graph(3, 7, true).is_err());
}

#[test]
#[ignore = "requires Python and NetworkX"]
fn left_right_agrees_with_networkx_graph_atlas() {
    let (version, atlas) = oracle::invoke("atlas", "");
    println!(
        "NetworkX {version}, {} graph atlas planarity checks",
        atlas.len()
    );
    assert_eq!(atlas.len(), 1253);
    for fields in atlas {
        let n = fields[0].parse().unwrap();
        let expected = fields[1] == "1";
        let edges = fields[2..]
            .chunks_exact(2)
            .map(|pair| (pair[0].parse().unwrap(), pair[1].parse().unwrap(), ()));
        verify_rotation(&graph(n, GraphKind::Undirected, edges), expected);
    }
}

#[test]
#[ignore = "requires Python and NetworkX"]
fn graph_structure_agrees_with_networkx() {
    let mut cases = Vec::new();
    let mut input = String::new();
    let mut cases_modularity = Vec::new();
    let mut quality = Vec::new();
    for trial in 0..210 {
        let n = 2 + trial % 19;
        let kind = if trial % 2 == 0 {
            GraphKind::Undirected
        } else {
            GraphKind::Directed
        };
        let base =
            gnp_random_graph(n, 0.05 * (1 + trial % 13) as f64, kind, 9209 + trial as u64).unwrap();
        let extra = (trial % 4 == 0).then_some((0, 0, ()));
        let parallel = (trial % 3 == 0)
            .then(|| base.edges().first().map(|e| (e.from, e.to, ())))
            .flatten();
        let left = graph(
            n,
            kind,
            base.edges()
                .iter()
                .map(|e| (e.from, e.to, ()))
                .chain(extra)
                .chain(parallel),
        );
        let shift = trial % n;
        let perm = |u| (n - 1 - u + shift) % n;
        let right = graph(
            n,
            kind,
            left.edges().iter().rev().map(|e| {
                let to = if trial % 5 == 0 && e.from == 0 {
                    (e.to + 1) % n
                } else {
                    e.to
                };
                (perm(e.from), perm(to), ())
            }),
        );
        let labels: Vec<_> = (0..n).map(|u| (u * (u % 3 + 1)) % 5).collect();
        let mut right_labels = vec![0; n];
        for u in 0..n {
            right_labels[perm(u)] = labels[u];
        }
        if trial % 7 == 0 {
            right_labels[0] = (right_labels[0] + 1) % 5;
        }
        let with_labels = trial % 3 != 0;
        let result = if with_labels {
            vf2pp_isomorphism_with_labels(&left, &right, &labels, &right_labels).unwrap()
        } else {
            vf2pp_isomorphism(&left, &right)
        };
        if let Some(ref mapping) = result {
            verify_bijection(&left, &right, mapping);
            if with_labels {
                assert!(
                    mapping
                        .iter()
                        .enumerate()
                        .all(|(u, &v)| labels[u] == right_labels[v])
                );
            }
        }
        input.push_str("I\n");
        oracle::graph(&mut input, &left, |_| 1.0);
        oracle::graph(&mut input, &right, |_| 1.0);
        writeln!(input, "{}", u8::from(with_labels)).unwrap();
        if with_labels {
            for label in labels.into_iter().chain(right_labels) {
                write!(input, "{label} ").unwrap();
            }
            input.push('\n');
        }
        cases.push(result.is_some());

        if trial % 4 == 0 {
            let graph = graph(
                n,
                kind,
                left.edges()
                    .iter()
                    .enumerate()
                    .map(|(id, e)| (e.from, e.to, (id % 8 + 1) as f64)),
            );
            let options = LouvainOptions {
                seed: trial as u64,
                resolution: 0.5 + (trial % 7) as f64 / 4.0,
                max_levels: (trial % 3 == 0).then_some(1),
                ..Default::default()
            };
            let result = louvain_communities(&graph, options).unwrap();
            let standalone = modularity(&graph, &result.community_of, options.resolution).unwrap();
            assert!((standalone - result.modularity).abs() < 1e-12);
            assert_eq!(result, louvain_communities(&graph, options).unwrap());
            input.push_str("L\n");
            oracle::graph(&mut input, &graph, |w| *w);
            writeln!(
                input,
                "{} {} {} {}",
                options.resolution,
                options.threshold,
                options.seed,
                options.max_levels.unwrap_or(0)
            )
            .unwrap();
            for id in result.community_of {
                write!(input, "{id} ").unwrap();
            }
            input.push('\n');
            cases.push(true);
            cases_modularity.push((cases.len() - 1, standalone));
        }
    }

    let named: [(&str, Graph<()>); 14] = [
        ("bull_graph", bull_graph()),
        ("chvatal_graph", chvatal_graph()),
        ("diamond_graph", diamond_graph()),
        ("frucht_graph", frucht_graph()),
        ("heawood_graph", heawood_graph()),
        ("hoffman_singleton_graph", hoffman_singleton_graph()),
        ("house_x_graph", house_x_graph()),
        ("icosahedral_graph", icosahedral_graph()),
        ("krackhardt_kite_graph", krackhardt_kite_graph()),
        ("pappus_graph", pappus_graph()),
        ("sedgewick_maze_graph", sedgewick_maze_graph()),
        ("truncated_cube_graph", truncated_cube_graph()),
        ("truncated_tetrahedron_graph", truncated_tetrahedron_graph()),
        ("tutte_graph", tutte_graph()),
    ];
    for (name, graph) in named {
        input.push_str("G\n");
        oracle::graph(&mut input, &graph, |_| 1.0);
        writeln!(input, "{name} 0").unwrap();
        cases.push(true);
    }
    for rows in 1..7 {
        for cols in 1..8 {
            for periodic in [false, true] {
                for (name, graph) in [
                    (
                        "triangular_lattice_graph",
                        triangular_lattice_graph(rows, cols, periodic),
                    ),
                    (
                        "hexagonal_lattice_graph",
                        hexagonal_lattice_graph(rows, cols, periodic),
                    ),
                ] {
                    let Ok(graph) = graph else { continue };
                    input.push_str("G\n");
                    oracle::graph(&mut input, &graph, |_| 1.0);
                    writeln!(input, "{name} 3 {rows} {cols} {}", u8::from(periodic)).unwrap();
                    cases.push(true);
                }
            }
        }
    }
    // Exact agreement on graphs with clear separated communities.
    for kind in [GraphKind::Undirected, GraphKind::Directed] {
        let base = caveman_graph(12, 6).unwrap();
        let graph = if kind == GraphKind::Directed {
            base.to_directed()
        } else {
            base
        };
        let options = LouvainOptions::default();
        let result = louvain_communities_by(&graph, |_| 1.0, options).unwrap();
        assert_eq!(result.community_count(), 12);
        input.push_str("L\n");
        oracle::graph(&mut input, &graph, |_| 1.0);
        writeln!(input, "1 1e-7 0 0").unwrap();
        for id in result.community_of {
            write!(input, "{id} ").unwrap();
        }
        input.push('\n');
        cases.push(true);
        cases_modularity.push((cases.len() - 1, result.modularity));
        quality.push((cases.len() - 1, result.modularity));
    }
    let query = format!("{}\n{input}", cases.len());
    let (version, results) = oracle::invoke("1", &query);
    assert_eq!(results.len(), cases.len());
    for &(index, expected) in &cases_modularity {
        let independently_evaluated: f64 = results[index][3].parse().unwrap();
        assert!(
            (expected - independently_evaluated).abs() < 1e-12,
            "modularity case {index}"
        );
    }
    for (id, expected) in quality {
        let score: f64 = results[id][1].parse().unwrap();
        assert!((expected - score).abs() < 1e-12);
    }
    let indices: BTreeSet<_> = cases_modularity.iter().map(|&(i, _)| i).collect();
    for (index, (expected, answer)) in cases.into_iter().zip(results).enumerate() {
        if !indices.contains(&index) {
            assert_eq!(expected, answer[1] == "1", "case {index}");
        }
    }
    println!(
        "NetworkX {version}, 210 VF2++ comparisons, {} Louvain modularity checks, {} total cases",
        indices.len(),
        query.lines().next().unwrap()
    );
}

#[test]
#[ignore = "requires Python and NetworkX"]
fn left_right_agrees_on_larger_projected_multigraphs() {
    let mut samples = Vec::new();
    for trial in 0..180 {
        let n = 8 + trial % 40;
        let kind = if trial % 2 == 0 {
            GraphKind::Directed
        } else {
            GraphKind::Undirected
        };
        let base = if trial % 3 == 0 {
            triangular_lattice_graph(n / 4, 5 + n % 6, false).unwrap()
        } else {
            gnp_random_graph(n, 0.025 * (1 + trial % 8) as f64, kind, 4391 + trial as u64).unwrap()
        };
        let n = base.node_count();
        let permute = |u| n - 1 - (u + trial) % n;
        let mut builder = Graph::builder(n, kind);
        for (i, edge) in base.edges().iter().enumerate().rev() {
            if i % 41 == 0 {
                continue;
            }
            let (u, v) = (permute(edge.from), permute(edge.to));
            builder.add_edge(u, v, ()).unwrap();
            if i % 7 == 0 {
                builder.add_edge(v, u, ()).unwrap();
            }
            if i % 19 == 0 {
                builder.add_edge(u, u, ()).unwrap();
            }
        }
        samples.push(builder.build());
    }
    let mut query = format!("{}\n", samples.len());
    for graph in &samples {
        query.push_str("P\n");
        oracle::graph(&mut query, graph, |_| 1.0);
    }
    let (version, results) = oracle::invoke("1", &query);
    assert_eq!(samples.len(), results.len());
    for (graph, expected) in samples.iter().zip(results) {
        verify_rotation(graph, expected[1] == "1");
    }
    println!(
        "NetworkX {version}, {} larger directed/undirected multigraph planarity checks",
        samples.len()
    );
}

#[test]
fn degree_sequences_and_seeded_small_world_models() {
    let degrees = [1, 4, 1, 2, 1, 1];
    let tree = degree_sequence_tree(&degrees).unwrap();
    assert_eq!(
        tree.nodes()
            .map(|u| tree.degree(u).unwrap())
            .collect::<Vec<_>>(),
        degrees
    );
    assert_eq!(weakly_connected_components(&tree).count, 1);
    assert!(degree_sequence_tree(&[1, 1, 1]).is_none());
    let ring = watts_strogatz_graph(55, 6, 0.0, 91).unwrap();
    let shortcuts = newman_watts_strogatz_graph(55, 6, 1.0, 91).unwrap();
    assert_eq!(shortcuts.edge_count(), ring.edge_count() * 2);
    assert_eq!(&shortcuts.edges()[..ring.edge_count()], ring.edges());
    assert!(
        connected_watts_strogatz_graph(55, 6, 0.4, 8, 91)
            .unwrap()
            .is_some()
    );
    assert!(
        connected_watts_strogatz_graph(55, 1, 0.4, 8, 91)
            .unwrap()
            .is_none()
    );
    for (iterations, n, m) in [(0, 2, 1), (1, 3, 3), (2, 6, 9), (3, 15, 27)] {
        let graph = dorogovtsev_goltsev_mendes_graph(iterations).unwrap();
        assert_eq!((graph.node_count(), graph.edge_count()), (n, m));
    }
}
