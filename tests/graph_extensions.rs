use std::collections::HashSet;

use graph_library::clustering::{clustering_coefficients, k_truss, transitivity, triangles};
use graph_library::connectivity::{
    DisjointSet, biconnected_components, condensation, cycle_basis, edge_connectivity,
    eulerian_circuit, immediate_dominators, lexicographical_topological_sort, minimum_st_edge_cut,
    reverse_topological_sort, st_connected, topological_generations, transitive_closure_dag,
    transitive_reduction,
};
use graph_library::functions::{AttributeMap, edge_subgraph, induced_subgraph, restricted_view};
use graph_library::generators::*;
use graph_library::shortest_path::DagShortestPathIndex;
use graph_library::{Graph, GraphError, GraphKind};

fn build<W>(
    n: usize,
    k: GraphKind,
    edges: impl IntoIterator<Item = (usize, usize, W)>,
) -> Graph<W> {
    let mut b = Graph::builder(n, k);
    for (u, v, w) in edges {
        b.add_edge(u, v, w).unwrap();
    }
    b.build()
}

#[test]
fn graph_functions_preserve_ids_and_parallel_payloads() {
    let mut builder = Graph::builder(5, GraphKind::Undirected);
    builder.add_path(&[0, 1, 2], 5_i64).unwrap();
    builder.add_cycle(&[2, 3, 4], 6).unwrap();
    builder.add_star(&[1, 2, 3], -3).unwrap();
    assert!(builder.add_path(&[0, 2, 7], 7).is_err());
    let loop_id = builder.add_edge(3, 3, 2).unwrap();
    let graph = builder.build();
    assert_eq!(graph.edge_count(), 8);
    assert_eq!(graph.selfloop_edges().count(), 1);
    assert_eq!(graph.nodes_with_selfloops().collect::<Vec<_>>(), [3]);
    assert_eq!(graph.degree(3), Ok(5));
    assert_eq!(graph.weighted_degree(3), Ok(13));
    assert!(graph.is_frozen() && graph.is_negatively_weighted());
    assert_eq!(graph.path_weight(&[0, 1, 2, 3]), Some(8));
    assert!(graph.is_path(&[2, 1, 3, 3]));
    assert_eq!(graph.common_neighbors(2, 3).unwrap(), [4, 1]);
    assert_eq!(
        graph.non_neighbors(0).unwrap().collect::<Vec<_>>(),
        [2, 3, 4]
    );
    assert!(
        graph
            .non_edges()
            .all(|(u, v)| u < v && !graph.is_path(&[u, v]))
    );
    assert_eq!(graph.to_directed().edge_count(), 15);
    assert_eq!(graph.to_directed().to_undirected().edge_count(), 15);

    let copy = induced_subgraph(&graph, [1, 2, 3]).unwrap().materialize();
    assert_eq!(copy.original_nodes, [1, 2, 3]);
    assert_eq!(copy.original_edges, [1, 2, 5, 6, 7]);
    let view = restricted_view(&graph, [3], [1]).unwrap();
    assert_eq!(
        view.edges().map(|(id, _)| id).collect::<Vec<_>>(),
        [0, 4, 5]
    );
    let selected = edge_subgraph(&graph, [loop_id, 4, loop_id])
        .unwrap()
        .materialize();
    assert_eq!(selected.original_nodes, [2, 3, 4]);
    assert_eq!(selected.original_edges, [4, loop_id]);
    assert!(edge_subgraph(&graph, [99]).is_err());
    let mut data = AttributeMap::for_edges(&graph);
    assert_eq!(data.set(loop_id, "loop"), Ok(None));
    assert_eq!(data.get(loop_id), Some(&"loop"));
    assert!(matches!(
        data.set(99, "bad"),
        Err(GraphError::InvalidEdge { .. })
    ));
}

#[test]
fn dag_queries_and_reduction_preserve_reachability() {
    let dag = build(
        7,
        GraphKind::Directed,
        [
            (5, 0, 4),
            (5, 0, 5),
            (5, 1, 9),
            (5, 3, 20),
            (0, 1, -10),
            (1, 3, 1),
            (2, 3, -4),
            (3, 4, 2),
        ],
    );
    let index = DagShortestPathIndex::new(&dag).unwrap();
    assert_eq!(index.between(5, 4).unwrap().unwrap().cost, -3);
    assert!(index.between(4, 5).unwrap().is_none());
    assert_eq!(
        index
            .from_sources([5, 2])
            .unwrap()
            .path_to(3)
            .unwrap()
            .nodes,
        [5, 0, 1, 3]
    );
    assert_eq!(index.minimum_path().unwrap().cost, -10);
    let rev = reverse_topological_sort(&dag).unwrap();
    let mut pos = [0; 7];
    for (index, node) in rev.into_iter().enumerate() {
        pos[node] = index;
    }
    assert!(dag.edges().iter().all(|edge| pos[edge.from] > pos[edge.to]));
    assert_eq!(
        lexicographical_topological_sort(&dag, |u| u).unwrap(),
        [2, 5, 0, 1, 3, 4, 6]
    );
    assert_eq!(
        topological_generations(&dag).unwrap(),
        [vec![2, 5, 6], vec![0], vec![1], vec![3], vec![4]]
    );
    let retained = transitive_reduction(&dag).unwrap();
    assert_eq!(retained, [0, 4, 5, 6, 7]);
    let reduced = build(
        7,
        GraphKind::Directed,
        retained.into_iter().map(|id| {
            let edge = dag.edge(id).unwrap();
            (edge.from, edge.to, edge.weight)
        }),
    );
    let closure = transitive_closure_dag(&dag).unwrap();
    assert_eq!(
        closure.edges(),
        transitive_closure_dag(&reduced).unwrap().edges()
    );
    for u in dag.nodes() {
        for v in dag.nodes() {
            assert_eq!(
                st_connected(&dag, u, v).unwrap(),
                u == v || closure.is_path(&[u, v])
            );
        }
    }
}

#[test]
fn connectivity_cuts_cycles_and_clustering_respect_multigraphs() {
    let graph = build(
        6,
        GraphKind::Undirected,
        [
            (0, 1, ()),
            (0, 1, ()),
            (1, 2, ()),
            (2, 0, ()),
            (2, 3, ()),
            (3, 4, ()),
            (4, 5, ()),
            (5, 3, ()),
            (5, 5, ()),
        ],
    );
    assert_eq!(edge_connectivity(&graph).unwrap(), 1);
    assert_eq!(minimum_st_edge_cut(&graph, 0, 4).unwrap().edges, [4]);
    assert_eq!(minimum_st_edge_cut(&graph, 0, 1).unwrap().edges.len(), 3);
    let blocks = biconnected_components(&graph).unwrap();
    assert_eq!(blocks.articulation_points, [2, 3]);
    let mut blocks = blocks.edge_blocks;
    blocks.iter_mut().for_each(|block| block.sort_unstable());
    blocks.sort_unstable();
    assert_eq!(blocks, [vec![0, 1, 2, 3], vec![4], vec![5, 6, 7], vec![8]]);
    let cycles = cycle_basis(&graph).unwrap();
    assert_eq!(cycles.len(), 4);
    for cycle in cycles {
        assert_eq!(cycle.nodes.first(), cycle.nodes.last());
        for (&id, pair) in cycle.edges.iter().zip(cycle.nodes.windows(2)) {
            let e = graph.edge(id).unwrap();
            assert!([e.from, e.to] == pair || [e.to, e.from] == pair);
        }
    }
    assert_eq!(triangles(&graph).unwrap(), [1; 6]);
    assert_eq!(k_truss(&graph, 3).unwrap(), [0, 2, 3, 5, 6, 7]);
    assert!(k_truss(&graph, 4).unwrap().is_empty());
    assert_eq!(clustering_coefficients(&graph).unwrap()[2], 1.0 / 3.0);
    assert_eq!(transitivity(&graph).unwrap(), 0.6);

    let graph = build(3, GraphKind::Directed, [(0, 1, ()), (1, 0, ()), (0, 0, ())]);
    let circuit = eulerian_circuit(&graph, Some(0)).unwrap().unwrap();
    assert_eq!(circuit.edges.len(), 3);
    assert_eq!(circuit.nodes.first(), circuit.nodes.last());
    assert_eq!(condensation(&graph).components.count, 2);
    let mut sets = DisjointSet::new(6);
    assert!(sets.union(1, 2).is_some());
    sets.union(4, 5);
    sets.union(2, 4);
    assert_eq!(sets.union(1, 5), None);
    assert_eq!(sets.component_size(5), 4);
    assert_eq!(sets.component_count(), 3);
}

#[test]
fn dominators_agree_with_vertex_deletion() {
    for seed in 0..100 {
        let graph = gnp_random_graph(12, 0.21, GraphKind::Directed, seed).unwrap();
        let source = 0;
        let result = immediate_dominators(&graph, source).unwrap();
        for target in 1..graph.node_count() {
            let reachable = st_connected(&graph, source, target).unwrap();
            assert_eq!(result[target].is_some(), reachable);
            if !reachable {
                continue;
            }
            let mut chain = HashSet::new();
            let mut u = result[target].unwrap();
            loop {
                assert!(chain.insert(u));
                if u == source {
                    break;
                }
                u = result[u].unwrap();
            }
            for candidate in 1..graph.node_count() {
                if candidate == target {
                    continue;
                }
                let pruned = build(
                    graph.node_count(),
                    GraphKind::Directed,
                    graph
                        .edges()
                        .iter()
                        .filter(|e| e.from != candidate && e.to != candidate)
                        .map(|e| (e.from, e.to, ())),
                );
                assert_eq!(
                    chain.contains(&candidate),
                    !st_connected(&pruned, source, target).unwrap()
                );
            }
        }
    }
}

#[test]
fn seeded_and_classic_generators_match_model_invariants() {
    assert_eq!(balanced_tree(3, 3).unwrap().edge_count(), 39);
    assert_eq!(binomial_tree(7).unwrap().degree(0).unwrap(), 7);
    assert_eq!(grid_2d_graph(3, 4, [true, false]).unwrap().edge_count(), 21);
    assert_eq!(complete_bipartite_graph(5, 3).unwrap().edge_count(), 15);
    assert_eq!(circular_ladder_graph(5).unwrap().edge_count(), 15);
    assert_eq!(petersen_graph().degree_histogram(), [0, 0, 0, 10]);
    assert_eq!(barbell_graph(4, 3).unwrap().edge_count(), 16);
    assert_eq!(
        line_graph(&cycle_graph(6, GraphKind::Undirected)).edge_count(),
        6
    );

    for seed in 0..30 {
        let graph = barabasi_albert_graph(50, 3, seed).unwrap();
        assert_eq!(graph.edge_count(), 3 + 46 * 3);
        assert_eq!(
            graph.edges(),
            barabasi_albert_graph(50, 3, seed).unwrap().edges()
        );
        let tree = random_labeled_tree(50, seed);
        assert_eq!(tree.edge_count(), 49);
        assert!(st_connected(&tree, 0, 49).unwrap());
        assert!(cycle_basis(&tree).unwrap().is_empty());
        let world = watts_strogatz_graph(60, 14, 0.8, seed).unwrap();
        assert_eq!(world.edge_count(), 60 * 7);
        assert_eq!(
            world
                .edges()
                .iter()
                .map(|e| (e.from.min(e.to), e.from.max(e.to)))
                .collect::<HashSet<_>>()
                .len(),
            world.edge_count()
        );
    }

    for m in [0, 11, 150, 190] {
        let graph = gnm_random_graph(20, m, GraphKind::Undirected, 6503).unwrap();
        assert_eq!(graph.edge_count(), m);
        assert_eq!(
            graph
                .edges()
                .iter()
                .map(|e| (e.from, e.to))
                .collect::<HashSet<_>>()
                .len(),
            m
        );
    }
    let degrees = [3, 2, 3, 2, 0];
    for graph in [
        configuration_model(&degrees, 94).unwrap(),
        havel_hakimi_graph(&degrees).unwrap(),
    ] {
        assert_eq!(
            graph
                .nodes()
                .map(|u| graph.degree(u).unwrap())
                .collect::<Vec<_>>(),
            degrees
        );
    }
    let partition = planted_partition_graph(4, 5, 1.0, 0.0, GraphKind::Undirected, 29).unwrap();
    assert_eq!(partition.edge_count(), 40);
    assert!(partition.edges().iter().all(|e| e.from / 5 == e.to / 5));

    let geometric = random_geometric_graph(100, 0.1, 83112).unwrap();
    let projected: HashSet<_> = geometric
        .graph
        .edges()
        .iter()
        .map(|e| (e.from, e.to))
        .collect();
    let mut expected = HashSet::new();
    for (i, a) in geometric.positions.iter().enumerate() {
        for (j, b) in geometric.positions.iter().enumerate().skip(i + 1) {
            if (a[0] - b[0]).hypot(a[1] - b[1]) <= 0.1 {
                expected.insert((i, j));
            }
        }
    }
    assert_eq!(projected, expected);
    assert_eq!(
        geometric_edges(&[[0.0, -0.0], [-0.0, 0.0]], 0.0).unwrap(),
        [(0, 1)]
    );
}

#[test]
fn dynamic_connectivity_handles_duplicates_and_rollback() {
    use graph_library::connectivity::{
        ConnectivityAnswer as A, ConnectivityEvent as E, offline_dynamic_connectivity,
    };
    let events = [
        E::ComponentCount,
        E::Add(0, 1),
        E::Add(1, 0),
        E::Add(1, 2),
        E::Add(3, 3),
        E::Connected(0, 2),
        E::Remove(0, 1),
        E::Connected(0, 2),
        E::Remove(1, 0),
        E::Connected(0, 2),
        E::ComponentCount,
        E::Add(2, 3),
        E::Connected(1, 3),
    ];
    assert_eq!(
        offline_dynamic_connectivity(4, &events).unwrap(),
        [
            A::ComponentCount(4),
            A::Connected(true),
            A::Connected(true),
            A::Connected(false),
            A::ComponentCount(3),
            A::Connected(true),
        ]
    );
    assert!(offline_dynamic_connectivity(4, &[E::Remove(0, 2)]).is_err());
}

#[test]
fn distance_centralities_count_shortest_paths_and_orient_directed_edges() {
    use graph_library::centrality::{
        betweenness_centrality, closeness_centrality, harmonic_centrality,
    };
    let path = path_graph(4, GraphKind::Undirected);
    assert_eq!(betweenness_centrality(&path, false), [0.0, 2.0, 2.0, 0.0]);
    assert_eq!(closeness_centrality(&path), [0.5, 0.75, 0.75, 0.5]);
    assert_eq!(
        harmonic_centrality(&path),
        [1.0 + 0.5 + 1.0 / 3.0, 2.5, 2.5, 1.0 + 0.5 + 1.0 / 3.0]
    );
    let directed = build(
        4,
        GraphKind::Directed,
        [(0, 1, ()), (0, 1, ()), (0, 2, ()), (1, 3, ()), (2, 3, ())],
    );
    let scores = betweenness_centrality(&directed, false);
    assert_eq!(scores, [0.0, 2.0 / 3.0, 1.0 / 3.0, 0.0]);
    assert_eq!(harmonic_centrality(&directed), [0.0, 1.0, 1.0, 2.5]);
}

#[test]
fn k_shortest_paths_agree_with_full_simple_path_enumeration() {
    use graph_library::shortest_path::k_shortest_simple_paths;
    for seed in 0..70 {
        let kind = if seed % 2 == 0 {
            GraphKind::Undirected
        } else {
            GraphKind::Directed
        };
        let graph = gnp_random_graph(7, 0.42, kind, seed).unwrap();
        let graph = build(
            7,
            kind,
            graph
                .edges()
                .iter()
                .enumerate()
                .map(|(id, e)| (e.from, e.to, (id % 5) as u64))
                .chain([(0, 1, 1), (0, 1, 0), (2, 2, 0)]),
        );
        let mut paths = vec![(0, 0, 1_u16, Vec::new())];
        let mut exhaustive = Vec::new();
        while let Some((node, cost, visited, edges)) = paths.pop() {
            if node == 6 {
                exhaustive.push((cost, edges));
                continue;
            }
            for arc in graph.outgoing(node).unwrap() {
                if visited & (1 << arc.node) != 0 {
                    continue;
                }
                let mut next = edges.clone();
                next.push(arc.edge);
                paths.push((
                    arc.node,
                    cost + u128::from(graph.edges()[arc.edge].weight),
                    visited | (1 << arc.node),
                    next,
                ));
            }
        }
        exhaustive.sort_unstable();
        let result = k_shortest_simple_paths(&graph, 0, 6, usize::MAX).unwrap();
        let mut result: Vec<_> = result
            .into_iter()
            .map(|path| (path.cost, path.edges))
            .collect();
        assert!(result.is_sorted_by_key(|path| path.0));
        result.sort_unstable();
        assert_eq!(result, exhaustive);
    }
}

#[test]
fn disjoint_paths_and_vertex_cuts_agree_with_deletion() {
    use graph_library::connectivity::{
        edge_disjoint_paths, minimum_st_vertex_cut, node_disjoint_paths,
    };
    for seed in 0..160 {
        let kind = if seed % 2 == 0 {
            GraphKind::Undirected
        } else {
            GraphKind::Directed
        };
        let sample = gnp_random_graph(8, 0.39, kind, seed).unwrap();
        let graph = build(
            8,
            kind,
            sample.edges().iter().map(|e| (e.from, e.to, ())).chain([
                (0, 1, ()),
                (0, 1, ()),
                (7, 7, ()),
            ]),
        );
        let source = 0;
        let sink = 7;
        let family = edge_disjoint_paths(&graph, source, sink).unwrap();
        assert_eq!(
            family.len(),
            minimum_st_edge_cut(&graph, source, sink)
                .unwrap()
                .edges
                .len()
        );
        let verify = |path: &graph_library::connectivity::EdgeWalk| {
            assert_eq!(path.nodes.first(), Some(&source));
            assert_eq!(path.nodes.last(), Some(&sink));
            assert_eq!(
                path.nodes.len(),
                path.nodes.iter().collect::<HashSet<_>>().len()
            );
            for (&id, pair) in path.edges.iter().zip(path.nodes.windows(2)) {
                let e = graph.edges()[id];
                assert!(
                    [e.from, e.to] == pair
                        || kind == GraphKind::Undirected && [e.to, e.from] == pair
                );
            }
        };
        let mut logical = HashSet::new();
        for path in &family {
            verify(path);
            assert!(path.edges.iter().all(|&id| logical.insert(id)));
        }
        let family = node_disjoint_paths(&graph, source, sink).unwrap();
        let mut interior = HashSet::new();
        for path in &family {
            verify(path);
            assert!(
                path.nodes[1..path.nodes.len() - 1]
                    .iter()
                    .all(|&node| interior.insert(node))
            );
        }
        let cut = minimum_st_vertex_cut(&graph, source, sink).unwrap();
        if graph.is_path(&[source, sink]) {
            assert!(cut.is_none());
        } else {
            let cut = cut.unwrap();
            assert_eq!(cut.len(), family.len());
            let pruned = build(
                8,
                kind,
                graph
                    .edges()
                    .iter()
                    .filter(|e| !cut.contains(&e.from) && !cut.contains(&e.to))
                    .map(|e| (e.from, e.to, ())),
            );
            assert!(!st_connected(&pruned, source, sink).unwrap());
        }
    }
}
