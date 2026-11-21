use graph_library::connectivity::*;
use graph_library::flow::{max_flow, max_flow_capacity_scaling, max_flow_ford_fulkerson};
use graph_library::generators::*;
use graph_library::shortest_path::{dijkstra, dijkstra_fibonacci};
use graph_library::structures::*;
use graph_library::traversal::{direction_optimizing_bfs, multi_source_bfs};
use graph_library::{Graph, GraphError, GraphKind};
use std::collections::BTreeSet;

#[test]
fn indices_and_distances_match_searches() {
    for kind in [GraphKind::Directed, GraphKind::Undirected] {
        for seed in 0..48 {
            let g = gnp_random_graph(35, (seed % 10) as f64 / 80., kind, seed).unwrap();
            let mut b = Graph::builder(35, kind);
            for (id, e) in g.edges().iter().enumerate() {
                b.add_edge(e.from, e.to, (id % 9) as u64).unwrap();
            }
            let g = b.build();
            let chains =
                ChainReachabilityIndex::with_budget(&g, 35 * 35 * size_of::<usize>()).unwrap();
            let weak = weakly_connected_components(&g);
            let scc = strongly_connected_components(&g);
            for dims in [0, 1, 4] {
                let index = ReachabilityIndex::new(&g, dims, seed);
                assert!(index.reachable(35, 0).is_err());
                for u in g.nodes() {
                    let reference = multi_source_bfs(&g, &[u]).unwrap();
                    let bfs = direction_optimizing_bfs(&g, u).unwrap();
                    let f = dijkstra_fibonacci(&g, u).unwrap();
                    let d = dijkstra(&g, u).unwrap();
                    for v in g.nodes() {
                        let reaches = reference.distance(v).is_some();
                        assert_eq!(index.reachable(u, v).unwrap(), reaches);
                        assert_eq!(chains.reachable(u, v).unwrap(), reaches);
                        assert_eq!(reference.distance(v), bfs.distance(v));
                        assert_eq!(d.distance(v), f.distance(v));
                        assert_eq!(
                            scc.strongly_connected(u, v).unwrap(),
                            reaches && index.reachable(v, u).unwrap()
                        );
                        if !g.is_directed() {
                            assert_eq!(weak.connected(u, v).unwrap(), reaches);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn sparse_index_handles_a_long_dag() {
    let n = 100_000;
    let g = path_graph(n, GraphKind::Directed);
    let bytes = n * size_of::<usize>();
    assert!(ChainReachabilityIndex::with_budget(&g, bytes - 1).is_none());
    let chains = ChainReachabilityIndex::with_budget(&g, bytes).unwrap();
    assert_eq!(chains.chain_count(), 1);
    let index = ReachabilityIndex::new(&g, 3, 71);
    for k in (1..n).step_by(317) {
        assert!(chains.reachable(0, k).unwrap());
        assert!(!chains.reachable(k, 0).unwrap());
        assert!(index.reachable(0, k).unwrap());
        assert!(!index.reachable(k, 0).unwrap());
    }
}

#[test]
fn augmentations_agree_on_integral_multigraphs() {
    for seed in 0..150 {
        let sample = gnp_random_graph(12, 0.27, GraphKind::Directed, seed).unwrap();
        let mut b = Graph::builder(12, GraphKind::Directed);
        for (i, e) in sample.edges().iter().enumerate() {
            b.add_edge(e.from, e.to, (i % 15) as u64).unwrap();
            b.add_edge(e.from, e.to, 1).unwrap();
        }
        let g = b.build();
        let result = max_flow_ford_fulkerson(&g, 0, 11).unwrap();
        assert_eq!(result.value, max_flow(&g, 0, 11).unwrap().value);
        let mut balance = [0i128; 12];
        let mut cut = 0u128;
        for (e, f) in g.edges().iter().zip(result.edge_flows) {
            assert!(f <= e.weight);
            balance[e.from] += i128::from(f);
            balance[e.to] -= i128::from(f);
            if result.source_side[e.from] && !result.source_side[e.to] {
                cut += u128::from(e.weight);
            }
        }
        assert_eq!(cut, result.value);
        assert!(balance[1..11].iter().all(|&x| x == 0));
        assert_eq!(balance[0], result.value as i128);
        assert_eq!(balance[11], -balance[0]);
    }
    let g = path_graph(2, GraphKind::Undirected);
    let mut b = Graph::builder(g.node_count(), g.kind());
    b.add_edge(0, 1, 1).unwrap();
    assert_eq!(
        max_flow_capacity_scaling(&b.build(), 0, 1),
        Err(GraphError::RequiresDirected)
    );
}

#[test]
fn indexed_heaps_match_an_ordered_reference() {
    let n = 97;
    let mut state = 51_u64;
    let mut random = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut fib = FibonacciHeap::new(n);
    let mut radix = RadixQueue::new(n);
    let mut reference: BTreeSet<(u128, usize)> = BTreeSet::new();
    let mut keys = vec![None; n];
    let mut floor = 0;
    for _ in 0..35_000 {
        if random() % 5 == 0 && !reference.is_empty() {
            let min = reference.first().unwrap().0;
            let (key, id) = fib.pop().unwrap();
            assert_eq!(key, min);
            assert_eq!(radix.pop(), Some((min, id))); // generated keys encode the ID
            assert!(reference.remove(&(key, id)));
            keys[id] = None;
            floor = min;
        } else {
            let id = random() as usize % n;
            let key = ((floor >> 8) + u128::from(random() % 1000) + 1) * 256 + id as u128;
            let accepted = keys[id].is_none_or(|old| key < old);
            assert_eq!(fib.decrease(id, key), accepted);
            assert_eq!(radix.decrease(id, key), accepted);
            if accepted {
                if let Some(old) = keys[id].replace(key) {
                    reference.remove(&(old, id));
                }
                reference.insert((key, id));
            }
        }
        assert_eq!(fib.len(), reference.len());
        assert_eq!(radix.len(), reference.len());
        assert_eq!(fib.peek(), reference.first().copied());
    }
    assert!(!radix.decrease(0, floor - 1));
    radix.clear();
    assert!(radix.decrease(0, u128::MAX));
    assert!(radix.decrease(0, u128::MAX - 1));
    assert!(radix.decrease(1, u128::MAX));
    assert_eq!(radix.pop(), Some((u128::MAX - 1, 0)));
    assert_eq!(radix.pop(), Some((u128::MAX, 1)));
    assert_eq!(radix.pop(), None);
}

#[test]
fn ordered_ranges_and_tree_boundaries() {
    let combine = |a: &String, b: &String| format!("{a}{b}");
    for n in 0..18 {
        let mut words: Vec<_> = (0..n).map(|i| format!("[{i}]")).collect();
        let mut tree = SegmentTree::new(&words, String::new(), combine);
        for i in 0..=n {
            for j in i..=n {
                assert_eq!(tree.query(i..j).unwrap(), words[i..j].concat());
            }
            if i < n {
                assert!(tree.set(i, "x".into()));
                words[i] = "x".into();
            }
        }
        assert!(!tree.set(n, "".into()));
        assert!(tree.query(0..n + 1).is_none());
    }
    let mut t = KaryTree::new(2, 'a');
    assert_eq!(t.add_child(0, 'b'), Ok(1));
    assert_eq!(t.add_child(0, 'c'), Ok(2));
    assert_eq!(t.add_child(0, 'd'), Err('d'));
    assert_eq!(t.add_child(30, 'z'), Err('z'));
    assert_eq!(t.children(0), Some(&[1, 2][..]));
    assert_eq!(t.parent(2), Some(0));
    assert_eq!(t.value(1), Some(&'b'));
}

#[test]
fn generated_families_and_sampled_forests() {
    for seed in 0..35 {
        for d in [0, 1, 3, 12] {
            let g = random_regular_bipartite_multigraph(10, d, seed).unwrap();
            assert!(g.nodes().all(|v| g.degree(v).unwrap() == d));
            assert_eq!(
                g.edges(),
                random_regular_bipartite_multigraph(10, d, seed)
                    .unwrap()
                    .edges()
            );
            let b = graph_library::matching::BipartiteGraph::from_edges(
                10,
                10,
                g.edges().iter().map(|e| (e.from, e.to - 10)),
            )
            .unwrap();
            assert_eq!(
                b.random_walk_maximum_matching(seed).len(),
                if d == 0 { 0 } else { 10 }
            );
            let ids = random_spanning_forest(&g, seed).unwrap();
            let mut sets = UnionFind::new(g.node_count());
            for id in ids {
                let e = g.edges()[id];
                assert!(sets.union(e.from, e.to).is_some());
            }
            assert_eq!(
                sets.component_count(),
                weakly_connected_components(&g).count
            );
        }
    }
    assert_eq!(
        bipartite_gnp_random_graph(8, 7, 1., 0)
            .unwrap()
            .edge_count(),
        56
    );
    assert_eq!(
        bipartite_gnp_random_graph(8, 0, 0.5, 0)
            .unwrap()
            .edge_count(),
        0
    );
    assert!(bipartite_gnp_random_graph(3, 4, f64::NAN, 0).is_err());
    assert_eq!(
        random_intersection_graph(9, 1, 1., 42)
            .unwrap()
            .edge_count(),
        36
    );
    assert_eq!(
        random_intersection_graph(9, 0, 1., 42)
            .unwrap()
            .edge_count(),
        0
    );
    let g = random_intersection_graph(100, 20, 0.11, 15).unwrap();
    assert_eq!(
        g.edges(),
        random_intersection_graph(100, 20, 0.11, 15)
            .unwrap()
            .edges()
    );
    assert!(g.edges().iter().all(|e| e.from < e.to));
    assert_eq!(
        g.edges()
            .iter()
            .map(|e| (e.from, e.to))
            .collect::<BTreeSet<_>>()
            .len(),
        g.edge_count()
    );
}

#[test]
fn structural_and_statistical_analytics() {
    use graph_library::{
        centrality::*, clustering::sampled_average_clustering,
        community::label_propagation_communities, isomorphism::rooted_tree_is_isomorphic,
        planarity::check_planarity,
    };
    let clique = complete_graph(7, GraphKind::Undirected).unwrap();
    let tree = path_graph(17, GraphKind::Undirected);
    assert_eq!(sampled_average_clustering(&clique, 500, 13).unwrap(), 1.);
    assert_eq!(sampled_average_clustering(&tree, 500, 13).unwrap(), 0.);
    assert!(sampled_average_clustering(&clique, 0, 0).is_err());
    assert!(rooted_tree_is_isomorphic(&tree, 0, &tree, 16).unwrap());
    assert!(!rooted_tree_is_isomorphic(&tree, 0, &tree, 1).unwrap());
    assert!(!rooted_tree_is_isomorphic(&clique, 0, &clique, 0).unwrap());
    let communities = caveman_graph(4, 5).unwrap();
    let groups = label_propagation_communities(&communities, 91).unwrap();
    assert_eq!(groups.len(), 20);
    assert_eq!(groups.iter().copied().collect::<BTreeSet<_>>().len(), 4);
    assert!(groups.chunks(5).all(|c| c.iter().all(|l| *l == c[0])));
    let exact = betweenness_centrality(&tree, false);
    assert_eq!(
        sampled_betweenness_centrality(&tree, false, 17, 0).unwrap(),
        exact
    );
    let mut mean = [0.; 17];
    for seed in 0..1000 {
        for (x, y) in mean
            .iter_mut()
            .zip(sampled_betweenness_centrality(&tree, false, 2, seed).unwrap())
        {
            *x += y / 1000.;
        }
    }
    assert!(mean.iter().zip(exact).all(|(a, b)| (a - b).abs() < 3.));
    let mut b = Graph::builder(8, GraphKind::Undirected);
    b.add_cycle(&[0, 1, 2], ()).unwrap();
    b.add_edge(2, 3, ()).unwrap();
    b.add_cycle(&[4, 5, 6], ()).unwrap();
    let dual = check_planarity(&b.build()).unwrap().dual_graph();
    assert_eq!(dual.node_count(), 3);
    assert_eq!(dual.edge_count(), 7);
    let bridges: Vec<_> = dual
        .edges()
        .iter()
        .filter(|e| e.from == e.to)
        .map(|e| e.weight)
        .collect();
    assert_eq!(bridges, [3]);
    assert_eq!(
        dual.edges()
            .iter()
            .map(|e| e.weight)
            .collect::<BTreeSet<_>>(),
        (0..7).collect()
    );
}

#[test]
fn grid_connectivity_is_a_snapshot() {
    use graph_library::grid::{Cell, Grid, Movement};
    let mut g = Grid::new(3, 3, Movement::Eight).unwrap();
    let a = Cell::new(0, 0);
    let b = Cell::new(2, 2);
    let before = g.connectivity_index();
    for row in 0..3 {
        g.set_blocked(Cell::new(row, 1), true).unwrap();
    }
    let after = g.connectivity_index();
    assert!(before.connected(a, b).unwrap());
    assert!(!after.connected(a, b).unwrap());
    assert!(after.connected(a, Cell::new(0, 1)).is_err());
    assert!(after.connected(Cell::new(3, 0), a).is_err());
}

#[test]
fn frequencies_and_rollback_queries() {
    let mut f = FenwickTree::new(5);
    assert_eq!(f.select(0), None);
    assert_eq!(f.replace(1, 2), Some(0));
    assert_eq!(f.replace(4, 3), Some(0));
    assert_eq!(
        (0..5).map(|i| f.select(i).unwrap()).collect::<Vec<_>>(),
        [1, 1, 4, 4, 4]
    );
    assert_eq!(f.replace(1, usize::MAX), None);
    assert_eq!(f.total(), 5);
    assert_eq!(f.replace(4, 0), Some(3));
    assert_eq!(f.prefix_sum(5), Some(2));
    let mut sets = RollbackUnionFind::new(8);
    sets.union(1, 2);
    let checkpoint = sets.checkpoint();
    sets.union(2, 3);
    sets.union(3, 2);
    assert_eq!(sets.find(3), sets.find(1));
    sets.rollback(checkpoint);
    assert_eq!(sets.find(2), sets.find(1));
    assert_ne!(sets.find(3), sets.find(1));
    assert_eq!(sets.component_count(), 7);
}

#[test]
fn simple_regular_sampler_and_complement() {
    assert!(random_regular_graph(5, 3, 0, 5).is_err());
    assert!(random_regular_graph(5, 5, 0, 5).is_err());
    assert!(random_regular_graph(8, 3, 0, 0).unwrap().is_none());
    for d in [0, 1, 2, 3, 7, 8, 9] {
        let g = random_regular_graph(10, d, 851, 500).unwrap().unwrap();
        assert_eq!(
            g.edges(),
            random_regular_graph(10, d, 851, 500)
                .unwrap()
                .unwrap()
                .edges()
        );
        assert!(g.nodes().all(|u| g.degree(u).unwrap() == d));
        let pairs: BTreeSet<_> = g
            .edges()
            .iter()
            .map(|e| {
                assert_ne!(e.from, e.to);
                (e.from.min(e.to), e.from.max(e.to))
            })
            .collect();
        assert_eq!(pairs.len(), 5 * d);
    }
}

#[test]
fn static_minima_and_weight_thresholds() {
    let values: Vec<_> = (0..411).map(|i| (i * 19 + i * i) % 47).collect();
    let rmq = RangeMinimumIndex::new(values.clone());
    for i in 0..values.len() {
        let mut best = i;
        for j in i..values.len() {
            if values[j] < values[best] {
                best = j;
            }
            assert_eq!(rmq.arg_min(i..j + 1), Some(best));
        }
    }
    assert_eq!(rmq.arg_min(0..0), None);
    for seed in 0..12 {
        let sample = gnp_random_graph(90, 0.04, GraphKind::Undirected, seed).unwrap();
        let mut b = Graph::builder(90, GraphKind::Undirected);
        for (i, e) in sample.edges().iter().enumerate() {
            b.add_edge(e.from, e.to, (i % 13) as i64 - 5).unwrap();
        }
        let g = b.build();
        let index = ThresholdConnectivityIndex::new(&g).unwrap();
        for threshold in [-6, -2, 1, 5, 8] {
            let mut sets = UnionFind::new(90);
            for e in g.edges().iter().filter(|e| e.weight <= threshold) {
                sets.union(e.from, e.to);
            }
            for u in g.nodes() {
                for v in g.nodes() {
                    assert_eq!(
                        index.connected_at(u, v, threshold).unwrap(),
                        sets.connected(u, v)
                    );
                }
            }
        }
    }
}

#[test]
fn filtered_forests_match_kruskal_including_ties() {
    for seed in 0..60 {
        let n = 200;
        let mut b = Graph::builder(n, GraphKind::Undirected);
        // Cheap disjoint pairs and quadruples force nontrivial sampled subproblems.
        for v in 0..n - 1 {
            b.add_edge(v, v + 1, ((v + 1).trailing_zeros()) as u64)
                .unwrap();
        }
        let extra = gnp_random_graph(n, 0.075, GraphKind::Undirected, seed).unwrap();
        for (id, e) in extra.edges().iter().enumerate() {
            b.add_edge(e.from, e.to, 5 + (id % 8) as u64).unwrap();
            if id % 9 == 0 {
                b.add_edge(e.from, e.to, 5 + (id % 8) as u64).unwrap();
            }
        }
        let g = b.build();
        assert_eq!(
            minimum_spanning_forest_filtered(&g, seed).unwrap(),
            minimum_spanning_forest(&g).unwrap()
        );
    }
    for n in [0, 1, 10] {
        let g = Graph::<u64>::builder(n, GraphKind::Undirected).build();
        assert_eq!(
            minimum_spanning_forest_filtered(&g, 0).unwrap(),
            minimum_spanning_forest(&g).unwrap()
        );
    }
}
