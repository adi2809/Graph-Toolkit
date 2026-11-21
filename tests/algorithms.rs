use graph_library::centrality::{PageRankError, PageRankOptions, pagerank};
use graph_library::connectivity::{
    bridges_and_articulation_points, core_numbers, minimum_spanning_forest,
    strongly_connected_components, topological_sort, weakly_connected_components,
};
use graph_library::flow::{
    global_min_cut, max_flow, max_flow_boykov_kolmogorov, max_flow_capacity_scaling,
    max_flow_dinic, max_flow_edmonds_karp, max_flow_ford_fulkerson, max_flow_network_simplex,
};
use graph_library::grid::{Cell, Grid, GridError, Movement};
use graph_library::matching::{
    BipartiteGraph, maximum_general_matching, maximum_general_matching_gabow,
};
use graph_library::shortest_path::{
    LandmarkRouter, Path, astar, bidirectional_dijkstra, dijkstra, zero_one_bfs,
};
use graph_library::traversal::{bfs, dfs, multi_source_bfs};
use graph_library::{Graph, GraphError, GraphKind};

fn graph<W>(
    n: usize,
    kind: GraphKind,
    edges: impl IntoIterator<Item = (usize, usize, W)>,
) -> Graph<W> {
    let mut builder = Graph::builder(n, kind);
    for (from, to, weight) in edges {
        builder.add_edge(from, to, weight).unwrap();
    }
    builder.build()
}

struct Samples(u64);

impl Samples {
    fn below(&mut self, upper: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 as usize) % upper
    }
}

#[test]
fn csr_and_traversals_preserve_order_and_edge_identity() {
    let mut builder = Graph::builder(6, GraphKind::Directed);
    assert!(matches!(
        builder.add_edge(1, 6, "invalid"),
        Err(GraphError::InvalidNode { .. })
    ));
    for (from, to) in [(0, 2), (0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (2, 2)] {
        builder.add_edge(from, to, "payload").unwrap();
    }
    let graph = builder.build();
    assert_eq!(graph.edge_count(), 7);
    assert_eq!(graph.edge(0).unwrap().weight, "payload");
    assert_eq!(
        graph
            .outgoing(0)
            .unwrap()
            .iter()
            .map(|arc| arc.edge)
            .collect::<Vec<_>>(),
        [0, 1, 2]
    );
    assert_eq!(
        graph
            .incoming(2)
            .unwrap()
            .iter()
            .map(|arc| arc.node)
            .collect::<Vec<_>>(),
        [0, 0, 2]
    );
    assert_eq!(bfs(&graph, 0).unwrap(), [0, 2, 1, 3, 4]);
    assert_eq!(dfs(&graph, 0).unwrap(), [0, 2, 3, 4, 1]);
    let forest = multi_source_bfs(&graph, &[1, 2, 2]).unwrap();
    assert_eq!(forest.order(), [1, 2, 3, 4]);
    assert_eq!(forest.source(3), Some(1));
    assert_eq!(forest.distance(4), Some(1));
    assert_eq!(forest.path_to(3), Some(vec![1, 3]));
    assert_eq!(forest.path_to(5), None);
    assert!(multi_source_bfs(&graph, &[]).unwrap().order().is_empty());
}

#[test]
fn grid_routes_agree_with_distance_fields() {
    let mut rng = Samples(101);
    for movement in [Movement::Four, Movement::Eight] {
        let mut grid = Grid::new(9, 11, movement).unwrap();
        for row in 0..9 {
            for col in 0..11 {
                grid.set_blocked(Cell::new(row, col), rng.below(4) == 0)
                    .unwrap();
            }
        }
        let source = Cell::new(0, 0);
        let other = Cell::new(8, 10);
        grid.set_blocked(source, false).unwrap();
        grid.set_blocked(other, false).unwrap();
        let a = grid.multi_source_bfs(&[source]).unwrap();
        let b = grid.multi_source_bfs(&[other]).unwrap();
        let multi = grid.multi_source_bfs(&[source, other]).unwrap();
        for row in 0..9 {
            for col in 0..11 {
                let cell = Cell::new(row, col);
                if grid.is_blocked(cell).unwrap() {
                    assert_eq!(
                        grid.shortest_path(source, cell),
                        Err(GridError::Blocked(cell))
                    );
                    assert_eq!(multi.distance(cell), None);
                    continue;
                }
                let route = grid.shortest_path(source, cell).unwrap();
                assert_eq!(route.map(|nodes| nodes.len() - 1), a.distance(cell));
                let expected = a.distance(cell).into_iter().chain(b.distance(cell)).min();
                assert_eq!(multi.distance(cell), expected);
            }
        }
        assert!(!grid.is_blocked(other).unwrap());
    }
    let mut grid = Grid::new(2, 2, Movement::Eight).unwrap();
    let (start, end) = (Cell::new(0, 0), Cell::new(1, 1));
    assert_eq!(grid.shortest_path(start, end).unwrap().unwrap().len(), 2);
    grid.set_blocked(Cell::new(0, 1), true).unwrap();
    assert_eq!(grid.shortest_path(start, end).unwrap().unwrap().len(), 3);
    grid.set_blocked(Cell::new(1, 0), true).unwrap();
    assert_eq!(grid.shortest_path(start, end).unwrap(), None);
}

fn check_path(graph: &Graph<u64>, path: &Path, source: usize, target: usize) {
    assert_eq!(path.nodes.first(), Some(&source));
    assert_eq!(path.nodes.last(), Some(&target));
    assert_eq!(path.edges.len() + 1, path.nodes.len());
    let mut seen = vec![false; graph.node_count()];
    let mut cost = 0;
    for &node in &path.nodes {
        assert!(!seen[node], "cycle in {path:?}");
        seen[node] = true;
    }
    for (nodes, &id) in path.nodes.windows(2).zip(&path.edges) {
        let edge = graph.edge(id).unwrap();
        assert!(
            (edge.from == nodes[0] && edge.to == nodes[1])
                || (graph.kind() == GraphKind::Undirected
                    && edge.to == nodes[0]
                    && edge.from == nodes[1])
        );
        cost += u128::from(edge.weight);
    }
    assert_eq!(path.cost, cost);
}

#[test]
fn shortest_paths_agree_with_floyd_warshall_on_multigraphs() {
    let mut random = Samples(1009);
    for case in 0..160 {
        let n = 2 + random.below(9);
        let kind = if case % 2 == 0 {
            GraphKind::Directed
        } else {
            GraphKind::Undirected
        };
        let weights = [0, 0, 0, 1, 1, 17, 50, u64::MAX];
        let edges = (0..random.below(5 * n)).map(|_| {
            (
                random.below(n),
                random.below(n),
                weights[random.below(weights.len())],
            )
        });
        let graph = graph(n, kind, edges);
        let router = LandmarkRouter::new(&graph, &[0, n / 2, n / 2, n - 1]).unwrap();
        let mut dist = vec![vec![u128::MAX; n]; n];
        for (i, row) in dist.iter_mut().enumerate() {
            row[i] = 0;
            for arc in graph.outgoing(i).unwrap() {
                row[arc.node] = row[arc.node].min(graph.edges()[arc.edge].weight.into());
            }
        }
        for k in 0..n {
            for i in 0..n {
                for j in 0..n {
                    dist[i][j] = dist[i][j].min(dist[i][k].saturating_add(dist[k][j]));
                }
            }
        }
        for (source, row) in dist.iter().enumerate() {
            let paths = dijkstra(&graph, source).unwrap();
            for (target, &cost) in row.iter().enumerate() {
                let expected = (cost != u128::MAX).then_some(cost);
                assert_eq!(paths.distance(target), expected);
                if let Some(path) = paths.path_to(target) {
                    check_path(&graph, &path, source, target);
                }
                let h: Vec<_> = (0..n)
                    .map(|id| {
                        if random.below(2) == 0 {
                            0
                        } else {
                            dist[id][target].min(u128::from(u64::MAX)) as u64
                        }
                    })
                    .collect();
                for path in [
                    astar(&graph, source, target, |node| h[node]).unwrap(),
                    bidirectional_dijkstra(&graph, source, target).unwrap(),
                    router.shortest_path(source, target).unwrap(),
                ] {
                    assert_eq!(
                        path.as_ref().map(|path| path.cost),
                        expected,
                        "case {case} {kind:?} from {source} to {target} on {:?}",
                        graph.edges()
                    );
                    if let Some(path) = path {
                        check_path(&graph, &path, source, target);
                    }
                }
            }
        }
    }
}

#[test]
fn zero_one_search_and_astar_handle_improvements_after_discovery() {
    let binary = graph(
        4,
        GraphKind::Directed,
        [(0, 1, 1), (0, 2, 0), (2, 1, 0), (1, 3, 1)],
    );
    let paths = zero_one_bfs(&binary, 0).unwrap();
    assert_eq!(paths.path_to(3).unwrap().nodes, [0, 2, 1, 3]);
    assert_eq!(paths.distance(3), Some(1));

    let weighted = graph(
        4,
        GraphKind::Directed,
        [(0, 1, 4), (0, 2, 1), (2, 1, 1), (1, 3, 5)],
    );
    let path = astar(&weighted, 0, 3, |node| [0, 0, 6, 0][node])
        .unwrap()
        .unwrap();
    assert_eq!(path.nodes, [0, 2, 1, 3]);
    assert_eq!(path.cost, 7);
    assert_eq!(
        zero_one_bfs(&weighted, 0).unwrap_err(),
        GraphError::NonBinaryWeight { edge: 0, weight: 4 }
    );
    let wide = graph(3, GraphKind::Directed, [(0, 1, u64::MAX), (1, 2, u64::MAX)]);
    assert_eq!(
        bidirectional_dijkstra(&wide, 0, 2).unwrap().unwrap().cost,
        2 * u128::from(u64::MAX)
    );
}

fn brute_matching(adjacency: &[usize], remaining: usize, memo: &mut [Option<usize>]) -> usize {
    if remaining == 0 {
        return 0;
    }
    if let Some(value) = memo[remaining] {
        return value;
    }
    let node = remaining.trailing_zeros() as usize;
    let rest = remaining & !(1 << node);
    let mut best = brute_matching(adjacency, rest, memo);
    let mut choices = adjacency[node] & rest;
    while choices != 0 {
        let other = choices & choices.wrapping_neg();
        best = best.max(1 + brute_matching(adjacency, rest ^ other, memo));
        choices ^= other;
    }
    memo[remaining] = Some(best);
    best
}

fn expected_matching(n: usize, pairs: &[(usize, usize)]) -> usize {
    let mut adjacency = vec![0; n];
    for &(a, b) in pairs {
        adjacency[a] |= 1 << b;
        adjacency[b] |= 1 << a;
    }
    let mut memo = vec![None; 1 << n];
    brute_matching(&adjacency, (1 << n) - 1, &mut memo)
}

#[test]
fn general_matching_agrees_with_exhaustive_search() {
    let n = 6;
    let possible: Vec<_> = (0..n)
        .flat_map(|u| (u + 1..n).map(move |v| (u, v)))
        .collect();
    for mask in 0..1 << possible.len() {
        let pairs: Vec<_> = possible
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(index, pair)| (mask & (1 << index) != 0).then_some(pair))
            .collect();
        check_general_matching(n, &pairs);
    }

    let mut rng = Samples(11351);
    for _ in 0..1_200 {
        let n = 7 + rng.below(6);
        let pairs: Vec<_> = (0..rng.below(n * n))
            .map(|_| (rng.below(n), rng.below(n)))
            .collect();
        check_general_matching(n, &pairs);
    }
}

fn check_general_matching(n: usize, pairs: &[(usize, usize)]) {
    let graph = graph(
        n,
        GraphKind::Undirected,
        pairs.iter().map(|&(a, b)| (a, b, ())),
    );
    let expected = expected_matching(n, pairs);
    for matching in [
        maximum_general_matching(&graph).unwrap(),
        maximum_general_matching_gabow(&graph).unwrap(),
    ] {
        assert_eq!(matching.len(), expected, "{pairs:?}");
        let mut seen = vec![false; n];
        for edge in matching.edges() {
            let (a, b) = pairs[edge];
            assert_ne!(a, b);
            assert!(!seen[a] && !seen[b]);
            seen[a] = true;
            seen[b] = true;
            assert_eq!(matching.mate(a).unwrap().node, b);
            assert_eq!(matching.mate(b).unwrap().edge, edge);
        }
        assert_eq!(matching.edges().count(), matching.len());
        assert_eq!(
            (0..n).filter(|&node| matching.mate(node).is_some()).count(),
            2 * expected
        );
    }
}

#[test]
fn gabow_and_edmonds_agree_on_larger_seeded_multigraphs() {
    let mut rng = Samples(147_220_165);
    for case in 0..300 {
        let n = 20 + rng.below(120);
        let edges = match case % 3 {
            0 => 2 * n,
            1 => n * n / 4,
            _ => n / 4,
        };
        let graph = graph(
            n,
            GraphKind::Undirected,
            (0..edges).map(|_| (rng.below(n), rng.below(n), ())),
        );
        let expected = maximum_general_matching(&graph).unwrap().len();
        let result = maximum_general_matching_gabow(&graph).unwrap();
        assert_eq!(result.len(), expected, "case {case}");
        let mut count = 0;
        for node in 0..n {
            if let Some(mate) = result.mate(node) {
                assert_ne!(node, mate.node);
                let edge = graph.edge(mate.edge).unwrap();
                assert!(
                    (node, mate.node) == (edge.from, edge.to)
                        || (mate.node, node) == (edge.from, edge.to)
                );
                let reverse = result.mate(mate.node).unwrap();
                assert_eq!((reverse.node, reverse.edge), (node, mate.edge));
                count += 1;
            }
        }
        assert_eq!(count, 2 * expected);
    }
}

#[test]
fn gabow_lifts_a_long_path_with_sparse_vertex_ids() {
    let active = 30_000;
    // Visit inner path vertices before the free endpoints in the greedy scan.
    let edges = (0..active - 3)
        .map(|node| (node, node + 1, ()))
        .chain([(active - 2, 0, ()), (active - 3, active - 1, ())]);
    let graph = graph(4 * active, GraphKind::Undirected, edges);
    let matching = maximum_general_matching_gabow(&graph).unwrap();
    assert_eq!(matching.len(), active / 2);
    assert_eq!(matching.mate(active - 2).unwrap().edge, active - 3);
    assert_eq!(matching.mate(active - 1).unwrap().edge, active - 2);
    assert!(matching.mate(4 * active - 1).is_none());
}

#[test]
fn bipartite_strategies_return_exact_matchings_and_covers() {
    let (left, right) = (3, 4);
    for mask in 0..1 << (left * right) {
        let pairs: Vec<_> = (0..left * right)
            .filter(|&bit| (mask & (1 << bit)) != 0)
            .map(|bit| (bit / right, bit % right))
            .collect();
        check_bipartite_matching(left, right, &pairs, mask);
    }

    let mut rng = Samples(201531);
    for seed in 0..300 {
        let (left, right) = (1 + rng.below(7), 1 + rng.below(7));
        let pairs: Vec<_> = (0..rng.below(left * right))
            .map(|_| (rng.below(left), rng.below(right)))
            .collect();
        check_bipartite_matching(left, right, &pairs, seed);
    }

    let complete =
        BipartiteGraph::from_edges(6, 6, (0..6).flat_map(|u| (0..6).map(move |v| (u, v)))).unwrap();
    let expected: Vec<_> = complete.randomized_maximum_matching(0).pairs().collect();
    assert_eq!(
        complete
            .randomized_maximum_matching(0)
            .pairs()
            .collect::<Vec<_>>(),
        expected
    );
    assert!((1..10).any(|seed| {
        !complete
            .randomized_maximum_matching(seed)
            .pairs()
            .eq(expected.iter().copied())
    }));
}

fn check_bipartite_matching(left: usize, right: usize, pairs: &[(usize, usize)], seed: u64) {
    let graph = BipartiteGraph::from_edges(left, right, pairs.iter().copied()).unwrap();
    let shifted: Vec<_> = pairs.iter().map(|&(u, v)| (u, left + v)).collect();
    let expected = expected_matching(left + right, &shifted);
    for matching in [
        graph.maximum_matching(),
        graph.randomized_maximum_matching(seed),
        graph.maximum_matching_kuhn(),
        graph.random_walk_maximum_matching(seed),
    ] {
        assert_eq!(matching.len(), expected, "{pairs:?}");
        assert_eq!(matching.pairs().count(), expected);
        for (u, v) in matching.pairs() {
            assert!(pairs.contains(&(u, v)));
            assert_eq!(matching.left_of(v), Some(u));
            assert_eq!(matching.right_of(u), Some(v));
        }
        let cover = matching.minimum_vertex_cover();
        assert_eq!(cover.left.len() + cover.right.len(), expected);
        for &(u, v) in pairs {
            assert!(cover.left.contains(&u) || cover.right.contains(&v));
        }
    }
}

fn cut_capacity(graph: &Graph<u64>, side: impl Fn(usize) -> bool) -> u128 {
    graph
        .edges()
        .iter()
        .filter(|edge| match graph.kind() {
            GraphKind::Directed => side(edge.from) && !side(edge.to),
            GraphKind::Undirected => side(edge.from) != side(edge.to),
        })
        .map(|edge| u128::from(edge.weight))
        .sum()
}

#[test]
fn flow_and_min_cut_agree_with_exhaustive_cuts() {
    let mut rng = Samples(17_152_319);
    for kind in [GraphKind::Directed, GraphKind::Undirected] {
        for case in 0..600 {
            let n = 2 + rng.below(8);
            let count = rng.below(n * n * 2);
            let capacities = [0, 1, 2, 3, 4, 7, 200, u64::MAX];
            let graph = graph(
                n,
                kind,
                (0..count).map(|_| {
                    (
                        rng.below(n),
                        rng.below(n),
                        capacities[rng.below(capacities.len())],
                    )
                }),
            );
            let source = rng.below(n);
            let sink = (source + 1 + rng.below(n - 1)) % n;
            let expected = (1..(1 << n) - 1)
                .filter(|&mask| {
                    kind == GraphKind::Undirected
                        || (mask & (1 << source) != 0 && mask & (1 << sink) == 0)
                })
                .map(|mask| cut_capacity(&graph, |node| mask & (1 << node) != 0))
                .min()
                .unwrap();
            if kind == GraphKind::Directed {
                let flows = [
                    max_flow,
                    max_flow_capacity_scaling,
                    max_flow_dinic,
                    max_flow_edmonds_karp,
                    max_flow_boykov_kolmogorov,
                ]
                .into_iter()
                .map(|solve| solve(&graph, source, sink).unwrap())
                .chain(std::iter::once_with(|| {
                    max_flow_network_simplex(&graph, source, sink).unwrap()
                }));
                for flow in flows {
                    assert_eq!(flow.value, expected, "directed case {case}");
                    assert!(flow.source_side[source] && !flow.source_side[sink]);
                    assert_eq!(
                        flow.value,
                        cut_capacity(&graph, |node| flow.source_side[node])
                    );
                    let mut balances = vec![0_i128; n];
                    for (&sent, edge) in flow.edge_flows.iter().zip(graph.edges()) {
                        assert!(sent <= edge.weight);
                        balances[edge.from] += i128::from(sent);
                        balances[edge.to] -= i128::from(sent);
                    }
                    for (node, &balance) in balances.iter().enumerate() {
                        let outward = if node == source {
                            flow.value as i128
                        } else if node == sink {
                            -(flow.value as i128)
                        } else {
                            0
                        };
                        assert_eq!(balance, outward);
                    }
                }
            } else {
                let cut = global_min_cut(&graph).unwrap().unwrap();
                assert_eq!(cut.value, expected, "undirected case {case}");
                assert!(cut.side.iter().any(|&set| set) && cut.side.iter().any(|&set| !set));
                assert_eq!(cut.value, cut_capacity(&graph, |node| cut.side[node]));
            }
        }
    }
}

#[test]
fn connectivity_and_spanning_forest_respect_graph_semantics() {
    let directed = graph(
        5,
        GraphKind::Directed,
        [(0, 1, ()), (1, 0, ()), (1, 2, ()), (2, 3, ()), (3, 2, ())],
    );
    let cc = strongly_connected_components(&directed);
    assert_eq!(cc.count, 3);
    assert_eq!(cc.component_of[0], cc.component_of[1]);
    assert_eq!(cc.component_of[2], cc.component_of[3]);
    assert!(cc.component_of[1] < cc.component_of[2]);
    assert_eq!(topological_sort(&directed), Err(GraphError::Cycle));

    let dag = graph(3, GraphKind::Directed, [(2, 0, ()), (2, 0, ()), (0, 1, ())]);
    assert_eq!(topological_sort(&dag).unwrap(), [2, 0, 1]);

    let weighted = graph(
        5,
        GraphKind::Undirected,
        [(1, 1, 0), (0, 1, 3), (1, 2, 2), (0, 1, 1), (3, 4, 4)],
    );
    let forest = minimum_spanning_forest(&weighted).unwrap();
    assert_eq!(forest.total_weight, 7);
    assert_eq!(forest.edges, [3, 2, 4]);
    assert_eq!(forest.component_count, 2);
}

#[test]
fn long_paths_use_iterative_searches() {
    let n = 50_000;
    let chain = graph(
        n,
        GraphKind::Directed,
        (1..n).map(|node| (node - 1, node, 1_u64)),
    );
    assert_eq!(dfs(&chain, 0).unwrap().len(), n);
    assert_eq!(strongly_connected_components(&chain).count, n);
    for solve in [
        max_flow,
        max_flow_dinic,
        max_flow_edmonds_karp,
        max_flow_boykov_kolmogorov,
    ] {
        assert_eq!(solve(&chain, 0, n - 1).unwrap().value, 1);
    }
    assert_eq!(max_flow_ford_fulkerson(&chain, 0, n - 1).unwrap().value, 1);
    // The greedy pass leaves an augmenting path that spans both partitions.
    let edges = (0..n - 1)
        .flat_map(|i| [(i, i), (i, i + 1)])
        .chain([(n - 1, 0)]);
    assert_eq!(
        BipartiteGraph::from_edges(n, n, edges)
            .unwrap()
            .maximum_matching()
            .len(),
        n
    );
}

#[cfg(feature = "multicommodity")]
mod fractional {
    use super::*;
    use graph_library::multicommodity::{Commodity, ConcurrentFlowError, maximum_concurrent_flow};

    #[test]
    fn optimizes_joint_routing_and_conserves_fractional_allocations() {
        let graph = graph(
            5,
            GraphKind::Directed,
            [(0, 1, 6), (0, 2, 5), (1, 4, 6), (2, 4, 5)],
        );
        let demands = [
            Commodity {
                source: 0,
                sink: 4,
                demand: 4,
            },
            Commodity {
                source: 1,
                sink: 4,
                demand: 4,
            },
        ];
        let result = maximum_concurrent_flow(&graph, &demands).unwrap();
        assert!((result.fraction - 1.375).abs() < 1e-8, "{result:?}");
        let amount = |commodity: usize, edge: usize| {
            result.commodity_flows[commodity]
                .iter()
                .find(|flow| flow.edge == edge)
                .unwrap()
                .amount
        };
        assert!((amount(0, 2) - 0.5).abs() < 1e-8, "{result:?}");
        assert!((amount(0, 3) - 5.0).abs() < 1e-8);
        assert!((amount(1, 2) - 5.5).abs() < 1e-8);

        let unreachable = [
            demands[0],
            Commodity {
                source: 3,
                sink: 4,
                demand: 1,
            },
        ];
        let zero = maximum_concurrent_flow(&graph, &unreachable).unwrap();
        assert_eq!(zero.fraction, 0.0);
        assert!(zero.commodity_flows.iter().all(Vec::is_empty));
        assert!(matches!(
            maximum_concurrent_flow(
                &graph,
                &[Commodity {
                    demand: 0,
                    ..demands[0]
                }]
            ),
            Err(ConcurrentFlowError::ZeroDemand { commodity: 0 })
        ));
    }

    #[test]
    fn shared_terminals_reduce_to_single_commodity_flow() {
        let mut rng = Samples(19_923_791);
        for case in 0..75 {
            let n = 2 + rng.below(8);
            let count = rng.below(n * n * 2);
            let graph = graph(
                n,
                GraphKind::Directed,
                (0..count).map(|_| (rng.below(n), rng.below(n), rng.below(30) as u64)),
            );
            let source = 0;
            let sink = n - 1;
            let demands = [1, 7, 23].map(|demand| Commodity {
                source,
                sink,
                demand,
            });
            let flow = max_flow(&graph, source, sink).unwrap();
            let concurrent = maximum_concurrent_flow(&graph, &demands).unwrap_or_else(|error| {
                panic!("case {case} ({n} vertices), {:?}: {error}", graph.edges())
            });
            assert!(
                (concurrent.fraction - flow.value as f64 / 31.0).abs() < 1e-7,
                "case {case}: {concurrent:?}"
            );
        }
    }
}

#[test]
fn bridges_and_cut_vertices_agree_with_deletion_on_multigraphs() {
    let mut rng = Samples(61_922);
    for _ in 0..180 {
        let n = 1 + rng.below(8);
        let edges: Vec<_> = (0..rng.below(4 * n))
            .map(|_| (rng.below(n), rng.below(n), ()))
            .collect();
        let full = graph(n, GraphKind::Undirected, edges.iter().copied());
        let full_count = weakly_connected_components(&full).count;
        let cuts = bridges_and_articulation_points(&full).unwrap();
        let expected_bridges: Vec<_> = (0..edges.len())
            .filter(|&removed| {
                let remaining = edges
                    .iter()
                    .enumerate()
                    .filter_map(|(id, edge)| (id != removed).then_some(*edge));
                weakly_connected_components(&graph(n, GraphKind::Undirected, remaining)).count
                    > full_count
            })
            .collect();
        let expected_points: Vec<_> = (0..n)
            .filter(|&removed| {
                let remaining = edges
                    .iter()
                    .copied()
                    .filter(|&(a, b, ())| a != removed && b != removed);
                // The removed vertex is still represented as an isolated slot.
                weakly_connected_components(&graph(n, GraphKind::Undirected, remaining)).count
                    > full_count + 1
            })
            .collect();
        assert_eq!(cuts.bridges, expected_bridges);
        assert_eq!(cuts.articulation_points, expected_points);
    }
    let disconnected = graph(6, GraphKind::Directed, [(2, 0, ()), (2, 1, ()), (3, 4, ())]);
    let weak = weakly_connected_components(&disconnected);
    assert_eq!(weak.component_of, [0, 0, 0, 1, 1, 2]);
    assert_eq!(weak.count, 3);
    assert_eq!(
        bridges_and_articulation_points(&disconnected),
        Err(GraphError::RequiresUndirected)
    );
}

#[test]
fn core_numbers_agree_with_repeated_degree_peeling() {
    let mut rng = Samples(6_191);
    for kind in [GraphKind::Directed, GraphKind::Undirected] {
        for _ in 0..120 {
            let n = 1 + rng.below(10);
            let edges: Vec<_> = (0..rng.below(5 * n))
                .map(|_| (rng.below(n), rng.below(n), ()))
                .collect();
            let graph = graph(n, kind, edges.iter().copied());
            let mut expected = vec![0; n];
            for k in 1..=edges.len() {
                let mut active = vec![true; n];
                loop {
                    let mut degree = vec![0; n];
                    for &(a, b, ()) in &edges {
                        if a != b && active[a] && active[b] {
                            degree[a] += 1;
                            degree[b] += 1;
                        }
                    }
                    let last_size = active.iter().filter(|&&included| included).count();
                    for (node, included) in active.iter_mut().enumerate() {
                        *included &= degree[node] >= k;
                    }
                    if active.iter().filter(|&&included| included).count() == last_size {
                        break;
                    }
                }
                for (node, included) in active.into_iter().enumerate() {
                    if included {
                        expected[node] = k;
                    }
                }
            }
            assert_eq!(core_numbers(&graph), expected, "{kind:?} {edges:?}");
        }
    }
}

#[test]
fn personalized_weighted_pagerank_agrees_with_networkx() {
    let graph = graph(
        5,
        GraphKind::Directed,
        [
            (0, 1, 2),
            (0, 1, 1),
            (0, 2, 1),
            (1, 2, 4),
            (2, 0, 2),
            (2, 3, 1),
            (3, 3, 1),
            (4, 1, 0),
        ],
    );
    let options = PageRankOptions {
        personalization: Some(&[2.0, 0.0, 0.0, 1.0, 1.0]),
        tolerance: 1e-12,
        ..PageRankOptions::default()
    };
    let expected = [
        0.16634869806153663,
        0.10604729501423395,
        0.12548929910018186,
        0.5544956602049995,
        0.04761904761904763,
    ];
    let rank = pagerank(&graph, options).unwrap();
    for (actual, expected) in rank.iter().zip(expected) {
        assert!((actual - expected).abs() < 1e-10);
    }
    assert!((rank.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert!(matches!(
        pagerank(
            &graph,
            PageRankOptions {
                max_iterations: 1,
                ..options
            }
        ),
        Err(PageRankError::DidNotConverge { iterations: 1, .. })
    ));
    let uniform = pagerank(
        &graph,
        PageRankOptions {
            personalization: Some(&[f64::MAX; 5]),
            damping: 0.0,
            ..options
        },
    )
    .unwrap();
    assert_eq!(uniform, [0.2; 5]);
}

#[test]
fn gomory_hu_queries_and_edge_partitions_agree_with_original_cuts() {
    use graph_library::flow::gomory_hu_tree;

    let mut rng = Samples(241_163);
    for _ in 0..95 {
        let n = 2 + rng.below(7);
        let count = rng.below(n * n);
        let graph = graph(
            n,
            GraphKind::Undirected,
            (0..count).map(|_| {
                (
                    rng.below(n),
                    rng.below(n),
                    [0, 1, 2, 3, 7, u64::MAX][rng.below(6)],
                )
            }),
        );
        let index = gomory_hu_tree(&graph).unwrap();
        let tree = index.tree();
        assert_eq!(tree.edge_count(), n - 1);

        for a in 0..n {
            for b in a + 1..n {
                let expected = (1..(1 << n) - 1)
                    .filter(|&mask| mask & (1 << a) != 0 && mask & (1 << b) == 0)
                    .map(|mask| cut_capacity(&graph, |node| mask & (1 << node) != 0))
                    .min()
                    .unwrap();
                assert_eq!(index.minimum_cut_value(a, b).unwrap(), expected);
                assert_eq!(index.minimum_cut_value(b, a).unwrap(), expected);
            }
        }
        for (id, cut) in tree.edges().iter().enumerate() {
            let mut reachable = vec![false; n];
            let mut stack = vec![cut.from];
            reachable[cut.from] = true;
            while let Some(node) = stack.pop() {
                for arc in tree.outgoing(node).unwrap() {
                    if arc.edge != id && !reachable[arc.node] {
                        reachable[arc.node] = true;
                        stack.push(arc.node);
                    }
                }
            }
            assert_eq!(cut.weight, cut_capacity(&graph, |node| reachable[node]));
        }
    }
}
