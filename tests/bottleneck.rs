use graph_library::connectivity::{
    minimum_bottleneck_arborescence, minimum_bottleneck_spanning_tree, minimum_spanning_tree,
};
use graph_library::matching::{
    is_matching, minimum_bottleneck_matching, minimum_bottleneck_matching_of_size,
    minimum_bottleneck_perfect_matching,
};
use graph_library::shortest_path::minimum_bottleneck_path;
use graph_library::{Graph, GraphError, GraphKind};

fn graph<W>(
    n: usize,
    kind: GraphKind,
    edges: impl IntoIterator<Item = (usize, usize, W)>,
) -> Graph<W> {
    let mut builder = Graph::builder(n, kind);
    for (a, b, cost) in edges {
        builder.add_edge(a, b, cost).unwrap();
    }
    builder.build()
}

struct Samples(u64);

impl Samples {
    fn below(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 as usize) % n
    }
}

fn reached(graph: &Graph<i64>, source: usize, mut include: impl FnMut(usize) -> bool) -> Vec<bool> {
    let mut visited = vec![false; graph.node_count()];
    let mut queue = vec![source];
    visited[source] = true;
    let mut head = 0;
    while let Some(&node) = queue.get(head) {
        head += 1;
        for arc in graph.outgoing(node).unwrap() {
            if !visited[arc.node] && include(arc.edge) {
                visited[arc.node] = true;
                queue.push(arc.node);
            }
        }
    }
    visited
}

#[test]
fn minimax_paths_and_rooted_trees_match_exhaustive_thresholds() {
    let mut rng = Samples(91_417);
    for kind in [GraphKind::Directed, GraphKind::Undirected] {
        for _ in 0..100 {
            let n = 1 + rng.below(8);
            let count = rng.below(n * 3);
            let costs = [i64::MIN, -12, -1, 0, 4, 8, i64::MAX];
            let graph = graph(
                n,
                kind,
                (0..count).map(|_| (rng.below(n), rng.below(n), costs[rng.below(costs.len())])),
            );
            let mut thresholds: Vec<_> = graph.edges().iter().map(|edge| edge.weight).collect();
            thresholds.sort_unstable();
            thresholds.dedup();
            for root in 0..n {
                let reachable: Vec<_> = thresholds
                    .iter()
                    .map(|&threshold| {
                        reached(&graph, root, |id| graph.edges()[id].weight <= threshold)
                    })
                    .collect();
                for dst in 0..n {
                    let result = minimum_bottleneck_path(&graph, root, dst).unwrap();
                    if root == dst {
                        let path = result.unwrap();
                        assert_eq!(path.bottleneck, None);
                        assert_eq!(path.nodes, [root]);
                        continue;
                    }
                    let expected = reachable
                        .iter()
                        .position(|row| row[dst])
                        .map(|pos| thresholds[pos]);
                    assert_eq!(result.as_ref().and_then(|path| path.bottleneck), expected);
                    if let Some(path) = result {
                        assert_eq!(path.nodes.first(), Some(&root));
                        assert_eq!(path.nodes.last(), Some(&dst));
                        assert_eq!(path.edges.len() + 1, path.nodes.len());
                        for (nodes, &id) in path.nodes.windows(2).zip(&path.edges) {
                            let edge = graph.edges()[id];
                            assert!(
                                (edge.from, edge.to) == (nodes[0], nodes[1])
                                    || (kind == GraphKind::Undirected
                                        && (edge.to, edge.from) == (nodes[0], nodes[1]))
                            );
                        }
                        assert_eq!(
                            path.edges.iter().map(|&id| graph.edges()[id].weight).max(),
                            expected
                        );
                    }
                }
                if kind == GraphKind::Directed {
                    let tree = minimum_bottleneck_arborescence(&graph, root).unwrap();
                    let expected = reachable
                        .iter()
                        .position(|row| row.iter().all(|&v| v))
                        .map(|pos| thresholds[pos]);
                    assert_eq!(tree.is_some(), n == 1 || expected.is_some());
                    if let Some(tree) = tree {
                        assert_eq!(tree.bottleneck, if n == 1 { None } else { expected });
                        assert_eq!(tree.edges.len(), n - 1);
                        assert!(
                            reached(&graph, root, |id| tree.edges.contains(&id))
                                .into_iter()
                                .all(|v| v)
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn spanning_tree_optima_agree_with_enumerated_edge_sets() {
    let mut rng = Samples(233_791);
    for _ in 0..175 {
        let n = 1 + rng.below(6);
        let count = rng.below(12);
        let graph = graph(
            n,
            GraphKind::Undirected,
            (0..count).map(|_| (rng.below(n), rng.below(n), rng.below(13) as i64 - 6)),
        );
        let mut optimal_total = None;
        let mut optimal_bottleneck = None;
        for mask in 0..1_usize << count {
            let contains = |id: usize| mask & (1 << id) != 0;
            if mask.count_ones() as usize != n - 1
                || !reached(&graph, 0, contains).into_iter().all(|v| v)
            {
                continue;
            }
            let weights = (0..count)
                .filter(|&id| contains(id))
                .map(|id| graph.edges()[id].weight);
            let total: i128 = weights.clone().map(i128::from).sum();
            optimal_total = Some(optimal_total.map_or(total, |old: i128| old.min(total)));
            if let Some(bottleneck) = weights.max() {
                optimal_bottleneck =
                    Some(optimal_bottleneck.map_or(bottleneck, |old: i64| old.min(bottleneck)));
            }
        }
        let spanning = minimum_spanning_tree(&graph).unwrap();
        let bottleneck = minimum_bottleneck_spanning_tree(&graph).unwrap();
        assert_eq!(
            spanning.as_ref().map(|tree| tree.total_weight),
            optimal_total
        );
        assert_eq!(bottleneck.is_some(), spanning.is_some());
        for (ids, reported) in spanning
            .iter()
            .map(|tree| (&tree.edges, tree.bottleneck))
            .chain(bottleneck.iter().map(|tree| (&tree.edges, tree.bottleneck)))
        {
            assert_eq!(ids.len(), n - 1);
            assert_eq!(reported, optimal_bottleneck);
            assert_eq!(
                ids.iter().map(|&id| graph.edges()[id].weight).max(),
                reported
            );
            assert!(
                reached(&graph, 0, |id| ids.contains(&id))
                    .into_iter()
                    .all(|v| v)
            );
        }
    }

    let empty = graph::<i64>(0, GraphKind::Undirected, []);
    assert!(
        minimum_bottleneck_spanning_tree(&empty)
            .unwrap()
            .unwrap()
            .edges
            .is_empty()
    );
    assert_eq!(
        minimum_spanning_tree(&empty).unwrap().unwrap().total_weight,
        0
    );
    let extremes = graph(
        5,
        GraphKind::Undirected,
        [
            (0, 1, i128::MIN),
            (1, 2, i128::MIN),
            (2, 3, i128::MAX),
            (3, 4, i128::MAX),
        ],
    );
    assert_eq!(
        minimum_spanning_tree(&extremes)
            .unwrap()
            .unwrap()
            .total_weight,
        -2
    );
    let overflow = graph(3, GraphKind::Undirected, [(0, 1, i128::MIN), (1, 2, -1)]);
    assert_eq!(
        minimum_spanning_tree(&overflow),
        Err(GraphError::WeightOverflow)
    );
}

#[test]
fn bottleneck_matchings_agree_with_every_feasible_edge_subset() {
    let mut rng = Samples(5_109);
    for _ in 0..230 {
        let n = 1 + rng.below(9);
        let count = rng.below(13);
        let costs = [i64::MIN, -8, -1, 0, 2, 3, i64::MAX];
        let graph = graph(
            n,
            GraphKind::Undirected,
            (0..count).map(|_| (rng.below(n), rng.below(n), costs[rng.below(costs.len())])),
        );
        let mut optimal: Vec<Option<i64>> = vec![None; n / 2 + 1];
        for mask in 1..1_usize << count {
            let ids: Vec<_> = (0..count).filter(|&id| mask & (1 << id) != 0).collect();
            if is_matching(&graph, &ids).unwrap() {
                let max = ids
                    .iter()
                    .map(|&id| graph.edges()[id].weight)
                    .max()
                    .unwrap();
                let value = &mut optimal[ids.len()];
                *value = Some(value.map_or(max, |old| old.min(max)));
            }
        }
        for (size, &expected) in optimal.iter().enumerate() {
            let result = minimum_bottleneck_matching_of_size(&graph, size).unwrap();
            assert_eq!(result.is_some(), size == 0 || expected.is_some());
            if let Some(result) = result {
                assert_eq!(result.matching.len(), size);
                assert_eq!(result.bottleneck, expected);
                let ids: Vec<_> = result.matching.edges().collect();
                assert!(is_matching(&graph, &ids).unwrap());
                assert_eq!(
                    ids.iter().map(|&id| graph.edges()[id].weight).max(),
                    expected
                );
            }
        }
        let size = optimal.iter().rposition(Option::is_some).unwrap_or(0);
        let result = minimum_bottleneck_matching(&graph).unwrap();
        assert_eq!(
            (result.matching.len(), result.bottleneck),
            (size, optimal[size])
        );
        assert_eq!(
            minimum_bottleneck_perfect_matching(&graph)
                .unwrap()
                .is_some(),
            n == 2 * size
        );
    }
}
