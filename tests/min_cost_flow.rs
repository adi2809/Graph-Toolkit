use graph_library::flow::{
    ArcCost, CostFlow, CostFlowError, flow_cost, min_cost_max_flow,
    min_cost_max_flow_successive_shortest_path, minimum_cost_flow,
    minimum_cost_flow_cycle_canceling, minimum_cost_flow_successive_shortest_path, network_simplex,
};
use graph_library::{Graph, GraphKind};

fn graph(n: usize, edges: &[(usize, usize, u64, i64)]) -> Graph<ArcCost> {
    let mut builder = Graph::builder(n, GraphKind::Directed);
    for &(from, to, capacity, unit_cost) in edges {
        builder
            .add_edge(
                from,
                to,
                ArcCost {
                    capacity,
                    unit_cost,
                },
            )
            .unwrap();
    }
    builder.build()
}

fn balances(graph: &Graph<ArcCost>, flow: &CostFlow) -> Vec<i128> {
    assert_eq!(flow_cost(graph, &flow.edge_flows).unwrap(), flow.total_cost);
    let mut balance = vec![0; graph.node_count()];
    for (edge, &amount) in graph.edges().iter().zip(&flow.edge_flows) {
        balance[edge.from] += i128::from(amount);
        balance[edge.to] -= i128::from(amount);
    }
    balance
}

fn exhaustive(
    graph: &Graph<ArcCost>,
    supply: &[i64],
    balance: &mut [i64],
    edge: usize,
    cost: i128,
    best: &mut Option<i128>,
    max_flow: &mut (u128, i128),
) {
    if let Some(arc) = graph.edge(edge) {
        for amount in 0..=arc.weight.capacity {
            let amount = amount as i64;
            balance[arc.from] += amount;
            balance[arc.to] -= amount;
            let next_cost = cost + i128::from(amount) * i128::from(arc.weight.unit_cost);
            exhaustive(graph, supply, balance, edge + 1, next_cost, best, max_flow);
            balance[arc.from] -= amount;
            balance[arc.to] += amount;
        }
        return;
    }
    if balance == supply {
        *best = Some(best.map_or(cost, |current| current.min(cost)));
    }
    let value = balance[0];
    if value >= 0
        && balance[1..balance.len() - 1].iter().all(|&x| x == 0)
        && balance[balance.len() - 1] == -value
    {
        let value = value as u128;
        if value > max_flow.0 || (value == max_flow.0 && cost < max_flow.1) {
            *max_flow = (value, cost);
        }
    }
}

#[test]
fn cost_scaling_agrees_with_exhaustive_integer_flows() {
    let mut seed = 40_169_u64;
    let mut sample = |n: usize| {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed as usize % n
    };
    for case in 0..280 {
        let n = 2 + sample(4);
        let edges: Vec<_> = (0..sample(9))
            .map(|_| {
                (
                    sample(n),
                    sample(n),
                    sample(3) as u64,
                    sample(13) as i64 - 6,
                )
            })
            .collect();
        let graph = graph(n, &edges);
        let mut supply = vec![0; n];
        if case % 3 == 0 {
            for &(a, b, capacity, _) in &edges {
                let amount = sample(capacity as usize + 1) as i64;
                supply[a] += amount;
                supply[b] -= amount;
            }
        } else {
            let amount = sample(4) as i64;
            supply[0] = amount;
            supply[n - 1] = -amount;
        }
        let mut best = None;
        let mut maximum = (0, 0);
        exhaustive(
            &graph,
            &supply,
            &mut vec![0; n],
            0,
            0,
            &mut best,
            &mut maximum,
        );
        let canceled = minimum_cost_flow_cycle_canceling(&graph, &supply).unwrap();
        assert_eq!(
            canceled.as_ref().map(|flow| flow.total_cost),
            best,
            "cycle cancel case {case}"
        );
        if let Some(result) = canceled {
            assert_eq!(
                balances(&graph, &result),
                supply.iter().copied().map(i128::from).collect::<Vec<_>>()
            );
        }
        let nonnegative = self::graph(
            n,
            &edges
                .iter()
                .map(|&(a, b, cap, cost)| (a, b, cap, cost.abs()))
                .collect::<Vec<_>>(),
        );
        let expected = minimum_cost_flow(&nonnegative, &supply).unwrap();
        let ssp = minimum_cost_flow_successive_shortest_path(&nonnegative, &supply).unwrap();
        assert_eq!(
            ssp.as_ref().map(|flow| flow.total_cost),
            expected.as_ref().map(|flow| flow.total_cost)
        );
        if let Some(result) = ssp {
            assert_eq!(
                balances(&nonnegative, &result),
                supply.iter().copied().map(i128::from).collect::<Vec<_>>()
            );
        }
        let expected = min_cost_max_flow(&nonnegative, 0, n - 1).unwrap();
        let ssp = min_cost_max_flow_successive_shortest_path(&nonnegative, 0, n - 1).unwrap();
        assert_eq!(
            (ssp.value, ssp.flow.total_cost),
            (expected.value, expected.flow.total_cost)
        );

        let simplex = network_simplex(&graph, &supply).unwrap();
        assert_eq!(
            simplex.as_ref().map(|flow| flow.total_cost),
            best,
            "simplex case {case}"
        );
        if let Some(result) = simplex {
            assert_eq!(
                balances(&graph, &result),
                supply.iter().copied().map(i128::from).collect::<Vec<_>>()
            );
        }
        let result = minimum_cost_flow(&graph, &supply).unwrap();
        assert_eq!(
            result.as_ref().map(|flow| flow.total_cost),
            best,
            "case {case}"
        );
        if let Some(result) = result {
            assert_eq!(
                balances(&graph, &result),
                supply.iter().copied().map(i128::from).collect::<Vec<_>>()
            );
        }

        let solved = min_cost_max_flow(&graph, 0, n - 1).unwrap();
        assert_eq!(
            (solved.value, solved.flow.total_cost),
            maximum,
            "case {case}"
        );
        let mut expected = vec![0; n];
        expected[0] = solved.value as i128;
        expected[n - 1] = -(solved.value as i128);
        assert_eq!(balances(&graph, &solved.flow), expected);
    }
}

#[test]
fn disconnected_negative_cycles_and_wide_arithmetic() {
    let network = graph(
        5,
        &[
            (0, 1, 3, 5),
            (0, 1, 5, -2),
            (2, 3, 3, -4),
            (3, 2, 2, 1),
            (4, 4, 4, -1),
        ],
    );
    let result = minimum_cost_flow(&network, &[6, -6, 0, 0, 0])
        .unwrap()
        .unwrap();
    assert_eq!(result.edge_flows, [1, 5, 2, 2, 4]);
    assert_eq!(result.total_cost, -15);
    assert!(
        minimum_cost_flow(&network, &[9, -9, 0, 0, 0])
            .unwrap()
            .is_none()
    );

    let wide = graph(2, &[(0, 1, u64::MAX, i64::MAX), (0, 1, u64::MAX, i64::MIN)]);
    let result = min_cost_max_flow(&wide, 0, 1).unwrap();
    assert_eq!(result.value, 2 * u128::from(u64::MAX));
    assert_eq!(result.flow.total_cost, -i128::from(u64::MAX));

    let overflow = graph(1, &[(0, 0, u64::MAX, i64::MIN), (0, 0, u64::MAX, i64::MIN)]);
    assert_eq!(
        minimum_cost_flow(&overflow, &[0]),
        Err(CostFlowError::ArithmeticOverflow)
    );
}
