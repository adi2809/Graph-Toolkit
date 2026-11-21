use graph_library::{Graph, flow::MaxFlow};

pub fn verify(g: &Graph<u64>, s: usize, t: usize, f: &MaxFlow) {
    let mut balance = vec![0i128; g.node_count()];
    let mut cut = 0u128;
    assert_eq!(f.edge_flows.len(), g.edge_count());
    assert_eq!(f.source_side.len(), g.node_count());
    assert!(f.source_side[s] && !f.source_side[t]);
    for (e, &used) in g.edges().iter().zip(&f.edge_flows) {
        assert!(used <= e.weight);
        balance[e.from] += i128::from(used);
        balance[e.to] -= i128::from(used);
        if f.source_side[e.from] && !f.source_side[e.to] {
            cut += u128::from(e.weight);
        }
    }
    assert_eq!(cut, f.value);
    assert_eq!(balance[s], i128::try_from(f.value).unwrap());
    assert_eq!(balance[t], -balance[s]);
    balance[s] = 0;
    balance[t] = 0;
    assert!(balance.iter().all(|&x| x == 0));
    let mut seen = vec![false; g.node_count()];
    let mut stack = vec![s];
    seen[s] = true;
    while let Some(u) = stack.pop() {
        for a in g.outgoing(u).unwrap() {
            if !seen[a.node] && f.edge_flows[a.edge] < g.edges()[a.edge].weight {
                seen[a.node] = true;
                stack.push(a.node);
            }
        }
        for a in g.incoming(u).unwrap() {
            if !seen[a.node] && f.edge_flows[a.edge] > 0 {
                seen[a.node] = true;
                stack.push(a.node);
            }
        }
    }
    assert_eq!(seen, f.source_side);
}
