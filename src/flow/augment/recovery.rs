//! Deterministic bit-by-bit rounding of a feasible dyadic flow.
use super::{MaxFlow, Network, validate_terminals};
use crate::{Graph, GraphError, NONE};

/// Rounds a feasible s–t flow with common denominator 2^bits. Every edge
/// rounds to its floor or ceiling, and total throughput rounds upwards.
/// Returns Some only if the resulting integral flow has a residual cut.
/// An invalid capacity, balance, numerator count or precision is an error.
/// Time O((V+E)(1+bits)), space O(V+E). Supports up to 64 fractional bits.
pub fn recover_dyadic_max_flow(
    graph: &Graph<u64>,
    source: usize,
    sink: usize,
    numerators: &[u128],
    bits: u32,
) -> Result<Option<MaxFlow>, GraphError> {
    validate_terminals(graph, source, sink)?;
    if bits > 64 || numerators.len() != graph.edge_count() {
        return Err(GraphError::InvalidModelParameters);
    }
    let unit = 1u128 << bits;
    // Keep the integer and fractional limbs separate. Aggregate throughput
    // can exceed u128 when expressed in units of the least significant bit.
    let mut balance = vec![(0i128, 0i128); graph.node_count()];
    for (edge, &x) in graph.edges().iter().zip(numerators) {
        if x > u128::from(edge.weight) * unit {
            return Err(GraphError::InvalidModelParameters);
        }
        let (q, r) = ((x >> bits) as i128, (x % unit) as i128);
        balance[edge.from].0 += q;
        balance[edge.from].1 += r;
        balance[edge.to].0 -= q;
        balance[edge.to].1 -= r;
    }
    for (v, (q, r)) in balance.into_iter().enumerate() {
        let whole = q + r.div_euclid(unit as i128);
        let fraction = r.rem_euclid(unit as i128);
        if (v == source && whole < 0) || (v != source && v != sink && (whole != 0 || fraction != 0))
        {
            return Err(GraphError::InvalidModelParameters);
        }
    }
    Ok(round_and_certify(
        graph,
        source,
        sink,
        numerators.to_vec(),
        bits,
    ))
}

pub(in crate::flow) fn round_and_certify(
    graph: &Graph<u64>,
    s: usize,
    t: usize,
    mut amounts: Vec<u128>,
    bits: u32,
) -> Option<MaxFlow> {
    let unit = 1u128 << bits;
    let mut fractional: Vec<_> = amounts
        .iter()
        .enumerate()
        .filter(|(_, x)| **x % unit != 0)
        .map(|(i, _)| i)
        .collect();
    let mut head = vec![NONE; graph.node_count()];
    let mut darts = Vec::new();
    let mut used = vec![false; graph.edge_count()];
    let mut touched = Vec::new();
    for k in 0..bits {
        if fractional.is_empty() {
            break;
        }
        let bit = 1u128 << k;
        darts.clear();
        touched.clear();
        let mut odd_source = false;
        for &id in &fractional {
            if amounts[id] & bit == 0 {
                continue;
            }
            let edge = &graph.edges()[id];
            used[id] = false;
            for u in [edge.from, edge.to] {
                if head[u] == NONE {
                    touched.push(u);
                }
                odd_source ^= u == s;
                let previous = std::mem::replace(&mut head[u], darts.len());
                darts.push((id, previous));
            }
        }
        // At this precision the odd-edge subgraph is Eulerian, except
        // possibly at s and t. Orient its open trail from s to t first.
        // All remaining maximal trails are closed and preserve throughput.
        let mut walk = |mut u: usize| {
            while head[u] != NONE {
                let (id, next) = darts[head[u]];
                head[u] = next;
                if std::mem::replace(&mut used[id], true) {
                    continue;
                }
                let e = &graph.edges()[id];
                if u == e.from {
                    amounts[id] += bit;
                    u = e.to;
                } else {
                    amounts[id] -= bit;
                    u = e.from;
                }
            }
            u
        };
        if odd_source {
            let end = walk(s);
            debug_assert_eq!(end, t);
        }
        for &start in &touched {
            let end = walk(start);
            debug_assert_eq!(end, start);
        }
        fractional.retain(|&id| amounts[id] % unit != 0);
    }
    let mut state = Network::new(graph);
    let mut incoming = 0u128;
    for ((e, pair), x) in graph.edges().iter().zip(&mut state.capacity).zip(amounts) {
        let f = (x >> bits) as u64;
        *pair = [e.weight - f, f];
        if e.from == s {
            state.value += u128::from(f);
        }
        if e.to == s {
            incoming += u128::from(f);
        }
    }
    state.value -= incoming;
    let flow = state.finish(s);
    (!flow.source_side[t]).then_some(flow)
}
