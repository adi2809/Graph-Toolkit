//! Exact integer reference for static cost-to-length cycle ratios.
//! Parametric search replaces c by q*c - p*length. A negative cycle
//! improves p/q strictly. The finite set of simple-cycle ratios terminates
//! the iteration. Each decision costs O(VE), with linear working storage.
use crate::{Graph, GraphError, GraphKind, NONE};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CycleRatioWeight {
    pub cost: i64,
    pub length: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct CycleStep {
    pub edge: usize,
    pub reversed: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RatioCycle {
    /// Original edge IDs in cyclic traversal order.
    pub steps: Vec<CycleStep>,
    pub total_cost: i128,
    pub total_length: u128,
}

#[derive(Clone, Copy)]
pub(super) struct RatioArc {
    pub u: usize,
    pub v: usize,
    pub cost: i128,
    pub length: u64,
}

fn transformed_cost(a: RatioArc, p: i128, q: i128) -> Result<i128, GraphError> {
    a.cost
        .checked_mul(q)
        .and_then(|x| x.checked_sub(p.checked_mul(i128::from(a.length))?))
        .ok_or(GraphError::WeightOverflow)
}

// With equality enabled, a symbolic negative infinitesimal per edge makes
// an exactly tight directed cycle detectable even when every cost is zero.
pub(super) fn cycle_at_threshold(
    n: usize,
    arcs: &[RatioArc],
    p: i128,
    q: i128,
    include_equal: bool,
) -> Result<Option<Vec<usize>>, GraphError> {
    if n == 0 || arcs.is_empty() {
        return Ok(None);
    }
    let costs: Vec<_> = arcs
        .iter()
        .map(|&a| transformed_cost(a, p, q))
        .collect::<Result<_, _>>()?;
    let mut distance = vec![(0i128, 0i128); n];
    let mut parent = vec![NONE; n];
    let mut changed = NONE;
    for _ in 0..n {
        changed = NONE;
        for (i, (a, &cost)) in arcs.iter().zip(&costs).enumerate() {
            let shifted = distance[a.u]
                .0
                .checked_add(cost)
                .ok_or(GraphError::WeightOverflow)?;
            let d = (shifted, distance[a.u].1 - i128::from(include_equal));
            if d < distance[a.v] {
                distance[a.v] = d;
                parent[a.v] = i;
                changed = a.v;
            }
        }
        if changed == NONE {
            return Ok(None);
        }
    }
    for _ in 0..n {
        changed = arcs[parent[changed]].u;
    }
    let first = changed;
    let mut cycle = Vec::new();
    loop {
        let i = parent[changed];
        cycle.push(i);
        changed = arcs[i].u;
        if changed == first {
            break;
        }
    }
    cycle.reverse();
    Ok(Some(cycle))
}

pub(super) fn best_below(
    n: usize,
    arcs: &[RatioArc],
    mut numerator: i128,
) -> Result<Option<Vec<usize>>, GraphError> {
    let mut denominator = 1;
    let mut best = None;
    while let Some(cycle) = cycle_at_threshold(n, arcs, numerator, denominator, false)? {
        let cost = cycle.iter().map(|&i| arcs[i].cost).sum::<i128>();
        let length = cycle
            .iter()
            .map(|&i| i128::from(arcs[i].length))
            .sum::<i128>();
        let (mut a, mut b) = (cost.unsigned_abs(), length as u128);
        while b != 0 {
            (a, b) = (b, a % b);
        }
        numerator = cost / a as i128;
        denominator = length / a as i128;
        best = Some(cycle);
    }
    Ok(best)
}

/// Smallest sum(cost)/sum(length) over directed cycles, including loops.
/// All lengths must be positive. A DAG returns None. Arithmetic exhaustion
/// in a parametric decision returns WeightOverflow.
pub fn minimum_ratio_cycle(
    graph: &Graph<CycleRatioWeight>,
) -> Result<Option<RatioCycle>, GraphError> {
    graph.require(GraphKind::Directed)?;
    ratio_cycle(graph, false)
}

/// Static undirected min-ratio circulation oracle with oriented gradients.
/// Traversing an input edge backwards negates its cost and preserves length.
/// Returns None when every circulation has zero gradient (including forests).
pub fn minimum_undirected_ratio_cycle(
    graph: &Graph<CycleRatioWeight>,
) -> Result<Option<RatioCycle>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    ratio_cycle(graph, true)
}

fn ratio_cycle(
    graph: &Graph<CycleRatioWeight>,
    undirected: bool,
) -> Result<Option<RatioCycle>, GraphError> {
    let mut arcs = Vec::new();
    let mut upper = 0;
    for e in graph.edges() {
        if e.weight.length == 0 {
            return Err(GraphError::InvalidModelParameters);
        }
        let cost = i128::from(e.weight.cost);
        upper = upper.max(cost + 1);
        arcs.push(RatioArc {
            u: e.from,
            v: e.to,
            cost,
            length: e.weight.length,
        });
        if undirected {
            arcs.push(RatioArc {
                u: e.to,
                v: e.from,
                cost: -cost,
                length: e.weight.length,
            });
        }
    }
    Ok(best_below(
        graph.node_count(),
        &arcs,
        if undirected { 0 } else { upper },
    )?
    .map(|cycle| RatioCycle {
        total_cost: cycle.iter().map(|&i| arcs[i].cost).sum(),
        total_length: cycle.iter().map(|&i| u128::from(arcs[i].length)).sum(),
        steps: cycle
            .into_iter()
            .map(|i| CycleStep {
                edge: if undirected { i / 2 } else { i },
                reversed: undirected && i % 2 == 1,
            })
            .collect(),
    }))
}
