//! Karp's minimum-mean cycle formula, evaluated independently on each SCC.
use super::ratio::{CycleStep, RatioArc, RatioCycle, cycle_at_threshold};
use crate::{Graph, GraphError, GraphKind, connectivity::strongly_connected_components};

#[derive(Clone, Copy)]
struct Fraction(i128, i128);
impl Fraction {
    fn less(self, other: Self) -> Result<bool, GraphError> {
        let left = self
            .0
            .checked_mul(other.1)
            .ok_or(GraphError::WeightOverflow)?;
        let right = other
            .0
            .checked_mul(self.1)
            .ok_or(GraphError::WeightOverflow)?;
        Ok(left < right)
    }
}

fn advance(
    arcs: &[RatioArc],
    dist: &mut Vec<i128>,
    next: &mut Vec<i128>,
) -> Result<(), GraphError> {
    next.fill(i128::MAX);
    for a in arcs {
        if dist[a.u] != i128::MAX {
            let candidate = dist[a.u]
                .checked_add(a.cost)
                .ok_or(GraphError::WeightOverflow)?;
            next[a.v] = next[a.v].min(candidate);
        }
    }
    std::mem::swap(dist, next);
    Ok(())
}

// F_k(v) is the cheapest walk of exactly k steps from a fixed start.
// Replay the recurrence after computing F_n, retaining just three rows.
fn component_mean(n: usize, arcs: &[RatioArc]) -> Result<Fraction, GraphError> {
    let mut dist = vec![i128::MAX; n];
    dist[0] = 0;
    let mut next = dist.clone();
    for _ in 0..n {
        advance(arcs, &mut dist, &mut next)?;
    }
    let last = dist.clone();
    dist.fill(i128::MAX);
    dist[0] = 0;
    let mut maxima = vec![None::<Fraction>; n];
    for k in 0..n {
        for v in 0..n {
            if last[v] == i128::MAX || dist[v] == i128::MAX {
                continue;
            }
            let delta = last[v]
                .checked_sub(dist[v])
                .ok_or(GraphError::WeightOverflow)?;
            let value = Fraction(delta, (n - k) as i128);
            if match maxima[v] {
                Some(old) => old.less(value)?,
                None => true,
            } {
                maxima[v] = Some(value);
            }
        }
        if k + 1 < n {
            advance(arcs, &mut dist, &mut next)?;
        }
    }
    let mut best = None;
    for value in maxima.into_iter().flatten() {
        if match best {
            Some(old) => value.less(old)?,
            None => true,
        } {
            best = Some(value);
        }
    }
    Ok(best.expect("a strongly connected component with an edge has a cycle"))
}

pub(super) fn mean_arc_ids(
    n: usize,
    arcs: &[RatioArc],
    negative_only: bool,
) -> Result<Option<Vec<usize>>, GraphError> {
    let mut b = Graph::builder(n, GraphKind::Directed);
    for a in arcs {
        b.add_edge(a.u, a.v, ())?;
    }
    let components = strongly_connected_components(&b.build());
    let mut size = vec![0; components.count];
    let mut local = vec![0; n];
    for (v, &c) in components.component_of.iter().enumerate() {
        local[v] = size[c];
        size[c] += 1;
    }
    let mut groups = vec![Vec::new(); size.len()];
    for (i, a) in arcs.iter().enumerate() {
        let c = components.component_of[a.u];
        if c == components.component_of[a.v] {
            groups[c].push(i);
        }
    }
    let mut best = None;
    let mut witness = None;
    for (c, ids) in groups.into_iter().enumerate() {
        if ids.is_empty() {
            continue;
        }
        let mapped: Vec<_> = ids
            .iter()
            .map(|&i| RatioArc {
                u: local[arcs[i].u],
                v: local[arcs[i].v],
                ..arcs[i]
            })
            .collect();
        let value = component_mean(size[c], &mapped)?;
        if negative_only && value.0 >= 0 {
            continue;
        }
        if match best {
            Some(old) => value.less(old)?,
            None => true,
        } {
            best = Some(value);
            witness = Some((size[c], mapped, ids));
        }
    }
    let (Some(value), Some((n, mapped, ids))) = (best, witness) else {
        return Ok(None);
    };
    let cycle = cycle_at_threshold(n, &mapped, value.0, value.1, true)?
        .expect("the minimum mean is attained");
    Ok(Some(cycle.into_iter().map(|i| ids[i]).collect()))
}

/// Exact minimum average edge cost of a directed cycle, with its original
/// edges in traversal order. Returns None on a DAG. The two-pass dynamic
/// program uses O(V+E) storage and O(V+E+sum_C |V_C| |E_C|) time over SCCs.
pub fn minimum_mean_cycle(graph: &Graph<i64>) -> Result<Option<RatioCycle>, GraphError> {
    graph.require(GraphKind::Directed)?;
    let arcs: Vec<_> = graph
        .edges()
        .iter()
        .map(|e| RatioArc {
            u: e.from,
            v: e.to,
            cost: i128::from(e.weight),
            length: 1,
        })
        .collect();
    Ok(
        mean_arc_ids(graph.node_count(), &arcs, false)?.map(|ids| RatioCycle {
            total_cost: ids.iter().map(|&i| arcs[i].cost).sum(),
            total_length: ids.len() as u128,
            steps: ids
                .into_iter()
                .map(|edge| CycleStep {
                    edge,
                    reversed: false,
                })
                .collect(),
        }),
    )
}
