//! A vertex potential cancels on every closed walk. The worst absolute
//! reduced gradient per unit length therefore bounds every cycle's ratio.
use super::{
    CycleRatioWeight, CycleStep, EmbeddedCycleHierarchy, Graph, GraphError, Metric, RatioCycle,
};
use crate::{
    flow::ratio::{RatioArc, cycle_at_threshold},
    structures::SegmentTree,
};

type MinimumTree = SegmentTree<Metric, fn(&Metric, &Metric) -> Metric>;
const ZERO: Metric = Metric {
    gradient: 0,
    length: 1,
};

pub(super) struct DualBound {
    potential: Vec<i128>,
    edges: MinimumTree,
}
impl DualBound {
    fn edge(p: &[i128], u: usize, v: usize, w: CycleRatioWeight) -> Metric {
        Metric {
            gradient: -(i128::from(w.cost) + p[u] - p[v]).abs(),
            length: u128::from(w.length),
        }
    }

    pub(super) fn new(graph: &Graph<CycleRatioWeight>, potential: Vec<i128>) -> Self {
        let scores: Vec<_> = graph
            .edges()
            .iter()
            .map(|e| Self::edge(&potential, e.from, e.to, e.weight))
            .collect();
        Self {
            potential,
            edges: SegmentTree::new(&scores, ZERO, |a, b| {
                if a.compare_ratio(*b).is_lt() { *a } else { *b }
            }),
        }
    }
    pub(super) fn update(&mut self, id: usize, u: usize, v: usize, w: CycleRatioWeight) {
        self.edges.set(id, Self::edge(&self.potential, u, v, w));
    }
}

fn threshold(c: &RatioCycle, factor: u32) -> Result<Metric, GraphError> {
    // Cancel before scaling to preserve precision on wide integer inputs.
    let mut a = c.total_cost.unsigned_abs();
    let mut b = c.total_length;
    while b != 0 {
        (a, b) = (b, a % b);
    }
    Ok(Metric {
        gradient: (c.total_cost / a as i128)
            .checked_mul(i128::from(factor))
            .ok_or(GraphError::WeightOverflow)?,
        length: c.total_length / a,
    })
}

impl EmbeddedCycleHierarchy {
    /// Globally certified minimum-ratio approximation on the current signed
    /// graph. If the minimum ratio is r < 0, the witness has ratio <= r/factor.
    /// A factor of one gives an exact minimum. None certifies no negative cycle.
    ///
    /// The maintained dual bound can certify a lifted cycle immediately.
    /// Otherwise each Bellman–Ford threshold test costs O(VE) time and O(V+E)
    /// space. Every additional witness strictly improves the threshold.
    /// Returns WeightOverflow if a threshold test exceeds signed 128-bit arithmetic.
    pub fn certified_cycle(&self, factor: u32) -> Result<Option<RatioCycle>, GraphError> {
        if factor == 0 {
            return Err(GraphError::InvalidModelParameters);
        }
        let bound = self
            .dual
            .edges
            .query(0..self.snapshot.edge_count())
            .unwrap();
        if bound.gradient == 0 {
            return Ok(None);
        }
        let mut candidate = self.improving_cycle();
        let mut target = candidate
            .as_ref()
            .map(|c| threshold(c, factor))
            .transpose()?
            .unwrap_or(ZERO);
        if !target.compare_ratio(bound).is_gt() {
            return Ok(candidate);
        }
        let mut arcs = Vec::with_capacity(self.snapshot.edge_count() * 2);
        for (id, e) in self.snapshot.edges().iter().enumerate() {
            let m = self.chain.arena.metric(self.chain.leaves[id].unwrap());
            for (u, v, cost) in [(e.from, e.to, m.gradient), (e.to, e.from, -m.gradient)] {
                arcs.push(RatioArc {
                    u,
                    v,
                    cost,
                    length: m.length as u64,
                });
            }
        }
        while target.compare_ratio(bound).is_gt() {
            let denominator =
                i128::try_from(target.length).map_err(|_| GraphError::WeightOverflow)?;
            let Some(ids) = cycle_at_threshold(
                self.snapshot.node_count(),
                &arcs,
                target.gradient,
                denominator,
                false,
            )?
            else {
                break;
            };
            let better = RatioCycle {
                total_cost: ids.iter().map(|&i| arcs[i].cost).sum(),
                total_length: ids.iter().map(|&i| u128::from(arcs[i].length)).sum(),
                steps: ids
                    .into_iter()
                    .map(|i| CycleStep {
                        edge: i / 2,
                        reversed: i % 2 == 1,
                    })
                    .collect(),
            };
            target = threshold(&better, factor)?;
            candidate = Some(better);
        }
        Ok(candidate)
    }
}
