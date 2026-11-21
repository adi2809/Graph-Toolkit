//! Dynamic metric repair on a deterministic forest-contraction hierarchy.
mod certify;
mod chain;
mod embedding;
mod paths;
mod solver;
use super::{CycleRatioWeight, CycleStep, RatioCycle};
use crate::{Graph, GraphError, GraphKind, NONE};
use chain::{Chain, Forest};
use embedding::Metric;
pub(super) use solver::finish_from_residual;
pub use solver::{CoreCycleFlow, max_flow_core_cycles};

/// Maintains signed cycle candidates for a fixed undirected multigraph.
/// Each level contracts disjoint forest stars and embeds parallel core edges
/// in a shortest parallel representative. A two-hop spanner bounds the local
/// load of replacement paths. The index retains closed walks from contractions
/// and both kinds of replacement. Batches repair dependent expressions.
/// The forest and spanners rebuild after a batch when at least
/// ceil(m/(1+floor(log2(m)))) distinct lengths have drifted beyond a
/// factor of two of the snapshot. Either operation can touch the hierarchy.
/// Use `certified_cycle` for a global ratio guarantee.
pub struct EmbeddedCycleHierarchy {
    chain: Chain,
    dual: certify::DualBound,
    drifted: usize,
    snapshot: Graph<CycleRatioWeight>,
}
impl EmbeddedCycleHierarchy {
    pub fn new(graph: &Graph<CycleRatioWeight>) -> Result<Self, GraphError> {
        Self::from_owned(graph.clone())
    }

    pub(in crate::flow) fn from_owned(graph: Graph<CycleRatioWeight>) -> Result<Self, GraphError> {
        graph.require(GraphKind::Undirected)?;
        if graph.edges().iter().any(|e| e.weight.length == 0) {
            return Err(GraphError::InvalidModelParameters);
        }

        let n = graph.node_count();
        let tree = crate::connectivity::low_stretch::build(
            n,
            graph
                .edges()
                .iter()
                .enumerate()
                .map(|(i, e)| (e.from, e.to, e.weight.length, i)),
        );
        let mut b = Graph::builder(n, GraphKind::Undirected);
        for id in tree {
            let e = &graph.edges()[id];
            b.add_edge(e.from, e.to, e.weight.cost)?;
        }
        let forest = b.build();
        let mut parent = vec![NONE; n];
        let mut depth = vec![0; n];
        let mut potential = vec![0i128; n];
        let mut stack = Vec::new();
        for root in 0..n {
            if parent[root] != NONE {
                continue;
            }
            parent[root] = root;
            stack.push(root);
            while let Some(u) = stack.pop() {
                for a in forest.out(u) {
                    if parent[a.node] == NONE {
                        parent[a.node] = u;
                        depth[a.node] = depth[u] + 1;
                        let e = &forest.edges()[a.edge];
                        let cost = i128::from(e.weight);
                        potential[a.node] = potential[u] + if e.from == u { cost } else { -cost };
                        stack.push(a.node);
                    }
                }
            }
        }
        let input = graph.edges().iter().enumerate().map(|(i, e)| {
            (
                i,
                e.from,
                e.to,
                Metric::arc(e.weight.cost, e.weight.length, true, true),
            )
        });
        Ok(Self {
            chain: Chain::build(Forest { parent, depth }, graph.edge_count(), input, None)?,
            dual: certify::DualBound::new(&graph, potential),
            drifted: 0,
            snapshot: graph,
        })
    }
    /// Replaces costs and positive lengths in original edge ID order supplied
    /// by the caller. Repeated IDs use their last update. Invalid batches are atomic.
    pub fn update(&mut self, changes: &[(usize, CycleRatioWeight)]) -> Result<(), GraphError> {
        for &(id, w) in changes {
            if id >= self.chain.leaves.len() {
                return Err(GraphError::InvalidEdge {
                    edge: id,
                    edge_count: self.chain.leaves.len(),
                });
            }
            if w.length == 0 {
                return Err(GraphError::InvalidModelParameters);
            }
        }
        // Coalesce duplicate IDs before deciding whether to pay for a
        // topology rebuild. One expensive rebuild thus accounts for at
        // least m/(1+log m) distinct edge changes since the previous one.
        let unique: std::collections::BTreeMap<_, _> = changes.iter().copied().collect();
        let mut drifted = self.drifted;
        for (&id, w) in &unique {
            let initial = self.snapshot.edges()[id].weight.length;
            let outside = |length: u64| {
                length > initial.saturating_mul(2) || initial > length.saturating_mul(2)
            };
            let current = self
                .chain
                .arena
                .metric(self.chain.leaves[id].unwrap())
                .length as u64;
            drifted -= usize::from(outside(current));
            drifted += usize::from(outside(w.length));
        }
        let m = self.snapshot.edge_count().max(1);
        if drifted >= m.div_ceil(1 + m.ilog2() as usize) {
            let mut builder = Graph::builder(self.snapshot.node_count(), GraphKind::Undirected);
            for (id, edge) in self.snapshot.edges().iter().enumerate() {
                let w = unique.get(&id).copied().unwrap_or_else(|| {
                    let metric = self.chain.arena.metric(self.chain.leaves[id].unwrap());
                    CycleRatioWeight {
                        cost: metric.gradient as i64,
                        length: metric.length as u64,
                    }
                });
                builder.add_edge(edge.from, edge.to, w)?;
            }
            *self = Self::from_owned(builder.build())?;
        } else {
            for (id, w) in unique {
                let edge = &self.snapshot.edges()[id];
                self.dual.update(id, edge.from, edge.to, w);
                self.chain.arena.update(
                    self.chain.leaves[id].unwrap(),
                    Metric::arc(w.cost, w.length, true, true),
                );
            }
            self.chain.arena.repair();
            self.drifted = drifted;
        }
        Ok(())
    }
    /// Returns an improving closed walk with signed original edge IDs. Walks
    /// may repeat edges. The denominator counts all traversals of the walk.
    pub fn improving_walk(&self) -> Option<RatioCycle> {
        let arena = &self.chain.arena;
        let best = arena.better_orientation(self.chain.winner?);
        let score = arena.metric(best);
        if score.gradient >= 0 {
            return None;
        }
        let mut walk = Vec::new();
        arena.lift(best, &mut walk);
        Some(RatioCycle {
            steps: walk
                .into_iter()
                .map(|(edge, reversed)| CycleStep { edge, reversed })
                .collect(),
            total_cost: score.gradient,
            total_length: score.length,
        })
    }

    /// Number of two-hop replacements across all levels of the current hierarchy.
    pub fn spanner_replacements(&self) -> usize {
        self.chain.shortcuts
    }

    /// Extracts the best simple cycle in the reported walk's loop-erasure
    /// decomposition. Its ratio is at least as favorable as the walk's.
    pub fn improving_cycle(&self) -> Option<RatioCycle> {
        let walk = self.improving_walk()?;
        let mut position = std::collections::BTreeMap::new();
        let mut path = Vec::<CycleStep>::new();
        let mut best: Option<(Metric, Vec<CycleStep>)> = None;
        let start_of = |s: CycleStep| {
            let e = &self.snapshot.edges()[s.edge];
            if s.reversed { e.to } else { e.from }
        };
        for step in walk.steps {
            let u = start_of(step);
            let v = start_of(CycleStep {
                reversed: !step.reversed,
                ..step
            });
            position.insert(u, path.len());
            path.push(step);
            if let Some(&index) = position.get(&v) {
                let cycle = path.split_off(index);
                let mut score = Metric {
                    gradient: 0,
                    length: 0,
                };
                for &a in &cycle {
                    position.remove(&start_of(a));
                    let leaf = self.chain.leaves[a.edge].unwrap();
                    let m = self
                        .chain
                        .arena
                        .metric(if a.reversed { leaf.reverse() } else { leaf });
                    score.gradient += m.gradient;
                    score.length += m.length;
                }
                if best
                    .as_ref()
                    .is_none_or(|(old, _)| score.compare_ratio(*old).is_lt())
                {
                    best = Some((score, cycle));
                }
            }
        }
        let (score, steps) = best?;
        Some(RatioCycle {
            steps,
            total_cost: score.gradient,
            total_length: score.length,
        })
    }
}
