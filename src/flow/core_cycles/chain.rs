//! Contract alternate depths of a rooted forest. Each lifted core edge
//! joins two ports through at most three walks from the preceding level.
use super::embedding::{Embedding, Metric, Walk};
use crate::{GraphError, NONE};

struct CoreEdge {
    u: usize,
    v: usize,
    walk: Walk,
}
impl CoreEdge {
    fn canonical(mut self) -> Self {
        if self.u > self.v {
            std::mem::swap(&mut self.u, &mut self.v);
            self.walk = self.walk.reverse();
        }
        self
    }
}
#[path = "spanner.rs"]
mod spanner;
struct Builder {
    arena: Embedding,
    cycles: Vec<Walk>,
    record_cycles: bool,
    shortcuts: usize,
}
impl Builder {
    fn close(&mut self, w: Walk) {
        if self.record_cycles {
            self.cycles.push(w);
        }
    }
    fn bundle(&mut self, walks: &[CoreEdge]) -> Walk {
        if walks.len() == 1 {
            return walks[0].walk;
        }
        let mid = walks.len() / 2;
        let a = self.bundle(&walks[..mid]);
        let b = self.bundle(&walks[mid..]);
        if self.record_cycles {
            let cycle = self.arena.join(a, b.reverse());
            self.close(cycle);
        }
        self.arena.parallel(a, b)
    }
    fn sparsify(&mut self, mut edges: Vec<CoreEdge>, forest: &Forest) -> Vec<CoreEdge> {
        // Stable ordering also fixes tie-breaking between parallel witnesses.
        edges.sort_by_key(|e| (e.u, e.v));
        let mut spanner = Vec::new();
        for group in edges.chunk_by(|a, b| a.u == b.u && a.v == b.v) {
            let first = &group[0];
            if first.u == first.v {
                for edge in group {
                    self.close(edge.walk);
                }
            } else {
                let walk = self.bundle(group);
                spanner.push(CoreEdge {
                    u: first.u,
                    v: first.v,
                    walk,
                });
            }
        }
        // The two terminal trees are joined only at a single pair of
        // depths. Each contraction preserves that property, so each
        // terminal core is bipartite and has no two-hop triangle shortcuts.
        if self.record_cycles {
            spanner::two_hop(self, spanner, forest)
        } else {
            spanner
        }
    }
}

pub(super) struct Chain {
    pub(super) shortcuts: usize,
    pub(super) arena: Embedding,
    pub(super) leaves: Vec<Option<Walk>>,
    pub(super) winner: Option<Walk>,
}
pub(super) struct Forest {
    pub(super) parent: Vec<usize>,
    pub(super) depth: Vec<usize>,
}
impl Chain {
    // In terminal mode the only negative cost is the closing edge t -> s.
    // Every cycle internal to a core has zero objective. Suppressing its
    // score saves space without changing the optimum among represented cycles.
    pub(super) fn build(
        mut forest: Forest,
        edge_count: usize,
        input: impl Iterator<Item = (usize, usize, usize, [Metric; 2])>,
        terminals: Option<(usize, usize)>,
    ) -> Result<Self, GraphError> {
        let mut b = Builder {
            arena: Embedding::default(),
            cycles: Vec::new(),
            shortcuts: 0,
            record_cycles: terminals.is_none(),
        };
        let n = forest.parent.len();
        let mut active: Vec<_> = (0..n).filter(|&u| forest.parent[u] != NONE).collect();
        let height = forest.depth.iter().copied().max().unwrap_or(0);
        // A lifted edge expands by at most three at each contraction.
        // Check a worst-case bound for all subsequent i64/u64 metric edits.
        if (u128::from(u64::MAX) * 3)
            .checked_mul(
                3u128
                    .checked_pow(usize::BITS - height.leading_zeros())
                    .ok_or(GraphError::GraphSizeOverflow)?,
            )
            .is_none_or(|v| v > i128::MAX as u128)
        {
            return Err(GraphError::GraphSizeOverflow);
        }
        let mut leaves = vec![None; edge_count];
        let mut initial = Vec::new();
        for (id, u, v, metric) in input {
            let walk = b.arena.edge(id, metric);
            leaves[id] = Some(walk);
            initial.push(CoreEdge { u, v, walk }.canonical());
        }
        let mut core = b.sparsify(initial, &forest);
        let mut height = height;
        let mut owner: Vec<_> = (0..n).collect();
        let mut port = vec![None; n];
        while height != 0 {
            for &u in &active {
                owner[u] = if forest.depth[u] % 2 == 1 {
                    forest.parent[u]
                } else {
                    u
                };
            }
            for e in &core {
                if owner[e.v] == e.u {
                    port[e.v] = Some(e.walk);
                } else if owner[e.u] == e.v {
                    port[e.u] = Some(e.walk.reverse());
                }
            }
            let mut lifted = Vec::with_capacity(core.len());
            for e in core {
                // Contract the chosen forest edge. Parallel alternatives
                // already have their closed walks recorded in the bundle.
                if owner[e.v] == e.u || owner[e.u] == e.v {
                    continue;
                }
                let (u, v) = (owner[e.u], owner[e.v]);
                if u == v && !b.record_cycles {
                    continue;
                }
                let mut walk = e.walk;
                if let Some(p) = port[e.u] {
                    walk = b.arena.join(p, walk);
                }
                if let Some(p) = port[e.v] {
                    walk = b.arena.join(walk, p.reverse());
                }
                lifted.push(CoreEdge { u, v, walk }.canonical());
            }

            active.retain(|&u| u == owner[u]);
            for &u in &active {
                forest.parent[u] = owner[forest.parent[u]];
                forest.depth[u] /= 2;
            }
            core = b.sparsify(lifted, &forest);
            height /= 2;
        }
        let winner = if let Some((s, t)) = terminals {
            core.into_iter()
                .find(|e| (e.u, e.v) == (s.min(t), s.max(t)))
                .map(|e| if e.u == s { e.walk } else { e.walk.reverse() })
        } else {
            b.arena.contest(b.cycles)
        };
        b.arena.seal();
        Ok(Self {
            shortcuts: b.shortcuts,
            arena: b.arena,
            leaves,
            winner,
        })
    }
}
