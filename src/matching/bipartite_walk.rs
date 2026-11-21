//! Alternating-path algorithms on the bipartite CSR representation.
use super::{BipartiteGraph, BipartiteMatching, NONE};
use crate::random::Random;

impl BipartiteGraph {
    /// Kuhn's one-root-at-a-time augmenting path search, O(L(E+1)+V).
    pub fn maximum_matching_kuhn(&self) -> BipartiteMatching {
        let mut m = self.empty_matching();
        let mut seen = vec![NONE; self.left_count()];
        let mut path = Vec::new();
        for root in 0..self.left_count() {
            path.clear();
            path.push((root, 0, NONE));
            seen[root] = root;
            while let Some((left, pos, _)) = path.last_mut() {
                let Some(&v) = self.neighbors(*left).get(*pos) else {
                    path.pop();
                    continue;
                };
                *pos += 1;
                let next = m.right[v];
                if next == NONE {
                    let mut end = v;
                    for (u, _, incoming) in path.drain(..).rev() {
                        m.pair(u, end);
                        end = incoming;
                    }
                    m.len += 1;
                    break;
                }
                if seen[next] != root {
                    seen[next] = root;
                    path.push((next, 0, v));
                }
            }
        }
        self.finish_matching(m, &(0..self.left_count()).collect::<Vec<_>>())
    }

    /// Goel–Kapralov–Khanna truncated, loop-erased random walks. Exact on
    /// arbitrary bipartite multigraphs. Failed attempts have a finite budget
    /// before Hopcroft–Karp finishes and certifies the remaining matching.
    pub fn random_walk_maximum_matching(&self, seed: u64) -> BipartiteMatching {
        let mut rng = Random(seed);
        let n = self.left_count();
        let mut free: Vec<_> = (0..n).filter(|&v| !self.neighbors(v).is_empty()).collect();
        let mut m = self.empty_matching();
        let mut position = vec![NONE; n];
        let mut trail = Vec::new();
        let mut misses = 0;
        let budget = 4 + 2 * (usize::BITS - n.leading_zeros());
        while !free.is_empty() && m.len < self.right_count && misses < budget {
            let slot = rng.below(free.len());
            let mut at = free[slot];
            let limit = n.div_ceil(free.len()).saturating_add(1).saturating_mul(2);
            let mut success = false;
            for _ in 0..limit {
                let neighbors = self.neighbors(at);
                let v = neighbors[rng.below(neighbors.len())];
                position[at] = trail.len();
                trail.push((at, v));
                at = m.right[v];
                if at == NONE {
                    success = true;
                    break;
                }
                // Re-entering a contracted matched pair closes a cycle.
                if position[at] != NONE {
                    for (u, _) in trail.drain(position[at]..) {
                        position[u] = NONE;
                    }
                }
            }
            for (u, v) in trail.drain(..) {
                position[u] = NONE;
                if success {
                    m.pair(u, v);
                }
            }
            if success {
                free.swap_remove(slot);
                m.len += 1;
                misses = 0;
            } else {
                misses += 1;
            }
        }
        if m.len == n {
            m.cover.left = (0..n).collect();
            m
        } else if m.len == self.right_count {
            m.cover.right = (0..self.right_count).collect();
            m
        } else {
            self.finish_matching(m, &(0..n).collect::<Vec<_>>())
        }
    }
}
