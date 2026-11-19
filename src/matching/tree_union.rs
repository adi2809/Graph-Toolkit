//! Gabow–Tarjan incremental tree set union. An ancestor mask replaces the
//! paper's microset lookup table on a word RAM with a highest-set-bit operation.

use crate::NONE;
use crate::connectivity::DisjointSet;

#[derive(Default)]
struct Microset {
    boundary: usize,
    vertices: Vec<usize>,
    marked: usize,
}

/// Links deactivate a vertex's mark. A find returns its nearest marked ancestor.
pub(super) struct TreeUnion {
    block_size: usize,
    microsets: Vec<Microset>,
    parent: Vec<usize>,
    block: Vec<usize>,
    slot: Vec<usize>,
    ancestors: Vec<usize>,
    macrosets: DisjointSet,
    macro_top: Vec<usize>,
}

impl TreeUnion {
    pub(super) fn new(vertices: usize) -> Self {
        // A permanent virtual root lets a forest share its top microset.
        let mut forest = Self {
            block_size: ((vertices + 1).ilog2() as usize).max(8),
            microsets: vec![Microset {
                boundary: NONE,
                ..Microset::default()
            }],
            parent: vec![NONE; vertices + 1],
            block: vec![NONE; vertices + 1],
            slot: vec![0; vertices + 1],
            ancestors: vec![0; vertices + 1],
            macrosets: DisjointSet::new(vertices + 1),
            macro_top: (0..=vertices).collect(),
        };
        forest.insert(0, vertices, true);
        forest
    }

    pub(super) fn add_root(&mut self, vertex: usize) {
        self.grow(self.parent.len() - 1, vertex);
    }

    /// Adds a previously absent leaf below an existing vertex, including a linked one.
    pub(super) fn grow(&mut self, parent: usize, vertex: usize) {
        debug_assert_eq!(self.block[vertex], NONE);
        self.parent[vertex] = parent;
        let block = self.block[parent];
        self.insert(block, vertex, true);
        if self.microsets[block].vertices.len() == self.block_size {
            self.split(block);
        }
    }

    pub(super) fn link(&mut self, vertex: usize) {
        let mask = 1 << self.slot[vertex];
        let marks = &mut self.microsets[self.block[vertex]].marked;
        debug_assert_ne!(*marks & mask, 0);
        *marks &= !mask;
    }

    pub(super) fn find(&mut self, vertex: usize) -> usize {
        let mut top = self.micro_find(vertex);
        if self.block[vertex] == self.block[top] {
            return top;
        }
        let mut macro_root = self.macrosets.find(top);
        top = self.macro_top[macro_root];
        loop {
            let next = self.micro_find(top);
            if self.block[top] == self.block[next] {
                return next;
            }
            let next_root = self.macrosets.find(next);
            top = self.macro_top[next_root];
            macro_root = self.macrosets.union(macro_root, next_root).unwrap();
            self.macro_top[macro_root] = top;
        }
    }

    fn micro_find(&self, vertex: usize) -> usize {
        let microset = &self.microsets[self.block[vertex]];
        let marked_ancestors = microset.marked & self.ancestors[vertex];
        if marked_ancestors == 0 {
            microset.boundary
        } else {
            microset.vertices[marked_ancestors.ilog2() as usize]
        }
    }

    fn insert(&mut self, block: usize, vertex: usize, marked: bool) {
        let parent = self.parent[vertex];
        let inherited = if parent == self.microsets[block].boundary {
            0
        } else {
            debug_assert_eq!(self.block[parent], block);
            self.ancestors[parent]
        };
        let microset = &mut self.microsets[block];
        let slot = microset.vertices.len();
        let bit = 1 << slot;
        self.block[vertex] = block;
        self.slot[vertex] = slot;
        self.ancestors[vertex] = inherited | bit;
        microset.vertices.push(vertex);
        if marked {
            microset.marked |= bit;
        }
    }

    fn split(&mut self, block: usize) {
        const WORD: usize = usize::BITS as usize;
        let old = std::mem::take(&mut self.microsets[block]);
        self.microsets[block].boundary = old.boundary;
        let len = old.vertices.len();
        let threshold = self.block_size / 4;
        // Slots are in parent-before-child order. Each pending chain contains
        // complete residual child subtrees, so a split has a single boundary.
        let mut next = [NONE; WORD];
        let mut pending = [Chain::EMPTY; WORD + 1];
        for slot in (0..len).rev() {
            let vertex = old.vertices[slot];
            let parent = self.parent[vertex];
            let parent_slot = if parent == old.boundary {
                len
            } else {
                self.slot[parent]
            };
            let chain = pending[slot].prepend(slot, &mut next);
            pending[parent_slot].append(chain, &mut next);
            if pending[parent_slot].len >= threshold {
                let chain = std::mem::replace(&mut pending[parent_slot], Chain::EMPTY);
                let id = self.microsets.len();
                self.microsets.push(Microset {
                    boundary: parent,
                    ..Microset::default()
                });
                self.rebuild(id, chain, &old, &next);
            }
        }
        self.rebuild(block, pending[len], &old, &next);
    }

    fn rebuild(&mut self, block: usize, chain: Chain, old: &Microset, next: &[usize]) {
        let mut slot = chain.first;
        while slot != NONE {
            self.insert(block, old.vertices[slot], old.marked & (1 << slot) != 0);
            slot = next[slot];
        }
    }
}

#[derive(Clone, Copy)]
struct Chain {
    first: usize,
    last: usize,
    len: usize,
}

impl Chain {
    const EMPTY: Self = Self {
        first: NONE,
        last: NONE,
        len: 0,
    };

    fn prepend(self, node: usize, next: &mut [usize]) -> Self {
        next[node] = self.first;
        Self {
            first: node,
            last: if self.len == 0 { node } else { self.last },
            len: self.len + 1,
        }
    }

    fn append(&mut self, other: Self, next: &mut [usize]) {
        if self.len == 0 {
            *self = other;
        } else {
            next[self.last] = other.first;
            self.last = other.last;
            self.len += other.len;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn incremental_links_agree_with_parent_walks_across_splits() {
        let n = 5_000;
        for chain in [false, true] {
            let mut forest = TreeUnion::new(n);
            let mut parents = vec![NONE; n];
            let mut marked = vec![true; n];
            let mut seed = 239_u64;
            let mut random = |upper: usize| {
                seed ^= seed << 13;
                seed ^= seed >> 7;
                seed ^= seed << 17;
                seed as usize % upper
            };
            for vertex in 0..n {
                if vertex == 0 || (!chain && random(17) == 0) {
                    forest.add_root(vertex);
                } else {
                    let parent = if chain { vertex - 1 } else { random(vertex) };
                    forest.grow(parent, vertex);
                    parents[vertex] = parent;
                }
                for _ in 0..4 {
                    let candidate = random(vertex + 1);
                    if parents[candidate] != NONE && marked[candidate] {
                        forest.link(candidate);
                        marked[candidate] = false;
                    }
                    let query = random(vertex + 1);
                    let mut expected = query;
                    while !marked[expected] {
                        expected = parents[expected];
                    }
                    assert_eq!(forest.find(query), expected);
                }
            }
        }
    }
}
