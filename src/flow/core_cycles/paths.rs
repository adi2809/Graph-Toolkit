//! During a terminal epoch, each tree edge is used in one direction only.
//! Heavy paths and lazy segment minima apply a push without listing its
//! edges. We expose a changed leaf only when its length bucket or threshold
//! eligibility changes. At epoch end a linear pass materializes the flow.
use super::{Chain, embedding::Walk, solver::residual_metric};
use crate::{Graph, NONE};
use std::ops::Range;

#[derive(Clone, Copy)]
struct Minimum {
    room: u128,
    until_change: u128,
    pending: u128,
}
const EMPTY: Minimum = Minimum {
    room: u128::MAX,
    until_change: u128::MAX,
    pending: 0,
};

struct Ranges {
    data: Vec<Minimum>,
    base: usize,
    threshold: u128,
}

fn budget(room: u128, threshold: u128) -> u128 {
    if room < threshold {
        u128::MAX
    } else {
        room - threshold.max(1u128 << room.ilog2()) + 1
    }
}

impl Ranges {
    fn new(rooms: &[u128], threshold: u128) -> Self {
        let base = rooms.len().max(1).next_power_of_two();
        let mut tree = Self {
            data: vec![EMPTY; 2 * base],
            base,
            threshold,
        };
        for (leaf, &room) in tree.data[base..].iter_mut().zip(rooms) {
            *leaf = Minimum {
                room,
                until_change: budget(room, threshold),
                pending: 0,
            };
        }
        for i in (1..base).rev() {
            tree.pull(i);
        }
        tree
    }
    fn pull(&mut self, i: usize) {
        self.data[i].room = self.data[2 * i].room.min(self.data[2 * i + 1].room);
        self.data[i].until_change = self.data[2 * i]
            .until_change
            .min(self.data[2 * i + 1].until_change);
    }
    fn apply(&mut self, i: usize, amount: u128) {
        let node = &mut self.data[i];
        node.room -= amount;
        node.until_change = node.until_change.saturating_sub(amount);
        if i < self.base {
            node.pending += amount;
        }
    }
    fn propagate(&mut self, i: usize) {
        let pending = std::mem::take(&mut self.data[i].pending);
        if pending != 0 {
            self.apply(2 * i, pending);
            self.apply(2 * i + 1, pending);
        }
    }
    fn bottleneck(&mut self, i: usize, range: Range<usize>, q: &Range<usize>) -> u128 {
        if q.start <= range.start && range.end <= q.end {
            return self.data[i].room;
        }
        self.propagate(i);
        let mid = range.start + (range.end - range.start) / 2;
        let mut room = u128::MAX;
        if q.start < mid {
            room = self.bottleneck(2 * i, range.start..mid, q);
        }
        if q.end > mid {
            room = room.min(self.bottleneck(2 * i + 1, mid..range.end, q));
        }
        room
    }
    fn subtract(&mut self, i: usize, range: Range<usize>, q: &Range<usize>, amount: u128) {
        if q.start <= range.start && range.end <= q.end {
            self.apply(i, amount);
            return;
        }
        self.propagate(i);
        let mid = range.start + (range.end - range.start) / 2;
        if q.start < mid {
            self.subtract(2 * i, range.start..mid, q, amount);
        }
        if q.end > mid {
            self.subtract(2 * i + 1, mid..range.end, q, amount);
        }
        self.pull(i);
    }
    fn detect(&mut self, i: usize, emit: &mut impl FnMut(usize, u128)) {
        if self.data[i].until_change != 0 {
            return;
        }
        if i >= self.base {
            let room = self.data[i].room;
            emit(i - self.base, room);
            self.data[i].until_change = budget(room, self.threshold);
        } else {
            self.propagate(i);
            self.detect(i * 2, emit);
            self.detect(i * 2 + 1, emit);
            self.pull(i);
        }
    }
}

pub(super) struct PathFlow {
    head: Vec<usize>,
    pos: Vec<usize>,
    parent: Vec<usize>,
    // Indexed by heavy-path position, excluding both roots.
    edges: Vec<(usize, usize)>,
    ranges: Ranges,
    marked: Vec<bool>,
    pieces: Vec<Range<usize>>,
}

impl PathFlow {
    pub(super) fn new(
        g: &Graph<u64>,
        parent: Vec<usize>,
        incoming: &[usize],
        level: &[usize],
        cap: &[[u128; 2]],
        threshold: u128,
        chain: &Chain,
    ) -> Self {
        let n = parent.len();
        let mut first = vec![NONE; n];
        let mut next = vec![NONE; n];
        let mut order = Vec::new();
        for (v, &p) in parent.iter().enumerate() {
            if p == v {
                order.push(v);
            } else if p != NONE {
                next[v] = first[p];
                first[p] = v;
            }
        }
        let mut i = 0;
        while i < order.len() {
            let mut child = first[order[i]];
            while child != NONE {
                order.push(child);
                child = next[child];
            }
            i += 1;
        }
        let mut size = vec![1; n];
        let mut heavy = vec![NONE; n];
        for &v in order.iter().rev() {
            let p = parent[v];
            if p != v {
                size[p] += size[v];
                if heavy[p] == NONE || size[v] > size[heavy[p]] {
                    heavy[p] = v;
                }
            }
        }
        let mut head = vec![NONE; n];
        let mut pos = vec![NONE; n];
        let mut edges = Vec::new();
        let mut rooms = Vec::new();
        for &v in &order {
            let p = parent[v];
            if p == v || (heavy[p] == v && parent[p] != p) {
                continue;
            }
            let mut child = v;
            while child != NONE {
                head[child] = v;
                pos[child] = edges.len();
                let id = incoming[child];
                let e = &g.edges()[id];
                let dir = usize::from(level[e.from] > level[e.to]);
                edges.push((id, dir));
                rooms.push(cap[id][dir]);
                child = heavy[child];
            }
        }
        let marked = chain.arena.mark_leaves(|id| {
            let e = &g.edges()[id];
            incoming[e.from] != id && incoming[e.to] != id
        });
        Self {
            head,
            pos,
            parent,
            edges,
            ranges: Ranges::new(&rooms, threshold),
            marked,
            pieces: Vec::new(),
        }
    }
    fn collect(&mut self, mut v: usize) {
        while self.parent[v] != v {
            let h = self.head[v];
            self.pieces.push(self.pos[h]..self.pos[v] + 1);
            v = self.parent[h];
        }
    }
    pub(super) fn push(
        &mut self,
        g: &Graph<u64>,
        chain: &mut Chain,
        best: Walk,
        cap: &mut [[u128; 2]],
    ) -> u128 {
        let (id, back) = chain.arena.marked_leaf(best, &self.marked);
        let dir = usize::from(back);
        let e = &g.edges()[id];
        self.pieces.clear();
        self.collect(e.from);
        self.collect(e.to);
        let mut amount = cap[id][dir];
        let base = self.ranges.base;
        for q in &self.pieces {
            amount = amount.min(self.ranges.bottleneck(1, 0..base, q));
        }
        cap[id][dir] -= amount;
        cap[id][dir ^ 1] += amount;
        let threshold = self.ranges.threshold;
        let mut notify = |id: usize, pair: [u128; 2], dir: usize| {
            let mut metric = residual_metric(pair, threshold);
            metric[dir ^ 1] = super::embedding::CLOSED;
            chain.arena.update(chain.leaves[id].unwrap(), metric);
        };
        notify(id, cap[id], dir);
        for q in &self.pieces {
            self.ranges.subtract(1, 0..base, q, amount);
        }
        self.ranges.detect(1, &mut |pos, room| {
            let (id, dir) = self.edges[pos];
            let mut pair = cap[id];
            pair[dir ^ 1] = pair[0] + pair[1] - room;
            pair[dir] = room;
            notify(id, pair, dir);
        });
        amount
    }
    pub(super) fn flush(mut self, cap: &mut [[u128; 2]]) {
        for i in 1..self.ranges.base {
            self.ranges.propagate(i);
        }
        for ((id, dir), leaf) in self
            .edges
            .into_iter()
            .zip(&self.ranges.data[self.ranges.base..])
        {
            let total = cap[id][0] + cap[id][1];
            cap[id][dir] = leaf.room;
            cap[id][dir ^ 1] = total - leaf.room;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn range_pushes_report_exactly_the_changed_buckets() {
        let mut seed = 335532u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut rooms: Vec<_> = (0..137).map(|_| u128::from(next()) << 32).collect();
        let threshold = 1u128 << 70;
        let mut index = Ranges::new(&rooms, threshold);
        let base = index.base;
        let class = |r: u128| (r >= threshold).then(|| r.ilog2());
        for _ in 0..200 {
            let before = rooms.clone();
            for _ in 0..3 {
                let start = next() as usize % rooms.len();
                let end = start + 1 + next() as usize % (rooms.len() - start);
                let min = *rooms[start..end].iter().min().unwrap();
                assert_eq!(index.bottleneck(1, 0..base, &(start..end)), min);
                let delta = (u128::from(next()) << 32) % (min + 1);
                index.subtract(1, 0..base, &(start..end), delta);
                rooms[start..end].iter_mut().for_each(|r| *r -= delta);
            }
            let changed: Vec<_> = before
                .iter()
                .zip(&rooms)
                .enumerate()
                .filter_map(|(i, (&a, &b))| (class(a) != class(b)).then_some((i, b)))
                .collect();
            let mut reported = vec![];
            index.detect(1, &mut |i, room| reported.push((i, room)));
            assert_eq!(changed, reported);
        }
    }
}
