use crate::NONE;

/// Indexed monotone radix heap for u128 priorities. Keys live in the caller's
/// distance array. Intrusive bucket lists make a decrease an O(1) relink.
pub(crate) struct RadixHeap {
    head: [usize; 129],
    links: Vec<[usize; 2]>,
    last: u128,
    occupied: [u64; 3],
}

impl RadixHeap {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            head: [NONE; 129],
            links: vec![[NONE; 2]; n],
            last: 0,
            occupied: [0; 3],
        }
    }

    pub(crate) fn reset(&mut self) {
        self.head.fill(NONE);
        self.last = 0;
        self.occupied.fill(0);
    }

    fn bucket(&self, key: u128) -> usize {
        (u128::BITS - (self.last ^ key).leading_zeros()) as usize
    }

    /// Inserts an undiscovered vertex (old == u128::MAX) or decreases a queued key.
    pub(crate) fn decrease(&mut self, node: usize, old: u128, key: u128) {
        debug_assert!(key >= self.last && key < old);
        if old != u128::MAX {
            self.unlink(self.bucket(old), node);
        }
        self.insert(self.bucket(key), node);
    }

    fn insert(&mut self, bucket: usize, node: usize) {
        let next = self.head[bucket];
        self.links[node] = [NONE, next];
        if next != NONE {
            self.links[next][0] = node;
        }
        self.head[bucket] = node;
        self.occupied[bucket / 64] |= 1 << (bucket % 64);
    }

    fn unlink(&mut self, bucket: usize, node: usize) {
        let [previous, next] = self.links[node];
        if previous == NONE {
            debug_assert_eq!(self.head[bucket], node);
            self.head[bucket] = next;
        } else {
            self.links[previous][1] = next;
        }
        if next != NONE {
            self.links[next][0] = previous;
        }
        if self.head[bucket] == NONE {
            self.occupied[bucket / 64] &= !(1 << (bucket % 64));
        }
    }

    /// Advances the monotone floor to the next key, which remains queued.
    pub(crate) fn peek_key(&mut self, keys: &[u128]) -> Option<u128> {
        if self.head[0] == NONE {
            let bucket = self
                .occupied
                .iter()
                .enumerate()
                .find_map(|(index, &mask)| {
                    (mask != 0).then(|| index * 64 + mask.trailing_zeros() as usize)
                })?;
            let mut node = self.head[bucket];
            let mut minimum = u128::MAX;
            while node != NONE {
                minimum = minimum.min(keys[node]);
                node = self.links[node][1];
            }
            self.last = minimum;

            // The new floor shares all higher bits with this bucket's keys.
            // Every redistributed entry moves down and higher buckets are unaffected.
            let mut node = std::mem::replace(&mut self.head[bucket], NONE);
            self.occupied[bucket / 64] &= !(1 << (bucket % 64));
            while node != NONE {
                let next = self.links[node][1];
                self.insert(self.bucket(keys[node]), node);
                node = next;
            }
        }
        Some(self.last)
    }

    pub(crate) fn pop(&mut self, keys: &[u128]) -> Option<(u128, usize)> {
        let key = self.peek_key(keys)?;
        let node = self.head[0];
        self.unlink(0, node);
        Some((key, node))
    }
}

/// Indexed 4-ary min-heap with a single live entry per vertex.
pub(crate) struct IndexedHeap<K> {
    keys: Vec<(K, usize)>,
    position: Vec<usize>,
}

impl<K: Copy + Ord> IndexedHeap<K> {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            keys: Vec::new(),
            position: vec![NONE; n],
        }
    }

    pub(crate) fn clear(&mut self) {
        for (_, node) in self.keys.drain(..) {
            self.position[node] = NONE;
        }
    }

    /// Inserts a vertex or decreases its key.
    pub(crate) fn decrease(&mut self, node: usize, key: K) {
        let mut index = self.position[node];
        let entry = (key, node);
        if index == NONE {
            index = self.keys.len();
            self.keys.push(entry);
            self.position[node] = index;
        } else {
            debug_assert!(self.keys[index] >= entry);
            self.keys[index] = entry;
        }
        while index > 0 && self.keys[(index - 1) / 4] > self.keys[index] {
            let parent = (index - 1) / 4;
            self.swap(index, parent);
            index = parent;
        }
    }

    pub(crate) fn pop(&mut self) -> Option<(K, usize)> {
        let entry = *self.keys.first()?;
        self.keys.swap_remove(0);
        self.position[entry.1] = NONE;
        if let Some(&(_, node)) = self.keys.first() {
            self.position[node] = 0;
        }
        self.sift_down(0);
        Some(entry)
    }

    pub(crate) fn peek(&self) -> Option<(K, usize)> {
        self.keys.first().copied()
    }

    /// Rebuilds from unique IDs in linear time.
    pub(crate) fn rebuild(&mut self, entries: Vec<(K, usize)>) {
        self.clear();
        self.keys = entries;
        for (position, &(_, id)) in self.keys.iter().enumerate() {
            self.position[id] = position;
        }
        for index in (0..self.keys.len().div_ceil(4)).rev() {
            self.sift_down(index);
        }
    }

    pub(crate) fn offer(&mut self, id: usize, key: K) {
        let pos = self.position[id];
        if pos == NONE || self.keys[pos].0 > key {
            self.decrease(id, key);
        }
    }

    fn sift_down(&mut self, mut index: usize) {
        while index * 4 + 1 < self.keys.len() {
            let start = index * 4 + 1;
            let end = (start + 4).min(self.keys.len());
            let child = (start..end).min_by_key(|&child| self.keys[child]).unwrap();
            if self.keys[index] <= self.keys[child] {
                break;
            }
            self.swap(index, child);
            index = child;
        }
    }

    fn swap(&mut self, a: usize, b: usize) {
        self.keys.swap(a, b);
        self.position[self.keys[a].1] = a;
        self.position[self.keys[b].1] = b;
    }
}

/// Indexed monotone priority queue with the full u128 key range. Insertions
/// and decreases must respect the last extracted minimum. Equal keys are allowed.
pub struct RadixQueue {
    heap: RadixHeap,
    keys: Vec<u128>,
    present: Vec<bool>,
    len: usize,
}
impl RadixQueue {
    pub fn new(capacity: usize) -> Self {
        Self {
            heap: RadixHeap::new(capacity),
            keys: vec![0; capacity],
            present: vec![false; capacity],
            len: 0,
        }
    }
    /// Returns false for an invalid ID, a key below the floor, or a
    /// non-decreasing update to an existing ID. Rejection preserves the queue.
    pub fn decrease(&mut self, id: usize, key: u128) -> bool {
        if id >= self.keys.len() || key < self.heap.last {
            return false;
        }
        if self.present[id] {
            if key >= self.keys[id] {
                return false;
            }
            self.heap.unlink(self.heap.bucket(self.keys[id]), id);
        } else {
            self.present[id] = true;
            self.len += 1;
        }
        self.heap.insert(self.heap.bucket(key), id);
        self.keys[id] = key;
        true
    }
    pub fn pop(&mut self) -> Option<(u128, usize)> {
        let result = self.heap.pop(&self.keys)?;
        self.present[result.1] = false;
        self.len -= 1;
        Some(result)
    }
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    pub fn clear(&mut self) {
        self.heap.reset();
        self.present.fill(false);
        self.len = 0;
    }
}
