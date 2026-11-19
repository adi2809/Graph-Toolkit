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
        let mut index = 0;
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
        Some(entry)
    }

    fn swap(&mut self, a: usize, b: usize) {
        self.keys.swap(a, b);
        self.position[self.keys[a].1] = a;
        self.position[self.keys[b].1] = b;
    }
}
