//! Data structures for dense vertex identifiers and range queries.
pub use crate::connectivity::{DisjointSet as UnionFind, RollbackUnionFind};
pub use crate::queue::RadixQueue;
mod link_cut;
pub use link_cut::LinkCutForest;
mod fibonacci;
pub use fibonacci::FibonacciHeap;
use std::ops::Range;

/// Rooted ordered tree of bounded arity. Node IDs remain stable as it grows.
pub struct KaryTree<T> {
    arity: usize,
    nodes: Vec<TreeNode<T>>,
}
struct TreeNode<T> {
    value: T,
    parent: Option<usize>,
    children: Vec<usize>,
}

impl<T> KaryTree<T> {
    pub fn new(arity: usize, root: T) -> Self {
        Self {
            arity,
            nodes: vec![TreeNode {
                value: root,
                parent: None,
                children: vec![],
            }],
        }
    }
    /// Returns the new ID, or the supplied value if the parent is invalid or full.
    pub fn add_child(&mut self, parent: usize, value: T) -> Result<usize, T> {
        let id = self.nodes.len();
        let Some(p) = self.nodes.get_mut(parent) else {
            return Err(value);
        };
        if p.children.len() == self.arity {
            return Err(value);
        }
        p.children.push(id);
        self.nodes.push(TreeNode {
            value,
            parent: Some(parent),
            children: vec![],
        });
        Ok(id)
    }
    pub fn value(&self, id: usize) -> Option<&T> {
        self.nodes.get(id).map(|n| &n.value)
    }
    pub fn parent(&self, id: usize) -> Option<usize> {
        self.nodes.get(id)?.parent
    }
    pub fn children(&self, id: usize) -> Option<&[usize]> {
        Some(&self.nodes.get(id)?.children)
    }
}

/// Ordered range folds with an associative operation and a two-sided identity.
/// Construction is linear. Point replacement and half-open queries cost O(log n).
pub struct SegmentTree<T, F> {
    data: Vec<T>,
    n: usize,
    identity: T,
    combine: F,
}
impl<T: Clone, F: Fn(&T, &T) -> T> SegmentTree<T, F> {
    pub fn new(values: &[T], identity: T, combine: F) -> Self {
        let n = values.len();
        let mut data = vec![identity.clone(); 2 * n];
        data[n..].clone_from_slice(values);
        for i in (1..n).rev() {
            data[i] = combine(&data[i * 2], &data[i * 2 + 1]);
        }
        Self {
            data,
            n,
            identity,
            combine,
        }
    }
    pub fn set(&mut self, index: usize, value: T) -> bool {
        if index >= self.n {
            return false;
        }
        let mut p = self.n + index;
        self.data[p] = value;
        while p > 1 {
            p /= 2;
            self.data[p] = (self.combine)(&self.data[p * 2], &self.data[p * 2 + 1]);
        }
        true
    }
    pub fn query(&self, range: Range<usize>) -> Option<T> {
        if range.start > range.end || range.end > self.n {
            return None;
        }
        let (mut l, mut r) = (range.start + self.n, range.end + self.n);
        let (mut prefix, mut suffix) = (self.identity.clone(), self.identity.clone());
        while l < r {
            if l % 2 == 1 {
                prefix = (self.combine)(&prefix, &self.data[l]);
                l += 1;
            }
            if r % 2 == 1 {
                r -= 1;
                suffix = (self.combine)(&self.data[r], &suffix);
            }
            l /= 2;
            r /= 2;
        }
        Some((self.combine)(&prefix, &suffix))
    }
}

/// Fenwick's binary indexed tree for nonnegative frequencies. Prefix sums,
/// point replacement and inverse prefix search each take O(log n) time.
pub struct FenwickTree {
    buckets: Vec<usize>,
    total: usize,
}
impl FenwickTree {
    pub fn new(n: usize) -> Self {
        Self {
            buckets: vec![0; n + 1],
            total: 0,
        }
    }
    pub fn total(&self) -> usize {
        self.total
    }
    /// Sum of entries in 0..end.
    pub fn prefix_sum(&self, mut end: usize) -> Option<usize> {
        if end >= self.buckets.len() {
            return None;
        }
        let mut sum = 0;
        while end != 0 {
            sum += self.buckets[end];
            end &= end - 1;
        }
        Some(sum)
    }
    /// Sets a frequency and returns its old value. Invalid indices and total
    /// overflow return None before any buckets change.
    pub fn replace(&mut self, index: usize, weight: usize) -> Option<usize> {
        let mut i = index.checked_add(1)?;
        let old = self.prefix_sum(i)? - self.prefix_sum(index)?;
        self.total = (self.total - old).checked_add(weight)?;
        while i < self.buckets.len() {
            self.buckets[i] = (self.buckets[i] - old) + weight;
            i += i & i.wrapping_neg();
        }
        Some(old)
    }
    /// Index containing the zero-based rank in the expanded frequency table.
    pub fn select(&self, mut rank: usize) -> Option<usize> {
        if rank >= self.total {
            return None;
        }
        let mut bit = 1 << (self.buckets.len() - 1).ilog2();
        let mut before = 0;
        while bit != 0 {
            let candidate = before | bit;
            if let Some(&count) = self.buckets.get(candidate) {
                if count <= rank {
                    before = candidate;
                    rank -= count;
                }
            }
            bit >>= 1;
        }
        Some(before)
    }
}

mod rmq;
pub use rmq::RangeMinimumIndex;
