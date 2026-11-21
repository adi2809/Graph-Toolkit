//! Immutable range minima. Within a machine-word block, a monotone stack of
//! candidate minima fits in one word. Prefix snapshots allow any subrange of
//! that block to be answered with one mask and a trailing-zero count.
use std::ops::Range;

pub struct RangeMinimumIndex<T> {
    values: Vec<T>,
    candidates: Vec<usize>,
    levels: Vec<Vec<usize>>,
}
const BLOCK: usize = usize::BITS as usize;
impl<T: Ord> RangeMinimumIndex<T> {
    pub fn new(values: Vec<T>) -> Self {
        let mut candidates = Vec::with_capacity(values.len());
        let mut minima = Vec::new();
        for (block, chunk) in values.chunks(BLOCK).enumerate() {
            let mut stack = 0_usize;
            for (i, value) in chunk.iter().enumerate() {
                while stack != 0 {
                    let top = stack.ilog2() as usize;
                    if chunk[top] <= *value {
                        break;
                    }
                    stack ^= 1 << top;
                }
                stack |= 1 << i;
                candidates.push(stack);
            }
            minima.push(block * BLOCK + stack.trailing_zeros() as usize);
        }
        let mut levels = vec![minima];
        let mut span = 1;
        while span <= levels[0].len() / 2 {
            let previous = levels.last().unwrap();
            levels.push(
                previous
                    .iter()
                    .zip(&previous[span..])
                    .map(|(&a, &b)| if values[a] <= values[b] { a } else { b })
                    .collect(),
            );
            span *= 2;
        }
        Self {
            values,
            candidates,
            levels,
        }
    }

    fn within(&self, start: usize, last: usize) -> usize {
        let mask = self.candidates[last] & (usize::MAX << (start % BLOCK));
        last / BLOCK * BLOCK + mask.trailing_zeros() as usize
    }
    /// Leftmost minimum in a nonempty, half-open range. O(1) per query.
    pub fn arg_min(&self, range: Range<usize>) -> Option<usize> {
        if range.is_empty() || range.end > self.values.len() {
            return None;
        }
        let first = range.start;
        let last = range.end - 1;
        let (a, b) = (first / BLOCK, last / BLOCK);
        if a == b {
            return Some(self.within(first, last));
        }
        let earlier = |x: usize, y: usize| {
            if (&self.values[x], x) <= (&self.values[y], y) {
                x
            } else {
                y
            }
        };
        let mut best = earlier(
            self.within(first, (a + 1) * BLOCK - 1),
            self.within(b * BLOCK, last),
        );
        let whole = b - a - 1;
        if whole != 0 {
            let k = whole.ilog2() as usize;
            best = earlier(
                best,
                earlier(self.levels[k][a + 1], self.levels[k][b - (1 << k)]),
            );
        }
        Some(best)
    }
    pub fn minimum(&self, range: Range<usize>) -> Option<&T> {
        self.arg_min(range).map(|index| &self.values[index])
    }
}
