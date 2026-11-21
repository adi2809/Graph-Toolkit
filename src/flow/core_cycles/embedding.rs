//! Shared walk expressions for forest contractions and parallel spanner edges.
//! Parents always have larger IDs. Batched repairs run in topological order.
use std::{
    cmp::{Ordering, Reverse},
    collections::BinaryHeap,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct Walk(usize);
impl Walk {
    pub(super) fn reverse(self) -> Self {
        Self(self.0 ^ 1)
    }
    fn index(self) -> usize {
        self.0 >> 1
    }
    fn backwards(self) -> bool {
        self.0 & 1 == 1
    }
    fn oriented(self, back: bool) -> Self {
        if back { self.reverse() } else { self }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct Metric {
    pub(super) gradient: i128,
    pub(super) length: u128,
}
pub(super) const CLOSED: Metric = Metric {
    gradient: 0,
    length: u128::MAX,
};

impl Metric {
    pub(super) fn arc(gradient: i64, length: u64, forward: bool, backward: bool) -> [Self; 2] {
        let g = i128::from(gradient);
        let l = u128::from(length);
        [
            if forward {
                Self {
                    gradient: g,
                    length: l,
                }
            } else {
                CLOSED
            },
            if backward {
                Self {
                    gradient: -g,
                    length: l,
                }
            } else {
                CLOSED
            },
        ]
    }
    pub(super) fn open(self) -> bool {
        self != CLOSED
    }
    fn join(self, tail: Self) -> Self {
        if !self.open() || !tail.open() {
            return CLOSED;
        }
        Self {
            gradient: self.gradient + tail.gradient,
            length: self.length + tail.length,
        }
    }
    pub(super) fn compare_ratio(self, other: Self) -> Ordering {
        match (self.open(), other.open()) {
            (false, false) => return Ordering::Equal,
            (false, true) => return Ordering::Greater,
            (true, false) => return Ordering::Less,
            _ => {}
        }
        let sign = self.gradient.signum().cmp(&other.gradient.signum());
        if sign != Ordering::Equal {
            return sign;
        }
        let order = compare_fractions(
            self.gradient.unsigned_abs(),
            self.length,
            other.gradient.unsigned_abs(),
            other.length,
        );
        if self.gradient < 0 {
            order.reverse()
        } else {
            order
        }
    }
}
fn compare_fractions(mut a: u128, mut b: u128, mut c: u128, mut d: u128) -> Ordering {
    if let (Some(x), Some(y)) = (a.checked_mul(d), c.checked_mul(b)) {
        return x.cmp(&y);
    }
    // Continued fractions avoid overflowing a cross product.
    let mut reverse = false;
    loop {
        let mut order = (a / b).cmp(&(c / d));
        let (r, s) = (a % b, c % d);
        if order == Ordering::Equal && (r == 0 || s == 0) {
            order = r.cmp(&s);
        }
        if order != Ordering::Equal || r == 0 || s == 0 {
            return if reverse { order.reverse() } else { order };
        }
        (a, b, c, d) = (b, r, d, s);
        reverse = !reverse;
    }
}

#[derive(Clone, Copy)]
enum Expression {
    Edge(usize),
    Join(Walk, Walk),
    // An embedded parallel spanner uses its shortest available witness in
    // each direction. A contest chooses the better ratio among closed walks.
    Parallel(Walk, Walk),
    Contest(Walk, Walk),
}
impl Expression {
    fn children(self) -> Option<[Walk; 2]> {
        match self {
            Self::Edge(_) => None,
            Self::Join(a, b) | Self::Parallel(a, b) | Self::Contest(a, b) => Some([a, b]),
        }
    }
}
struct Node {
    expr: Expression,
    metric: [Metric; 2],
}

#[derive(Default)]
pub(super) struct Embedding {
    nodes: Vec<Node>,
    offsets: Vec<usize>,
    parents: Vec<usize>,
    pending: BinaryHeap<Reverse<usize>>,
    dirty: Vec<bool>,
}
impl Embedding {
    pub(super) fn metric(&self, w: Walk) -> Metric {
        self.nodes[w.index()].metric[usize::from(w.backwards())]
    }
    fn insert(&mut self, expr: Expression, metric: [Metric; 2]) -> Walk {
        let id = self.nodes.len();
        self.nodes.push(Node { expr, metric });
        Walk(id * 2)
    }
    pub(super) fn edge(&mut self, id: usize, metric: [Metric; 2]) -> Walk {
        self.insert(Expression::Edge(id), metric)
    }
    fn choose(&self, a: Walk, b: Walk, ratio: bool) -> Walk {
        let (x, y) = (self.metric(a), self.metric(b));
        let order = if ratio {
            x.compare_ratio(y)
        } else {
            x.length.cmp(&y.length)
        };
        if order.is_gt() { b } else { a }
    }
    fn evaluate(&self, expr: Expression, back: bool) -> Metric {
        match expr {
            Expression::Join(a, b) => self
                .metric(a.oriented(back))
                .join(self.metric(b.oriented(back))),
            Expression::Parallel(a, b) | Expression::Contest(a, b) => self.metric(self.choose(
                a.oriented(back),
                b.oriented(back),
                matches!(expr, Expression::Contest(..)),
            )),
            Expression::Edge(_) => unreachable!(),
        }
    }
    fn branch(&mut self, expr: Expression) -> Walk {
        let metric = [self.evaluate(expr, false), self.evaluate(expr, true)];
        self.insert(expr, metric)
    }
    pub(super) fn join(&mut self, a: Walk, b: Walk) -> Walk {
        self.branch(Expression::Join(a, b))
    }
    pub(super) fn parallel(&mut self, a: Walk, b: Walk) -> Walk {
        self.branch(Expression::Parallel(a, b))
    }
    pub(super) fn contest(&mut self, mut cycles: Vec<Walk>) -> Option<Walk> {
        while cycles.len() > 1 {
            let mut next = Vec::with_capacity(cycles.len().div_ceil(2));
            for pair in cycles.chunks(2) {
                next.push(match pair {
                    [a, b] => self.branch(Expression::Contest(*a, *b)),
                    [a] => *a,
                    _ => unreachable!(),
                });
            }
            cycles = next;
        }
        cycles.pop()
    }
    pub(super) fn better_orientation(&self, w: Walk) -> Walk {
        self.choose(w, w.reverse(), true)
    }
    pub(super) fn seal(&mut self) {
        let n = self.nodes.len();
        self.offsets = vec![0; n + 1];
        self.dirty = vec![false; n];
        for node in &self.nodes {
            if let Some(children) = node.expr.children() {
                for w in children {
                    self.offsets[w.index() + 1] += 1;
                }
            }
        }
        for i in 0..n {
            self.offsets[i + 1] += self.offsets[i];
        }
        let mut cursor = self.offsets.clone();
        self.parents = vec![0; self.offsets[n]];
        for (i, node) in self.nodes.iter().enumerate() {
            if let Some(children) = node.expr.children() {
                for w in children {
                    self.parents[cursor[w.index()]] = i;
                    cursor[w.index()] += 1;
                }
            }
        }
    }
    fn notify(&mut self, id: usize) {
        for &p in &self.parents[self.offsets[id]..self.offsets[id + 1]] {
            if !std::mem::replace(&mut self.dirty[p], true) {
                self.pending.push(Reverse(p));
            }
        }
    }
    pub(super) fn update(&mut self, leaf: Walk, metric: [Metric; 2]) {
        let i = leaf.index();
        if self.nodes[i].metric != metric {
            self.nodes[i].metric = metric;
            self.notify(i);
        }
    }
    pub(super) fn repair(&mut self) {
        while let Some(Reverse(i)) = self.pending.pop() {
            self.dirty[i] = false;
            let expr = self.nodes[i].expr;
            let metric = [self.evaluate(expr, false), self.evaluate(expr, true)];
            if self.nodes[i].metric != metric {
                self.nodes[i].metric = metric;
                self.notify(i);
            }
        }
    }
    pub(super) fn mark_leaves(&self, marked: impl Fn(usize) -> bool) -> Vec<bool> {
        let mut output = Vec::with_capacity(self.nodes.len());
        for node in &self.nodes {
            let flag = match node.expr {
                Expression::Edge(id) => marked(id),
                expr => {
                    let [a, b] = expr.children().unwrap();
                    output[a.index()] || output[b.index()]
                }
            };
            output.push(flag);
        }
        output
    }

    // A terminal walk contains exactly one crossing. Skip each entire
    // tree segment on the way to that crossing's original directed edge.
    pub(super) fn marked_leaf(&self, mut w: Walk, marked: &[bool]) -> (usize, bool) {
        loop {
            let back = w.backwards();
            w = match self.nodes[w.index()].expr {
                Expression::Edge(id) => return (id, back),
                Expression::Join(a, b) => (if marked[a.index()] { a } else { b }).oriented(back),
                Expression::Parallel(a, b) => {
                    self.choose(a.oriented(back), b.oriented(back), false)
                }
                Expression::Contest(..) => {
                    unreachable!("terminal routes have no internal contests")
                }
            };
        }
    }

    pub(super) fn lift(&self, walk: Walk, out: &mut Vec<(usize, bool)>) {
        out.clear();
        let mut stack = vec![walk];
        while let Some(w) = stack.pop() {
            let back = w.backwards();
            match self.nodes[w.index()].expr {
                Expression::Edge(id) => out.push((id, back)),
                Expression::Join(a, b) => {
                    if back {
                        stack.extend([a.reverse(), b.reverse()]);
                    } else {
                        stack.extend([b, a]);
                    }
                }
                Expression::Parallel(a, b) => {
                    stack.push(self.choose(a.oriented(back), b.oriented(back), false))
                }
                Expression::Contest(a, b) => {
                    stack.push(self.choose(a.oriented(back), b.oriented(back), true))
                }
            }
        }
    }
}
