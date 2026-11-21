//! Rooted Sleator–Tarjan link–cut forest with capacities on parent edges.
use crate::NONE;

const INFINITE: u128 = u128::MAX;
#[derive(Clone)]
struct Node {
    // An auxiliary tree uses in-order root-to-leaf order. At its root,
    // parent denotes the represented edge above the exposed path.
    parent: usize,
    children: [usize; 2],
    capacity: u128,
    minimum: (u128, usize),
    debit: u128,
}

/// A changing rooted forest over vertex IDs 0..n. Every edge has a u64
/// capacity, including zero. Root queries, links, cuts, and path operations
/// take amortized O(log n) time. Invalid vertex IDs panic.
pub struct LinkCutForest {
    nodes: Vec<Node>,
    stack: Vec<usize>,
}

impl LinkCutForest {
    pub fn new(n: usize) -> Self {
        Self {
            nodes: (0..n)
                .map(|i| Node {
                    parent: NONE,
                    children: [NONE; 2],
                    capacity: INFINITE,
                    minimum: (INFINITE, i),
                    debit: 0,
                })
                .collect(),
            stack: vec![],
        }
    }
    fn auxiliary_root(&self, x: usize) -> bool {
        let p = self.nodes[x].parent;
        p == NONE || !self.nodes[p].children.contains(&x)
    }
    fn pull(&mut self, x: usize) {
        let mut minimum = (self.nodes[x].capacity, x);
        for child in self.nodes[x].children {
            if child != NONE {
                minimum = minimum.min(self.nodes[child].minimum);
            }
        }
        self.nodes[x].minimum = minimum;
    }
    fn debit(&mut self, x: usize, amount: u128) {
        if x == NONE || self.nodes[x].minimum.0 == INFINITE {
            return;
        }
        let node = &mut self.nodes[x];
        if node.capacity != INFINITE {
            node.capacity -= amount;
        }
        node.minimum.0 -= amount;
        node.debit += amount;
    }
    fn push(&mut self, x: usize) {
        let amount = std::mem::take(&mut self.nodes[x].debit);
        if amount != 0 {
            for child in self.nodes[x].children {
                self.debit(child, amount);
            }
        }
    }
    fn rotate(&mut self, x: usize) {
        let p = self.nodes[x].parent;
        let g = self.nodes[p].parent;
        let side = usize::from(self.nodes[p].children[1] == x);
        let middle = self.nodes[x].children[1 - side];
        if !self.auxiliary_root(p) {
            let at = usize::from(self.nodes[g].children[1] == p);
            self.nodes[g].children[at] = x;
        }
        self.nodes[x].parent = g;
        self.nodes[p].children[side] = middle;
        if middle != NONE {
            self.nodes[middle].parent = p;
        }
        self.nodes[x].children[1 - side] = p;
        self.nodes[p].parent = x;
        self.pull(p);
        self.pull(x);
    }
    fn splay(&mut self, x: usize) {
        self.stack.push(x);
        let mut v = x;
        while !self.auxiliary_root(v) {
            v = self.nodes[v].parent;
            self.stack.push(v);
        }
        while let Some(v) = self.stack.pop() {
            self.push(v);
        }
        while !self.auxiliary_root(x) {
            let p = self.nodes[x].parent;
            if !self.auxiliary_root(p) {
                let g = self.nodes[p].parent;
                let aligned = (self.nodes[p].children[1] == x) == (self.nodes[g].children[1] == p);
                self.rotate(if aligned { p } else { x });
            }
            self.rotate(x);
        }
    }
    fn expose(&mut self, x: usize) {
        let mut v = x;
        let mut suffix = NONE;
        while v != NONE {
            self.splay(v);
            // The detached right child retains v as its path parent.
            self.nodes[v].children[1] = suffix;
            self.pull(v);
            suffix = v;
            v = self.nodes[v].parent;
        }
        self.splay(x);
    }

    pub fn root(&mut self, x: usize) -> usize {
        self.expose(x);
        let mut v = x;
        while self.nodes[v].children[0] != NONE {
            self.push(v);
            v = self.nodes[v].children[0];
        }
        self.splay(v);
        v
    }

    /// Attaches a root to a node in another tree. Returns false for cycles
    /// and for children that already have a parent.
    pub fn link(&mut self, child: usize, parent: usize, capacity: u64) -> bool {
        if self.root(parent) == child {
            return false;
        }
        self.expose(child);
        if self.nodes[child].children[0] != NONE {
            return false;
        }
        self.nodes[child].capacity = u128::from(capacity);
        self.pull(child);
        self.nodes[child].parent = parent;
        true
    }

    /// Removes the parent edge and returns its remaining capacity.
    pub fn cut(&mut self, x: usize) -> Option<u64> {
        self.expose(x);
        let ancestors = self.nodes[x].children[0];
        if ancestors == NONE {
            return None;
        }
        self.nodes[ancestors].parent = NONE;
        self.nodes[x].children[0] = NONE;
        let amount = std::mem::replace(&mut self.nodes[x].capacity, INFINITE);
        self.pull(x);
        Some(amount as u64)
    }

    /// Minimum remaining capacity and the child owning the minimizing edge.
    pub fn path_minimum(&mut self, x: usize) -> Option<(u64, usize)> {
        self.expose(x);
        let (cap, id) = self.nodes[x].minimum;
        (cap != INFINITE).then_some((cap as u64, id))
    }

    /// Subtracts the same amount from each edge up to the root. Insufficient
    /// capacity leaves the forest unchanged. An empty path accepts any amount.
    pub fn subtract_path(&mut self, x: usize, amount: u64) -> bool {
        self.expose(x);
        if self.nodes[x].minimum.0 < u128::from(amount) {
            return false;
        }
        self.debit(x, u128::from(amount));
        true
    }
}
