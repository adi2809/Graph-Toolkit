use crate::NONE;

struct Node<K> {
    key: Option<K>,
    parent: usize,
    first_child: usize,
    previous: usize,
    next: usize,
    degree: usize,
    marked: bool,
}
impl<K> Default for Node<K> {
    fn default() -> Self {
        Self {
            key: None,
            parent: NONE,
            first_child: NONE,
            previous: NONE,
            next: NONE,
            degree: 0,
            marked: false,
        }
    }
}

/// Indexed Fibonacci min-heap. IDs in 0..capacity can be reused after extraction.
/// Insert and decrease-key take O(1) amortized time. Extract-min takes O(log n)
/// amortized time, where n is the live heap size.
pub struct FibonacciHeap<K> {
    nodes: Vec<Node<K>>,
    roots: Vec<usize>,
    minimum: usize,
    len: usize,
}
impl<K: Copy + Ord> FibonacciHeap<K> {
    pub fn new(capacity: usize) -> Self {
        Self {
            nodes: (0..capacity).map(|_| Node::default()).collect(),
            roots: vec![],
            minimum: NONE,
            len: 0,
        }
    }
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    pub fn peek(&self) -> Option<(K, usize)> {
        Some((self.nodes.get(self.minimum)?.key?, self.minimum))
    }
    /// Inserts an absent ID or strictly decreases its key. Returns false for
    /// an out-of-range ID or a key that would leave the priority unchanged or larger.
    pub fn decrease(&mut self, id: usize, key: K) -> bool {
        let Some(node) = self.nodes.get_mut(id) else {
            return false;
        };
        if node.key.is_some_and(|old| key >= old) {
            return false;
        }
        let new = node.key.replace(key).is_none();
        if new {
            self.len += 1;
            self.roots.push(id);
        } else {
            let mut v = id;
            let mut p = self.nodes[v].parent;
            if p != NONE && self.nodes[v].key < self.nodes[p].key {
                loop {
                    self.cut(v, p);
                    v = p;
                    p = self.nodes[v].parent;
                    if p == NONE {
                        break;
                    }
                    if !self.nodes[v].marked {
                        self.nodes[v].marked = true;
                        break;
                    }
                }
            }
        }
        if self.minimum == NONE || self.nodes[id].key < self.nodes[self.minimum].key {
            self.minimum = id;
        }
        true
    }

    fn cut(&mut self, v: usize, p: usize) {
        let (a, b) = (self.nodes[v].previous, self.nodes[v].next);
        if a == NONE {
            self.nodes[p].first_child = b;
        } else {
            self.nodes[a].next = b;
        }
        if b != NONE {
            self.nodes[b].previous = a;
        }
        self.nodes[p].degree -= 1;
        self.make_root(v);
    }
    fn make_root(&mut self, v: usize) {
        let node = &mut self.nodes[v];
        node.parent = NONE;
        node.marked = false;
        node.previous = NONE;
        node.next = NONE;
        self.roots.push(v);
    }
    fn attach(&mut self, child: usize, parent: usize) {
        let head = self.nodes[parent].first_child;
        self.nodes[child].parent = parent;
        self.nodes[child].next = head;
        if head != NONE {
            self.nodes[head].previous = child;
        }
        self.nodes[parent].first_child = child;
        self.nodes[parent].degree += 1;
    }
    pub fn pop(&mut self) -> Option<(K, usize)> {
        let result = self.peek()?;
        let old = self.minimum;
        let mut c = self.nodes[old].first_child;
        while c != NONE {
            let next = self.nodes[c].next;
            self.make_root(c);
            c = next;
        }
        self.nodes[old] = Node::default();
        self.minimum = NONE;
        self.len -= 1;
        let mut ranks = Vec::new();
        while let Some(mut root) = self.roots.pop() {
            if root == old {
                continue;
            }
            loop {
                let degree = self.nodes[root].degree;
                if degree >= ranks.len() {
                    ranks.resize(degree + 1, NONE);
                }
                let mut other = std::mem::replace(&mut ranks[degree], NONE);
                if other == NONE {
                    ranks[degree] = root;
                    break;
                }
                if self.nodes[other].key < self.nodes[root].key {
                    std::mem::swap(&mut root, &mut other);
                }
                self.attach(other, root);
            }
        }
        for root in ranks.into_iter().filter(|&r| r != NONE) {
            self.roots.push(root);
            if self.minimum == NONE || self.nodes[root].key < self.nodes[self.minimum].key {
                self.minimum = root;
            }
        }
        Some(result)
    }
}
