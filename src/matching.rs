//! Exact cardinality and weighted matching, bipartite assignment, and edge covers.

use std::error::Error;
use std::fmt;

use crate::connectivity::DisjointSet;
use crate::{Arc, EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

mod assignment;
mod bottleneck;
mod gabow;
mod tree_union;
mod weighted;
pub use assignment::{
    bipartite_coloring, minimum_cost_bipartite_matching, minimum_cost_full_bipartite_matching,
};
pub use bottleneck::{
    BottleneckMatching, minimum_bottleneck_matching, minimum_bottleneck_matching_of_size,
    minimum_bottleneck_perfect_matching,
};
pub use gabow::maximum_general_matching_gabow;
pub use weighted::{
    WeightedMatching, maximum_weight_matching, minimum_cost_matching, minimum_cost_perfect_matching,
};

/// Left and right identifiers each start at zero.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum BipartiteError {
    InvalidLeft { vertex: usize, count: usize },
    InvalidRight { vertex: usize, count: usize },
}

impl fmt::Display for BipartiteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::InvalidLeft { vertex, count } => {
                write!(f, "left vertex {vertex} is outside 0..{count}")
            }
            Self::InvalidRight { vertex, count } => {
                write!(f, "right vertex {vertex} is outside 0..{count}")
            }
        }
    }
}

impl Error for BipartiteError {}

/// CSR bipartite graph with adjacency stored only for the left partition.
#[derive(Clone, Debug)]
pub struct BipartiteGraph {
    offsets: Vec<usize>,
    right: Vec<usize>,
    right_count: usize,
}

impl BipartiteGraph {
    /// Parallel edges are accepted. Neighbor order follows the input iterator.
    pub fn from_edges(
        left_count: usize,
        right_count: usize,
        edges: impl IntoIterator<Item = (usize, usize)>,
    ) -> Result<Self, BipartiteError> {
        let mut offsets = vec![0; left_count + 1];
        let mut pairs = Vec::new();
        for (left, right) in edges {
            if left >= left_count {
                return Err(BipartiteError::InvalidLeft {
                    vertex: left,
                    count: left_count,
                });
            }
            if right >= right_count {
                return Err(BipartiteError::InvalidRight {
                    vertex: right,
                    count: right_count,
                });
            }
            offsets[left + 1] += 1;
            pairs.push((left, right));
        }
        for node in 0..left_count {
            offsets[node + 1] += offsets[node];
        }
        let mut cursor = offsets[..left_count].to_vec();
        let mut right = vec![0; pairs.len()];
        for (left, node) in pairs {
            right[cursor[left]] = node;
            cursor[left] += 1;
        }
        Ok(Self {
            offsets,
            right,
            right_count,
        })
    }

    pub fn left_count(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn right_count(&self) -> usize {
        self.right_count
    }

    fn neighbors(&self, left: usize) -> &[usize] {
        &self.right[self.offsets[left]..self.offsets[left + 1]]
    }

    /// Deterministic Hopcroft–Karp with iterative blocking-path searches.
    /// After initialization, the matching phase takes O(E √V) time.
    pub fn maximum_matching(&self) -> BipartiteMatching {
        self.solve(None)
    }

    /// Randomizes the greedy warm start, then uses Hopcroft–Karp to finish exactly.
    /// The seed makes the chosen matching reproducible.
    pub fn randomized_maximum_matching(&self, seed: u64) -> BipartiteMatching {
        self.solve(Some(SplitMix64(seed)))
    }

    fn solve(&self, mut rng: Option<SplitMix64>) -> BipartiteMatching {
        let mut matching = BipartiteMatching {
            left: vec![NONE; self.left_count()],
            right: vec![NONE; self.right_count],
            cover: BipartiteVertexCover::default(),
            len: 0,
        };
        let mut order: Vec<_> = (0..self.left_count())
            .filter(|&u| !self.neighbors(u).is_empty())
            .collect();
        if let Some(rng) = &mut rng {
            for index in (1..order.len()).rev() {
                order.swap(index, rng.below(index + 1));
            }
        }
        for &u in &order {
            let neighbors = self.neighbors(u);
            let offset = rng.as_mut().map_or(0, |rng| rng.below(neighbors.len()));
            let candidates = neighbors[offset..].iter().chain(&neighbors[..offset]);
            if let Some(&v) = candidates.into_iter().find(|&&v| matching.right[v] == NONE) {
                matching.pair(u, v);
                matching.len += 1;
            }
        }

        let mut search = Layers {
            level: vec![NONE; self.left_count()],
            cursor: vec![0; self.left_count()],
            queue: Vec::new(),
            stack: Vec::new(),
        };
        while let Some(shortest) = search.build(self, &matching, &order) {
            for &left in &order {
                if matching.left[left] == NONE
                    && search.level[left] == 0
                    && search.augment(self, &mut matching, left, shortest)
                {
                    matching.len += 1;
                }
            }
        }
        // The final unsuccessful BFS gives the alternating reachable left vertices.
        let mut reached_right = vec![false; self.right_count];
        for &left in &search.queue {
            for &right in self.neighbors(left) {
                reached_right[right] = true;
            }
        }
        matching.cover = BipartiteVertexCover {
            left: (0..self.left_count())
                .filter(|&u| !self.neighbors(u).is_empty() && search.level[u] == NONE)
                .collect(),
            right: (0..self.right_count)
                .filter(|&v| reached_right[v])
                .collect(),
        };
        matching
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct BipartiteVertexCover {
    pub left: Vec<usize>,
    pub right: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct BipartiteMatching {
    left: Vec<usize>,
    right: Vec<usize>,
    cover: BipartiteVertexCover,
    len: usize,
}

impl BipartiteMatching {
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn right_of(&self, left: usize) -> Option<usize> {
        self.left.get(left).copied().filter(|&mate| mate != NONE)
    }

    pub fn left_of(&self, right: usize) -> Option<usize> {
        self.right.get(right).copied().filter(|&mate| mate != NONE)
    }

    pub fn pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.left
            .iter()
            .copied()
            .enumerate()
            .filter(|&(_, mate)| mate != NONE)
    }

    /// A minimum vertex cover obtained by Kőnig's theorem. Its size equals `len()`.
    pub fn minimum_vertex_cover(&self) -> &BipartiteVertexCover {
        &self.cover
    }

    fn pair(&mut self, left: usize, right: usize) {
        self.left[left] = right;
        self.right[right] = left;
    }
}

struct Layers {
    level: Vec<usize>,
    cursor: Vec<usize>,
    queue: Vec<usize>,
    stack: Vec<(usize, usize)>,
}

impl Layers {
    fn build(
        &mut self,
        graph: &BipartiteGraph,
        matching: &BipartiteMatching,
        order: &[usize],
    ) -> Option<usize> {
        self.queue.clear();
        for &left in order {
            self.cursor[left] = 0;
            self.level[left] = NONE;
            if matching.left[left] == NONE {
                self.level[left] = 0;
                self.queue.push(left);
            }
        }
        let mut shortest = NONE;
        let mut head = 0;
        while let Some(&left) = self.queue.get(head) {
            head += 1;
            let next_level = self.level[left] + 1;
            if next_level > shortest {
                break;
            }
            for &right in graph.neighbors(left) {
                let mate = matching.right[right];
                if mate == NONE {
                    shortest = next_level;
                } else if next_level < shortest && self.level[mate] == NONE {
                    self.level[mate] = next_level;
                    self.queue.push(mate);
                }
            }
        }
        (shortest != NONE).then_some(shortest)
    }

    fn augment(
        &mut self,
        graph: &BipartiteGraph,
        matching: &mut BipartiteMatching,
        root: usize,
        shortest: usize,
    ) -> bool {
        self.stack.push((root, NONE));
        while let Some(&(left, _)) = self.stack.last() {
            let cursor = &mut self.cursor[left];
            let Some(&right) = graph.neighbors(left).get(*cursor) else {
                self.level[left] = NONE;
                self.stack.pop();
                continue;
            };
            *cursor += 1;
            let mate = matching.right[right];
            let next_level = self.level[left] + 1;
            if mate == NONE && next_level == shortest {
                let mut right = right;
                while let Some((left, previous_right)) = self.stack.pop() {
                    self.level[left] = NONE;
                    matching.pair(left, right);
                    right = previous_right;
                }
                return true;
            }
            if mate != NONE && next_level < shortest && self.level[mate] == next_level {
                self.stack.push((mate, right));
            }
        }
        false
    }
}

struct SplitMix64(u64);

impl SplitMix64 {
    fn below(&mut self, upper: usize) -> usize {
        let bound = upper as u64;
        loop {
            self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
            let mut x = self.0;
            x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
            x ^= x >> 31;
            let product = u128::from(x) * u128::from(bound);
            let low = product as u64;
            if low >= bound || low >= bound.wrapping_neg() % bound {
                return (product >> 64) as usize;
            }
        }
    }
}

#[derive(Clone, Debug)]
pub struct GeneralMatching {
    mates: Vec<Arc>,
    len: usize,
}

const VACANT: Arc = Arc {
    node: NONE,
    edge: NONE,
};

impl GeneralMatching {
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Returns the matched neighbor and the logical edge ID used.
    pub fn mate(&self, node: NodeId) -> Option<Arc> {
        self.mates
            .get(node)
            .copied()
            .filter(|entry| *entry != VACANT)
    }

    pub fn edges(&self) -> impl Iterator<Item = EdgeId> + '_ {
        self.mates
            .iter()
            .enumerate()
            .filter_map(|(node, arc)| (arc.node != NONE && node < arc.node).then_some(arc.edge))
    }

    fn pair(&mut self, a: NodeId, b: NodeId, edge: EdgeId) {
        self.mates[a] = Arc { node: b, edge };
        self.mates[b] = Arc { node: a, edge };
    }
}

/// Greedy maximal matching in O(V + E) time. Self-loops are ignored.
pub fn maximal_matching<W>(graph: &Graph<W>) -> Result<GeneralMatching, GraphError> {
    graph.require(GraphKind::Undirected)?;
    Ok(greedy_matching(graph))
}

/// Tests whether every ID belongs to a loop-free, vertex-disjoint set of edges.
pub fn is_matching<W>(graph: &Graph<W>, edges: &[EdgeId]) -> Result<bool, GraphError> {
    graph.require(GraphKind::Undirected)?;
    Ok(matched_vertices(graph, edges).is_some())
}

/// Tests whether the matching covers every vertex.
pub fn is_perfect_matching<W>(graph: &Graph<W>, edges: &[EdgeId]) -> Result<bool, GraphError> {
    graph.require(GraphKind::Undirected)?;
    Ok(edges.len().checked_mul(2) == Some(graph.node_count())
        && matched_vertices(graph, edges).is_some())
}

/// Tests whether any additional non-loop edge could be added to the matching.
pub fn is_maximal_matching<W>(graph: &Graph<W>, edges: &[EdgeId]) -> Result<bool, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let Some(covered) = matched_vertices(graph, edges) else {
        return Ok(false);
    };
    Ok(graph
        .edges()
        .iter()
        .all(|edge| edge.from == edge.to || covered[edge.from] || covered[edge.to]))
}

fn matched_vertices<W>(graph: &Graph<W>, edges: &[EdgeId]) -> Option<Vec<bool>> {
    let mut covered = vec![false; graph.node_count()];
    for &id in edges {
        let edge = graph.edge(id)?;
        if edge.from == edge.to || covered[edge.from] || covered[edge.to] {
            return None;
        }
        covered[edge.from] = true;
        covered[edge.to] = true;
    }
    Some(covered)
}

/// Tests whether every vertex is incident to at least one supplied edge.
/// Loops cover their incident vertex and duplicate edge IDs are allowed.
pub fn is_edge_cover<W>(graph: &Graph<W>, edges: &[EdgeId]) -> bool {
    let mut covered = vec![false; graph.node_count()];
    for &id in edges {
        let Some(edge) = graph.edge(id) else {
            return false;
        };
        covered[edge.from] = true;
        covered[edge.to] = true;
    }
    covered.into_iter().all(|vertex| vertex)
}

/// Returns a smallest set of incident edges covering every vertex.
/// An isolated vertex makes the cover infeasible. A loop can cover its vertex.
pub fn minimum_edge_cover<W>(graph: &Graph<W>) -> Result<Option<Vec<EdgeId>>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    if (0..graph.node_count()).any(|node| graph.out(node).is_empty()) {
        return Ok(None);
    }
    let matching = maximum_general_matching_gabow(graph)?;
    let mut cover: Vec<_> = matching.edges().collect();
    for node in 0..graph.node_count() {
        if matching.mate(node).is_none() {
            cover.push(graph.out(node)[0].edge);
        }
    }
    Ok(Some(cover))
}

/// Edmonds' blossom algorithm with a greedy warm start, iterative path lifting,
/// and union-find contraction. Edge payloads and self-loops are ignored.
pub fn maximum_general_matching<W>(graph: &Graph<W>) -> Result<GeneralMatching, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    let mut matching = greedy_matching(graph);
    let mut search = BlossomSearch::new(n);
    for root in 0..n {
        if matching.mates[root] == VACANT && graph.out(root).iter().any(|arc| arc.node != root) {
            search.augment(graph, &mut matching, root);
        }
    }
    Ok(matching)
}

fn greedy_matching<W>(graph: &Graph<W>) -> GeneralMatching {
    let n = graph.node_count();
    let mut matching = GeneralMatching {
        mates: vec![VACANT; n],
        len: 0,
    };
    for node in 0..n {
        if matching.mates[node] != VACANT {
            continue;
        }
        if let Some(arc) = graph
            .out(node)
            .iter()
            .find(|arc| arc.node != node && matching.mates[arc.node] == VACANT)
        {
            matching.pair(node, arc.node, arc.edge);
            matching.len += 1;
        }
    }
    matching
}

struct BlossomSearch {
    sets: DisjointSet,
    base: Vec<NodeId>,
    parent: Vec<Arc>,
    outer: Vec<bool>,
    queue: Vec<NodeId>,
    touched: Vec<NodeId>,
    contracted: Vec<NodeId>,
    mark: Vec<usize>,
    epoch: usize,
}

impl BlossomSearch {
    fn new(n: usize) -> Self {
        Self {
            sets: DisjointSet::new(n),
            base: (0..n).collect(),
            parent: vec![VACANT; n],
            outer: vec![false; n],
            queue: Vec::new(),
            touched: Vec::new(),
            contracted: Vec::new(),
            mark: vec![0; n],
            epoch: 0,
        }
    }

    fn augment<W>(&mut self, graph: &Graph<W>, matching: &mut GeneralMatching, root: NodeId) {
        self.outer[root] = true;
        self.queue.push(root);
        self.touched.push(root);
        if let Some(end) = self.find_endpoint(graph, matching) {
            let mut node = end;
            while self.parent[node] != VACANT {
                let step = self.parent[node];
                let displaced = matching.mates[step.node].node;
                matching.pair(node, step.node, step.edge);
                if displaced == NONE {
                    break;
                }
                node = displaced;
            }
            matching.len += 1;
        }
        for node in self.touched.drain(..) {
            self.sets.reset(node);
            self.base[node] = node;
            self.parent[node] = VACANT;
            self.outer[node] = false;
        }
        self.queue.clear();
    }

    fn find_endpoint<W>(&mut self, graph: &Graph<W>, matching: &GeneralMatching) -> Option<NodeId> {
        let mut head = 0;
        while let Some(&node) = self.queue.get(head) {
            head += 1;
            for &arc in graph.out(node) {
                let next = arc.node;
                if self.base_of(node) == self.base_of(next) || matching.mates[node].node == next {
                    continue;
                }
                if self.outer[next] {
                    let base = self.common_base(node, next, matching);
                    self.mark_branch(node, base, arc, matching);
                    self.mark_branch(
                        next,
                        base,
                        Arc {
                            node,
                            edge: arc.edge,
                        },
                        matching,
                    );
                    for member in self.contracted.drain(..) {
                        if let Some(representative) = self.sets.union(base, member) {
                            self.base[representative] = base;
                        }
                    }
                } else if self.parent[next] == VACANT {
                    self.parent[next] = Arc {
                        node,
                        edge: arc.edge,
                    };
                    self.touched.push(next);
                    let mate = matching.mates[next].node;
                    if mate == NONE {
                        return Some(next);
                    }
                    self.outer[mate] = true;
                    self.touched.push(mate);
                    self.queue.push(mate);
                }
            }
        }
        None
    }

    fn base_of(&mut self, node: NodeId) -> NodeId {
        self.base[self.sets.find(node)]
    }

    fn common_base(&mut self, a: NodeId, b: NodeId, matching: &GeneralMatching) -> NodeId {
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.mark.fill(0);
            self.epoch = 1;
        }
        let (mut a, mut b) = (a, b);
        loop {
            if a != NONE {
                let base = self.base_of(a);
                if self.mark[base] == self.epoch {
                    return base;
                }
                self.mark[base] = self.epoch;
                let mate = matching.mates[base].node;
                a = if mate == NONE {
                    NONE
                } else {
                    self.parent[mate].node
                };
            }
            std::mem::swap(&mut a, &mut b);
        }
    }

    fn mark_branch(
        &mut self,
        mut node: NodeId,
        base: NodeId,
        mut bridge: Arc,
        matching: &GeneralMatching,
    ) {
        while self.base_of(node) != base {
            let mate = matching.mates[node].node;
            self.contracted.extend([node, mate]);
            if !self.outer[mate] {
                self.outer[mate] = true;
                self.queue.push(mate);
            }
            self.parent[node] = bridge;
            let step = self.parent[mate];
            bridge = Arc {
                node: mate,
                edge: step.edge,
            };
            node = step.node;
        }
    }
}
