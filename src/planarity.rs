//! Left–right planarity testing and combinatorial embeddings.

use crate::{Arc, EdgeId, Graph, GraphBuilder, GraphError, GraphKind, NONE, NodeId};

mod witness;
pub use witness::{
    KuratowskiKind, KuratowskiWitness, check_planarity_with_witness, kuratowski_witness,
};

/// Clockwise cyclic orders for the underlying simple undirected graph.
/// Each pair keeps the first original edge ID. Self-loops are discarded.
#[derive(Clone, Debug)]
pub struct PlanarEmbedding {
    offsets: Vec<usize>,
    rotation: Vec<Arc>,
    face_next: Vec<usize>,
}

impl PlanarEmbedding {
    pub fn node_count(&self) -> usize {
        self.offsets.len() - 1
    }

    pub fn edge_count(&self) -> usize {
        self.rotation.len() / 2
    }

    pub fn clockwise_neighbors(&self, node: NodeId) -> Result<&[Arc], GraphError> {
        crate::graph::validate_node(node, self.node_count())?;
        Ok(&self.rotation[self.offsets[node]..self.offsets[node + 1]])
    }

    /// Each face boundary walk, without repeating its first vertex at the
    /// end. Bridges occur twice on a walk. For a disconnected graph, the
    /// boundary of the outer face is represented separately per component.
    /// Isolated vertices have no edge boundary.
    pub fn face_boundaries(&self) -> Vec<Vec<NodeId>> {
        let mut seen = vec![false; self.rotation.len()];
        let mut faces = Vec::new();
        for dart in 0..self.rotation.len() {
            if seen[dart] {
                continue;
            }
            let mut walk = Vec::new();
            let mut current = dart;
            while !std::mem::replace(&mut seen[current], true) {
                walk.push(self.rotation[current].node);
                current = self.face_next[current];
            }
            faces.push(walk);
        }
        faces
    }
}

/// Tests planarity of the underlying simple undirected graph in O(V + E) time
/// and space. The iterative implementation handles deep DFS trees.
pub fn is_planar<W>(graph: &Graph<W>) -> bool {
    if graph.node_count() < 5 || graph.edge_count() < 9 {
        return true;
    }
    let Some(simple) = simple_projection(graph, false) else {
        return false;
    };
    is_simple_planar(&simple)
}

fn is_simple_planar<W>(simple: &Graph<W>) -> bool {
    if simple.node_count() < 5 || simple.edge_count() < 9 {
        return true;
    }
    if simple.edge_count() > planar_edge_bound(simple.node_count()) {
        return false;
    }
    let orientation = Orientation::new(simple);
    let ordered = orientation.order(false, &[]);
    Constraints::new(&orientation, simple.edge_count()).test(&ordered)
}

/// Returns a planar rotation system, or None if the graph is nonplanar.
/// Like `is_planar`, this uses the underlying simple undirected graph and
/// requires O(V + E) time and space.
pub fn check_planarity<W>(graph: &Graph<W>) -> Option<PlanarEmbedding> {
    let simple = simple_projection(graph, false)?;
    embed_simple(&simple)
}

fn embed_simple(simple: &Graph<EdgeId>) -> Option<PlanarEmbedding> {
    if simple.edge_count() > planar_edge_bound(simple.node_count()) {
        return None;
    }
    let orientation = Orientation::new(simple);
    let ordered = orientation.order(false, &[]);
    let mut constraints = Constraints::new(&orientation, simple.edge_count());
    if !constraints.test(&ordered) {
        return None;
    }
    constraints.resolve_sides();
    let ordered = orientation.order(true, &constraints.left);
    Some(orientation.embed(simple, &ordered, &constraints.left))
}

fn planar_edge_bound(n: usize) -> usize {
    if n >= 3 {
        3 * n - 6
    } else {
        n.saturating_sub(1)
    }
}

fn simple_projection<W>(graph: &Graph<W>, keep_dense_prefix: bool) -> Option<Graph<EdgeId>> {
    let n = graph.node_count();
    let bound = planar_edge_bound(n);
    let limit = bound + usize::from(keep_dense_prefix);
    let mut builder =
        GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, limit.min(graph.edge_count()));
    let mut best = vec![NONE; n];
    let mut touched = Vec::new();
    let mut count = 0;
    for u in 0..n {
        let incoming = if graph.is_directed() {
            graph.inc(u)
        } else {
            &[]
        };
        for arc in graph.out(u).iter().chain(incoming) {
            if arc.node <= u {
                continue;
            }
            let id = &mut best[arc.node];
            if *id == NONE {
                count += 1;
                if count > bound && !keep_dense_prefix {
                    return None;
                }
                touched.push(arc.node);
            }
            *id = (*id).min(arc.edge);
        }
        for v in touched.drain(..) {
            let id = builder
                .add_edge(u, v, std::mem::replace(&mut best[v], NONE))
                .unwrap();
            if keep_dense_prefix && id + 1 == limit {
                return Some(builder.build());
            }
        }
    }
    Some(builder.build())
}

#[derive(Clone, Copy)]
struct OrientedEdge {
    from: NodeId,
    to: NodeId,
    low: usize,
    next_low: usize,
    nesting: usize,
}

struct Orientation {
    height: Vec<usize>,
    parent: Vec<EdgeId>,
    roots: Vec<NodeId>,
    edges: Vec<OrientedEdge>,
}

impl Orientation {
    fn new<W>(graph: &Graph<W>) -> Self {
        let n = graph.node_count();
        let mut state = Self {
            height: vec![NONE; n],
            parent: vec![NONE; n],
            roots: Vec::new(),
            edges: vec![
                OrientedEdge {
                    from: NONE,
                    to: NONE,
                    low: 0,
                    next_low: 0,
                    nesting: 0,
                };
                graph.edge_count()
            ],
        };
        let mut stack = Vec::new();
        for root in 0..n {
            if state.height[root] != NONE {
                continue;
            }
            state.roots.push(root);
            state.height[root] = 0;
            stack.push((root, 0));
            while let Some((node, index)) = stack.last_mut() {
                let u = *node;
                if let Some(&arc) = graph.out(u).get(*index) {
                    *index += 1;
                    if state.edges[arc.edge].from != NONE {
                        continue;
                    }
                    let height = state.height[u];
                    let tree = state.height[arc.node] == NONE;
                    state.edges[arc.edge] = OrientedEdge {
                        from: u,
                        to: arc.node,
                        low: if tree { height } else { state.height[arc.node] },
                        next_low: height,
                        nesting: 0,
                    };
                    if tree {
                        state.height[arc.node] = height + 1;
                        state.parent[arc.node] = arc.edge;
                        stack.push((arc.node, 0));
                    } else {
                        state.finish_edge(arc.edge);
                    }
                } else {
                    let edge = state.parent[u];
                    stack.pop();
                    if edge != NONE {
                        state.finish_edge(edge);
                    }
                }
            }
        }
        state
    }

    fn finish_edge(&mut self, id: EdgeId) {
        let edge = self.edges[id];
        self.edges[id].nesting = 2 * edge.low + usize::from(edge.next_low < self.height[edge.from]);
        let parent = self.parent[edge.from];
        if parent != NONE {
            let parent = &mut self.edges[parent];
            match edge.low.cmp(&parent.low) {
                std::cmp::Ordering::Less => {
                    parent.next_low = parent.low.min(edge.next_low);
                    parent.low = edge.low;
                }
                std::cmp::Ordering::Equal => parent.next_low = parent.next_low.min(edge.next_low),
                std::cmp::Ordering::Greater => parent.next_low = parent.next_low.min(edge.low),
            }
        }
    }

    // Stable counting sort of all edges by nesting, then scatter into CSR
    // rows. No comparison sort is needed, even at a high-degree vertex.
    fn order(&self, signed: bool, left: &[bool]) -> Rows {
        let maximum = self
            .edges
            .iter()
            .map(|edge| edge.nesting)
            .max()
            .unwrap_or(0);
        let zero = if signed { maximum } else { 0 };
        let key = |id: EdgeId| {
            if signed && left[id] {
                zero - self.edges[id].nesting
            } else {
                zero + self.edges[id].nesting
            }
        };
        let mut buckets = vec![0; maximum + zero + 1];
        let mut offsets = vec![0; self.height.len() + 1];
        for (id, edge) in self.edges.iter().enumerate() {
            buckets[key(id)] += 1;
            offsets[edge.from + 1] += 1;
        }
        let mut end = 0;
        for count in &mut buckets {
            end += *count;
            *count = end;
        }
        let mut sorted = vec![0; self.edges.len()];
        for id in (0..self.edges.len()).rev() {
            let count = &mut buckets[key(id)];
            *count -= 1;
            sorted[*count] = id;
        }
        for node in 0..self.height.len() {
            offsets[node + 1] += offsets[node];
        }
        let mut cursor = offsets[..self.height.len()].to_vec();
        let mut ids = vec![0; self.edges.len()];
        for id in sorted {
            let next = &mut cursor[self.edges[id].from];
            ids[*next] = id;
            *next += 1;
        }
        Rows { offsets, ids }
    }

    fn embed(&self, graph: &Graph<EdgeId>, rows: &Rows, left_side: &[bool]) -> PlanarEmbedding {
        let n = graph.node_count();
        let m = graph.edge_count();
        // Dart 2e follows the DFS orientation and 2e+1 reverses it.
        let mut rotation = Rotation {
            first: vec![NONE; n],
            next: vec![NONE; 2 * m],
            previous: vec![NONE; 2 * m],
        };
        for u in 0..n {
            let mut predecessor = NONE;
            for &id in rows.at(u) {
                rotation.insert_after(u, 2 * id, predecessor);
                predecessor = 2 * id;
            }
        }

        let mut left = vec![NONE; n];
        let mut right = vec![NONE; n];
        let mut stack = Vec::new();
        for &root in &self.roots {
            stack.push((root, 0));
            while let Some((u, index)) = stack.last_mut() {
                let Some(&id) = rows.at(*u).get(*index) else {
                    stack.pop();
                    continue;
                };
                *index += 1;
                let v = self.edges[id].to;
                let reverse = 2 * id + 1;
                if self.parent[v] == id {
                    let first = rotation.first[v];
                    let last = if first == NONE {
                        NONE
                    } else {
                        rotation.previous[first]
                    };
                    rotation.insert_after(v, reverse, last);
                    rotation.first[v] = reverse;
                    left[*u] = 2 * id;
                    right[*u] = 2 * id;
                    stack.push((v, 0));
                } else if left_side[id] {
                    let before = rotation.previous[left[v]];
                    rotation.insert_after(v, reverse, before);
                    left[v] = reverse;
                } else {
                    rotation.insert_after(v, reverse, right[v]);
                }
            }
        }

        let mut arcs = Vec::with_capacity(2 * m);
        let mut offsets = Vec::with_capacity(n + 1);
        let mut at = vec![0; 2 * m];
        let mut darts = Vec::with_capacity(2 * m);
        for first in rotation.first {
            offsets.push(arcs.len());
            if first == NONE {
                continue;
            }
            let mut dart = first;
            loop {
                let edge = self.edges[dart / 2];
                at[dart] = arcs.len();
                darts.push(dart);
                arcs.push(Arc {
                    node: if dart % 2 == 0 { edge.to } else { edge.from },
                    edge: graph.edges()[dart / 2].weight,
                });
                dart = rotation.next[dart];
                if dart == first {
                    break;
                }
            }
        }
        offsets.push(arcs.len());
        let face_next = darts
            .into_iter()
            .map(|dart| at[rotation.previous[dart ^ 1]])
            .collect();
        PlanarEmbedding {
            offsets,
            rotation: arcs,
            face_next,
        }
    }
}

struct Rows {
    offsets: Vec<usize>,
    ids: Vec<EdgeId>,
}

impl Rows {
    fn at(&self, node: NodeId) -> &[EdgeId] {
        &self.ids[self.offsets[node]..self.offsets[node + 1]]
    }
}

struct Rotation {
    first: Vec<usize>,
    next: Vec<usize>,
    previous: Vec<usize>,
}

impl Rotation {
    fn insert_after(&mut self, node: NodeId, dart: usize, predecessor: usize) {
        if predecessor == NONE {
            debug_assert_eq!(self.first[node], NONE);
            self.first[node] = dart;
            self.next[dart] = dart;
            self.previous[dart] = dart;
            return;
        }
        let next = self.next[predecessor];
        self.next[predecessor] = dart;
        self.previous[dart] = predecessor;
        self.next[dart] = next;
        self.previous[next] = dart;
    }
}

#[derive(Clone, Copy)]
struct Interval {
    low: EdgeId,
    high: EdgeId,
}

impl Interval {
    const EMPTY: Self = Self {
        low: NONE,
        high: NONE,
    };

    fn empty(self) -> bool {
        self.low == NONE
    }
}

#[derive(Clone, Copy)]
struct Pair {
    left: Interval,
    right: Interval,
}

struct Constraints<'g> {
    orientation: &'g Orientation,
    stack: Vec<Pair>,
    base: Vec<usize>,
    lowest: Vec<EdgeId>,
    reference: Vec<EdgeId>,
    left: Vec<bool>,
}

impl<'g> Constraints<'g> {
    fn new(orientation: &'g Orientation, edges: usize) -> Self {
        Self {
            orientation,
            stack: Vec::new(),
            base: vec![0; edges],
            lowest: vec![NONE; edges],
            reference: vec![NONE; edges],
            left: vec![false; edges],
        }
    }

    fn test(&mut self, rows: &Rows) -> bool {
        let mut dfs = Vec::new();
        for &root in &self.orientation.roots {
            dfs.push((root, 0));
            while let Some((node, index)) = dfs.last_mut() {
                let u = *node;
                if let Some(&edge) = rows.at(u).get(*index) {
                    self.base[edge] = self.stack.len();
                    let first = *index == 0;
                    *index += 1;
                    let target = self.orientation.edges[edge].to;
                    if self.orientation.parent[target] == edge {
                        dfs.push((target, 0));
                    } else {
                        self.lowest[edge] = edge;
                        self.stack.push(Pair {
                            left: Interval::EMPTY,
                            right: Interval {
                                low: edge,
                                high: edge,
                            },
                        });
                        if !self.integrate(edge, first) {
                            return false;
                        }
                    }
                } else {
                    dfs.pop();
                    let parent = self.orientation.parent[u];
                    if parent != NONE {
                        self.trim(parent);
                        let (_, index) = dfs.last().unwrap();
                        if !self.integrate(parent, *index == 1) {
                            return false;
                        }
                    }
                }
            }
            debug_assert!(self.stack.is_empty());
        }
        true
    }

    fn low(&self, edge: EdgeId) -> usize {
        self.orientation.edges.get(edge).map_or(NONE, |e| e.low)
    }

    fn conflicts(&self, interval: Interval, edge: EdgeId) -> bool {
        !interval.empty() && self.low(interval.high) > self.low(edge)
    }

    fn integrate(&mut self, edge: EdgeId, first: bool) -> bool {
        let from = self.orientation.edges[edge].from;
        if self.low(edge) >= self.orientation.height[from] {
            return true;
        }
        let parent = self.orientation.parent[from];
        if first {
            self.lowest[parent] = self.lowest[edge];
            return true;
        }

        let mut merged = Pair {
            left: Interval::EMPTY,
            right: Interval::EMPTY,
        };
        while self.stack.len() > self.base[edge] {
            let mut pair = self.stack.pop().unwrap();
            if !pair.left.empty() {
                std::mem::swap(&mut pair.left, &mut pair.right);
            }
            if !pair.left.empty() {
                return false;
            }

            if self.low(pair.right.low) > self.low(parent) {
                self.append(&mut merged.right, pair.right);
            } else {
                self.reference[pair.right.low] = self.lowest[parent];
            }
        }

        while let Some(&last) = self.stack.last() {
            if !self.conflicts(last.left, edge) && !self.conflicts(last.right, edge) {
                break;
            }
            let mut pair = self.stack.pop().unwrap();
            if self.conflicts(pair.right, edge) {
                std::mem::swap(&mut pair.left, &mut pair.right);
            }
            if self.conflicts(pair.right, edge) {
                return false;
            }
            self.append(&mut merged.right, pair.right);
            self.append(&mut merged.left, pair.left);
        }
        if !merged.left.empty() || !merged.right.empty() {
            self.stack.push(merged);
        }
        true
    }

    fn append(&mut self, interval: &mut Interval, next: Interval) {
        if next.empty() {
            return;
        }
        if interval.empty() {
            *interval = next;
        } else {
            self.reference[interval.low] = next.high;
            interval.low = next.low;
        }
    }

    fn trim(&mut self, incoming: EdgeId) {
        let parent = self.orientation.edges[incoming].from;
        let height = self.orientation.height[parent];
        while self
            .stack
            .last()
            .is_some_and(|pair| self.low(pair.left.low).min(self.low(pair.right.low)) == height)
        {
            let pair = self.stack.pop().unwrap();
            if !pair.left.empty() {
                self.left[pair.left.low] = true;
            }
        }

        if let Some(mut pair) = self.stack.pop() {
            self.trim_interval(&mut pair.left, pair.right.low, parent);
            self.trim_interval(&mut pair.right, pair.left.low, parent);
            self.stack.push(pair);
        }
        if self.low(incoming) < height {
            let pair = self.stack.last().unwrap();
            let (l, r) = (pair.left.high, pair.right.high);
            self.reference[incoming] = if l != NONE && (r == NONE || self.low(l) > self.low(r)) {
                l
            } else {
                r
            };
        }
    }

    fn trim_interval(&mut self, interval: &mut Interval, opposite: EdgeId, parent: NodeId) {
        while interval.high != NONE && self.orientation.edges[interval.high].to == parent {
            interval.high = self.reference[interval.high];
        }
        if interval.high == NONE && !interval.empty() {
            self.reference[interval.low] = opposite;
            self.left[interval.low] = true;
            interval.low = NONE;
        }
    }

    fn resolve_sides(&mut self) {
        let mut chain = Vec::new();
        for start in 0..self.reference.len() {
            let mut edge = start;
            while self.reference[edge] != NONE {
                chain.push(edge);
                edge = std::mem::replace(&mut self.reference[edge], NONE);
            }
            let mut left = self.left[edge];
            while let Some(id) = chain.pop() {
                self.left[id] ^= left;
                left = self.left[id];
            }
        }
    }
}

impl PlanarEmbedding {
    /// Plane dual of the simple projection. One boundary per connected
    /// component is chosen as the shared exterior (face 0). A bridge becomes
    /// a loop. Each dual edge carries the original primal edge ID.
    pub fn dual_graph(&self) -> Graph<EdgeId> {
        let (count, face) = self.dart_face_ids();
        let mut first = std::collections::HashMap::new();
        let mut dual = Graph::builder(count, GraphKind::Undirected);
        for (arc, id) in self.rotation.iter().zip(face) {
            if let Some(other) = first.remove(&arc.edge) {
                dual.add_edge(other, id, arc.edge).unwrap();
            } else {
                first.insert(arc.edge, id);
            }
        }
        dual.build()
    }

    /// Returns dual face IDs in the same flattened vertex order as the
    /// clockwise neighbor lists. Face zero joins the component exteriors.
    pub fn dart_face_ids(&self) -> (usize, Vec<usize>) {
        let mut components = crate::connectivity::DisjointSet::new(self.node_count());
        for u in 0..self.node_count() {
            for a in &self.rotation[self.offsets[u]..self.offsets[u + 1]] {
                components.union(u, a.node);
            }
        }
        let mut exterior = vec![false; self.node_count()];
        let mut face = vec![NONE; self.rotation.len()];
        let mut count = 1;
        for dart in 0..face.len() {
            if face[dart] != NONE {
                continue;
            }
            let component = components.find(self.rotation[dart].node);
            let id = if std::mem::replace(&mut exterior[component], true) {
                count += 1;
                count - 1
            } else {
                0
            };
            let mut current = dart;
            while face[current] == NONE {
                face[current] = id;
                current = self.face_next[current];
            }
        }
        (count, face)
    }
}
