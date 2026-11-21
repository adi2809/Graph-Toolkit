//! Exact VF2++ whole-graph and induced subgraph isomorphism for multigraphs.

use std::cmp::Reverse;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::hash::Hash;

use crate::queue::IndexedHeap;
use crate::random::mix64;
use crate::{Graph, NONE, NodeId};

/// Finds the first vertex bijection, stored in left-vertex order. Edge
/// direction, multiplicity, and self-loops are preserved. Edge payloads
/// have no effect. Two empty graphs of the same kind are isomorphic.
pub fn vf2pp_isomorphism<A, B>(left: &Graph<A>, right: &Graph<B>) -> Option<Vec<NodeId>> {
    vf2pp_all_isomorphisms(left, right).next()
}

pub fn vf2pp_is_isomorphic<A, B>(left: &Graph<A>, right: &Graph<B>) -> bool {
    vf2pp_isomorphism(left, right).is_some()
}

/// Lazily enumerates vertex bijections. Preprocessing makes compact copies of
/// the two topologies, so the iterator can outlive its inputs. The search is
/// iterative. Highly symmetric graphs can have factorially many mappings.
pub fn vf2pp_all_isomorphisms<A, B>(left: &Graph<A>, right: &Graph<B>) -> Isomorphisms {
    Isomorphisms::new(left, right, false, |_, _| ())
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct IsomorphismLabelError {
    pub left_expected: usize,
    pub left_actual: usize,
    pub right_expected: usize,
    pub right_actual: usize,
}

impl fmt::Display for IsomorphismLabelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "expected {} left and {} right labels, got {} and {}",
            self.left_expected, self.right_expected, self.left_actual, self.right_actual
        )
    }
}

impl Error for IsomorphismLabelError {}

/// Finds a vertex bijection preserving both topology and the supplied labels.
pub fn vf2pp_isomorphism_with_labels<A, B, L: Eq + Hash>(
    left: &Graph<A>,
    right: &Graph<B>,
    left_labels: &[L],
    right_labels: &[L],
) -> Result<Option<Vec<NodeId>>, IsomorphismLabelError> {
    Ok(vf2pp_all_isomorphisms_with_labels(left, right, left_labels, right_labels)?.next())
}

pub fn vf2pp_all_isomorphisms_with_labels<A, B, L: Eq + Hash>(
    left: &Graph<A>,
    right: &Graph<B>,
    left_labels: &[L],
    right_labels: &[L],
) -> Result<Isomorphisms, IsomorphismLabelError> {
    validate_labels(left, right, left_labels, right_labels)?;
    Ok(Isomorphisms::new(left, right, false, |side, node| {
        if side == 0 {
            &left_labels[node]
        } else {
            &right_labels[node]
        }
    }))
}

/// Finds a vertex-induced copy of `pattern` in `host`. The returned vector
/// maps each pattern vertex to its host vertex. Edges and nonedges between
/// mapped vertices, direction, loops, and parallel multiplicities are exact.
/// Edges incident to host vertices outside the image are unrestricted.
pub fn vf2pp_subgraph_isomorphism<W, P>(
    host: &Graph<W>,
    pattern: &Graph<P>,
) -> Option<Vec<NodeId>> {
    vf2pp_all_subgraph_isomorphisms(host, pattern).next()
}

pub fn vf2pp_subgraph_is_isomorphic<W, P>(host: &Graph<W>, pattern: &Graph<P>) -> bool {
    vf2pp_subgraph_isomorphism(host, pattern).is_some()
}

/// Enumerates all pattern-to-host injections that preserve the induced
/// topology. Different injections with the same image are returned
/// separately. An empty pattern of the same graph kind has one empty mapping.
pub fn vf2pp_all_subgraph_isomorphisms<W, P>(host: &Graph<W>, pattern: &Graph<P>) -> Isomorphisms {
    Isomorphisms::new(pattern, host, true, |_, _| ())
}

pub fn vf2pp_subgraph_isomorphism_with_labels<W, P, L: Eq + Hash>(
    host: &Graph<W>,
    pattern: &Graph<P>,
    host_labels: &[L],
    pattern_labels: &[L],
) -> Result<Option<Vec<NodeId>>, IsomorphismLabelError> {
    Ok(
        vf2pp_all_subgraph_isomorphisms_with_labels(host, pattern, host_labels, pattern_labels)?
            .next(),
    )
}

pub fn vf2pp_all_subgraph_isomorphisms_with_labels<W, P, L: Eq + Hash>(
    host: &Graph<W>,
    pattern: &Graph<P>,
    host_labels: &[L],
    pattern_labels: &[L],
) -> Result<Isomorphisms, IsomorphismLabelError> {
    validate_labels(host, pattern, host_labels, pattern_labels)?;
    Ok(Isomorphisms::new(pattern, host, true, |side, node| {
        if side == 0 {
            &pattern_labels[node]
        } else {
            &host_labels[node]
        }
    }))
}

fn validate_labels<A, B, L>(
    left: &Graph<A>,
    right: &Graph<B>,
    left_labels: &[L],
    right_labels: &[L],
) -> Result<(), IsomorphismLabelError> {
    if left.node_count() != left_labels.len() || right.node_count() != right_labels.len() {
        return Err(IsomorphismLabelError {
            left_expected: left.node_count(),
            left_actual: left_labels.len(),
            right_expected: right.node_count(),
            right_actual: right_labels.len(),
        });
    }
    Ok(())
}

pub struct Isomorphisms {
    search: Option<Search>,
    remaining_empty: bool,
}

impl Isomorphisms {
    fn new<A, B, L: Eq + Hash>(
        left: &Graph<A>,
        right: &Graph<B>,
        induced: bool,
        label: impl Fn(usize, NodeId) -> L,
    ) -> Self {
        let mut result = Self {
            search: None,
            remaining_empty: false,
        };
        let n = left.node_count();
        // An equal-sized induced match is an ordinary isomorphism and can
        // use the stronger joint color-refinement prefilter.
        let induced = induced && n != right.node_count();
        let size_ok = if induced {
            n <= right.node_count() && left.edge_count() <= right.edge_count()
        } else {
            n == right.node_count() && left.edge_count() == right.edge_count()
        };
        if left.kind() != right.kind() || !size_ok {
            return result;
        }
        if n == 0 {
            result.remaining_empty = true;
            return result;
        }

        let mut keys = HashMap::new();
        let left_colors = seed_colors(left, |node| label(0, node), &mut keys, induced);
        let right_colors = seed_colors(right, |node| label(1, node), &mut keys, induced);
        let classes = keys.len();
        drop(keys);
        if !balanced([&left_colors, &right_colors], classes, induced) {
            return result;
        }
        let topology = [Topology::new(left), Topology::new(right)];
        let colors = [left_colors, right_colors];
        let refinement = if induced {
            Some((colors, classes))
        } else {
            refine(&topology, colors, classes)
        };
        if let Some((colors, classes)) = refinement {
            result.search = Some(Search::new(topology, colors, classes, induced));
        }
        result
    }
}

impl Iterator for Isomorphisms {
    type Item = Vec<NodeId>;

    fn next(&mut self) -> Option<Self::Item> {
        if std::mem::take(&mut self.remaining_empty) {
            return Some(Vec::new());
        }
        self.search.as_mut()?.next_mapping()
    }
}

impl std::iter::FusedIterator for Isomorphisms {}

type ColorKey<L> = (L, usize, usize, usize);

fn seed_colors<W, L: Eq + Hash>(
    graph: &Graph<W>,
    label: impl Fn(NodeId) -> L,
    keys: &mut HashMap<ColorKey<L>, usize>,
    induced: bool,
) -> Vec<usize> {
    graph
        .nodes()
        .map(|node| {
            let out = graph.out(node);
            let key = (
                label(node),
                if induced { 0 } else { out.len() },
                if induced { 0 } else { graph.inc(node).len() },
                out.iter().filter(|arc| arc.node == node).count(),
            );
            let id = keys.len();
            *keys.entry(key).or_insert(id)
        })
        .collect()
}

fn balanced(colors: [&[usize]; 2], classes: usize, induced: bool) -> bool {
    let mut counts = vec![0_isize; classes];
    for &color in colors[0] {
        counts[color] += 1;
    }
    for &color in colors[1] {
        counts[color] -= 1;
    }
    counts
        .into_iter()
        .all(|count| count == 0 || induced && count < 0)
}

// A few joint 1-WL refinement rounds tighten candidate colors before VF2++.
// Commutative fingerprints avoid sorting every neighborhood. A collision can
// only keep a larger class. The search still checks exact edge multiplicities.
fn refine(
    sides: &[Topology; 2],
    mut colors: [Vec<usize>; 2],
    mut classes: usize,
) -> Option<([Vec<usize>; 2], usize)> {
    let n = colors[0].len();
    let mut next = [vec![0; n], vec![0; n]];
    let mut table = HashMap::new();
    for _ in 0..6 {
        if classes == n {
            break;
        }
        table.clear();
        for side in 0..2 {
            let graph = &sides[side];
            for u in 0..n {
                let neighbors = |row: &[Neighbor]| {
                    row.iter().fold(0_u64, |sum, arc| {
                        let tag = (colors[side][arc.node] as u64)
                            .wrapping_mul(0x9e3779b97f4a7c15)
                            .wrapping_add(arc.count as u64);
                        sum.wrapping_add(mix64(tag))
                    })
                };
                let out = neighbors(graph.out.at(u));
                let incoming = graph
                    .incoming
                    .as_ref()
                    .map_or(0, |rows| neighbors(rows.at(u)));
                let id = table.len();
                next[side][u] = *table.entry((colors[side][u], out, incoming)).or_insert(id);
            }
        }
        let refined = table.len();
        if !balanced([&next[0], &next[1]], refined, false) {
            return None;
        }
        if refined == classes {
            break;
        }
        classes = refined;
        std::mem::swap(&mut colors, &mut next);
    }
    Some((colors, classes))
}

#[derive(Clone, Copy)]
struct Neighbor {
    node: NodeId,
    count: usize,
}

struct Rows {
    offsets: Vec<usize>,
    neighbors: Vec<Neighbor>,
}

impl Rows {
    fn new<W>(graph: &Graph<W>, backwards: bool) -> Self {
        let mut state = Self {
            offsets: Vec::with_capacity(graph.node_count() + 1),
            neighbors: Vec::new(),
        };
        let mut location = vec![NONE; graph.node_count()];
        for node in graph.nodes() {
            let begin = state.neighbors.len();
            state.offsets.push(begin);
            let arcs = if backwards {
                graph.inc(node)
            } else {
                graph.out(node)
            };
            for arc in arcs {
                if arc.node == node {
                    continue;
                }
                let at = &mut location[arc.node];
                if *at == NONE || *at < begin {
                    *at = state.neighbors.len();
                    state.neighbors.push(Neighbor {
                        node: arc.node,
                        count: 1,
                    });
                } else {
                    state.neighbors[*at].count += 1;
                }
            }
        }
        state.offsets.push(state.neighbors.len());
        state
    }

    fn at(&self, node: NodeId) -> &[Neighbor] {
        &self.neighbors[self.offsets[node]..self.offsets[node + 1]]
    }
}

struct Topology {
    out: Rows,
    incoming: Option<Rows>,
}

impl Topology {
    fn new<W>(graph: &Graph<W>) -> Self {
        Self {
            out: Rows::new(graph, false),
            incoming: graph.is_directed().then(|| Rows::new(graph, true)),
        }
    }

    fn inc(&self) -> &Rows {
        self.incoming.as_ref().unwrap_or(&self.out)
    }

    fn around(&self, node: NodeId) -> impl Iterator<Item = &Neighbor> {
        let incoming = self.incoming.as_ref().map_or(&[][..], |rows| rows.at(node));
        self.out.at(node).iter().chain(incoming)
    }

    fn degree(&self, node: NodeId) -> usize {
        self.out.at(node).len() + self.incoming.as_ref().map_or(0, |rows| rows.at(node).len())
    }
}

// Rare/high-degree roots, BFS layers, then the most constraints to earlier
// vertices. An indexed heap uses a single live entry per vertex.
fn matching_order(graph: &Topology, color: &[usize], counts: &[usize]) -> Vec<NodeId> {
    let n = color.len();
    let mut roots: Vec<_> = (0..n).collect();
    roots.sort_unstable_by_key(|&u| (counts[color[u]], Reverse(graph.degree(u)), u));
    let mut state = vec![0_u8; n];
    let mut prior = vec![0; n];
    let mut current = Vec::new();
    let mut next = Vec::new();
    let mut order = Vec::with_capacity(n);
    let mut heap = IndexedHeap::new(n);
    let priority = |u, score| Reverse((score, graph.degree(u), Reverse(counts[color[u]])));
    for root in roots {
        if state[root] != 0 {
            continue;
        }
        state[root] = 1;
        current.push(root);
        while !current.is_empty() {
            for u in current.drain(..) {
                state[u] = 2;
                heap.decrease(u, priority(u, prior[u]));
            }
            while let Some((_, u)) = heap.pop() {
                state[u] = 3;
                order.push(u);
                for arc in graph.around(u) {
                    let v = arc.node;
                    match state[v] {
                        0 => {
                            state[v] = 1;
                            next.push(v);
                            prior[v] += 1;
                        }
                        1 => prior[v] += 1,
                        2 => {
                            prior[v] += 1;
                            heap.decrease(v, priority(v, prior[v]));
                        }
                        _ => {}
                    }
                }
            }
            std::mem::swap(&mut current, &mut next);
        }
    }
    order
}

#[derive(Clone, Copy)]
enum Candidates {
    Color(usize),
    Leaves(NodeId),
    Enters(NodeId),
}

struct Frame {
    candidates: Candidates,
    cursor: usize,
}

struct Search {
    graph: [Topology; 2],
    colors: [Vec<usize>; 2],
    order: Vec<NodeId>,
    groups: FreeGroups,
    map: [Vec<NodeId>; 2],
    frontier: [Vec<[usize; 2]>; 2],
    frames: Vec<Frame>,
    scratch: CutScratch,
    balance: Option<FrontierBalance>,
}

impl Search {
    fn new(graph: [Topology; 2], colors: [Vec<usize>; 2], classes: usize, induced: bool) -> Self {
        let (n, m) = (colors[0].len(), colors[1].len());
        let balance = induced.then(|| FrontierBalance::new(&colors, classes));
        let groups = FreeGroups::new(&colors[1], classes);
        let order = matching_order(&graph[0], &colors[0], &groups.count);
        let mut state = Self {
            graph,
            colors,
            order,
            groups,
            map: [vec![NONE; n], vec![NONE; m]],
            frontier: [vec![[0; 2]; n], vec![[0; 2]; m]],
            frames: Vec::with_capacity(n),
            scratch: CutScratch {
                mapped_count: vec![0; m],
                groups: vec![[0; 2]; 4 * classes],
                mapped_touched: Vec::new(),
                group_touched: Vec::new(),
            },
            balance,
        };
        state.frames.push(state.frame(0));
        state
    }

    fn frame(&self, depth: usize) -> Frame {
        let u = self.order[depth];
        let color = self.colors[0][u];
        let mut candidates = Candidates::Color(color);
        let mut count = self.groups.count[color];
        if count > 1 {
            let patterns = [
                (
                    self.graph[0].out.at(u),
                    self.graph[1].inc(),
                    Candidates::Enters as fn(_) -> _,
                ),
                (
                    self.graph[0]
                        .incoming
                        .as_ref()
                        .map_or(&[][..], |rows| rows.at(u)),
                    &self.graph[1].out,
                    Candidates::Leaves as fn(_) -> _,
                ),
            ];
            for (adjacent, row, source) in patterns {
                for neighbor in adjacent {
                    let mapped = self.map[0][neighbor.node];
                    if mapped != NONE && row.at(mapped).len() < count {
                        candidates = source(mapped);
                        count = row.at(mapped).len();
                    }
                }
            }
        }
        let cursor = match candidates {
            Candidates::Color(color) => self.groups.head[color],
            _ => 0,
        };
        Frame { candidates, cursor }
    }

    fn next_mapping(&mut self) -> Option<Vec<NodeId>> {
        while let Some(depth) = self.frames.len().checked_sub(1) {
            let frame = &mut self.frames[depth];
            let target = match frame.candidates {
                Candidates::Color(_) => (frame.cursor != NONE).then_some(frame.cursor),
                Candidates::Enters(node) => self.graph[1]
                    .inc()
                    .at(node)
                    .get(frame.cursor)
                    .map(|n| n.node),
                Candidates::Leaves(node) => {
                    self.graph[1].out.at(node).get(frame.cursor).map(|n| n.node)
                }
            };
            let Some(v) = target else {
                self.frames.pop();
                if depth > 0 {
                    let u = self.order[depth - 1];
                    self.update(u, self.map[0][u], false);
                }
                continue;
            };
            frame.cursor = if matches!(frame.candidates, Candidates::Color(_)) {
                self.groups.next[v]
            } else {
                frame.cursor + 1
            };
            let u = self.order[depth];
            if self.map[1][v] != NONE
                || self.colors[0][u] != self.colors[1][v]
                || flags(self.frontier[0][u]) != flags(self.frontier[1][v])
                || !self.feasible(u, v)
            {
                continue;
            }
            self.update(u, v, true);
            if self
                .balance
                .as_ref()
                .is_some_and(|balance| balance.deficits != 0)
            {
                self.update(u, v, false);
                continue;
            }
            if depth + 1 == self.order.len() {
                let mapping = self.map[0].clone();
                self.update(u, v, false);
                return Some(mapping);
            }
            self.frames.push(self.frame(depth + 1));
        }
        None
    }

    fn feasible(&mut self, u: NodeId, v: NodeId) -> bool {
        let extra_vertices = self.map[1].len() - self.map[0].len();
        if !self.scratch.check(
            [self.graph[0].out.at(u), self.graph[1].out.at(v)],
            &self.map,
            &self.colors,
            &self.frontier,
            extra_vertices,
        ) {
            return false;
        }
        if self.graph[0].incoming.is_some() {
            self.scratch.check(
                [self.graph[0].inc().at(u), self.graph[1].inc().at(v)],
                &self.map,
                &self.colors,
                &self.frontier,
                extra_vertices,
            )
        } else {
            true
        }
    }

    fn update(&mut self, u: NodeId, v: NodeId, add: bool) {
        self.groups.set_free(v, self.colors[1][v], !add);
        self.map[0][u] = if add { v } else { NONE };
        self.map[1][v] = if add { u } else { NONE };
        for (side, node) in [u, v].into_iter().enumerate() {
            let graph = &self.graph[side];
            let frontier = &mut self.frontier[side];
            let map = &self.map[side];
            let colors = &self.colors[side];
            let weight = if side == 0 { -1 } else { 1 };
            let balance = &mut self.balance;
            if let Some(balance) = balance {
                let delta = if add { -weight } else { weight };
                balance.change(4 * colors[node] + flags(frontier[node]), delta);
            }
            let mut change = |rows: &Rows, direction| {
                for adjacent in rows.at(node) {
                    let v = adjacent.node;
                    let before = flags(frontier[v]);
                    let count = &mut frontier[v][direction];
                    if add {
                        *count += 1;
                    } else {
                        *count -= 1;
                    }
                    if let Some(balance) = balance.as_mut().filter(|_| map[v] == NONE) {
                        let after = flags(frontier[v]);
                        if before != after {
                            balance.change(4 * colors[v] + before, -weight);
                            balance.change(4 * colors[v] + after, weight);
                        }
                    }
                }
            };
            change(&graph.out, 0);
            if let Some(incoming) = &graph.incoming {
                change(incoming, 1);
            }
        }
    }
}

// Every unmatched pattern vertex still needs a host vertex with the same
// label and adjacency-to-mapped frontier bits. Maintain these necessary
// inequalities incrementally, including vertices outside both frontiers.
struct FrontierBalance {
    surplus: Vec<isize>,
    deficits: usize,
}

impl FrontierBalance {
    fn new(colors: &[Vec<usize>; 2], classes: usize) -> Self {
        let mut state = Self {
            surplus: vec![0; 4 * classes],
            deficits: 0,
        };
        for (color, weight) in colors.iter().zip([-1, 1]) {
            for &id in color {
                state.change(4 * id, weight);
            }
        }
        state
    }

    fn change(&mut self, class: usize, delta: isize) {
        let surplus = &mut self.surplus[class];
        self.deficits -= usize::from(*surplus < 0);
        *surplus += delta;
        self.deficits += usize::from(*surplus < 0);
    }
}

fn flags(count: [usize; 2]) -> usize {
    usize::from(count[0] != 0) + 2 * usize::from(count[1] != 0)
}

// Intrusive free lists make very large color classes cheap to enumerate.
// Backtracking restores removed nodes in the reverse order of removal.
struct FreeGroups {
    head: Vec<NodeId>,
    count: Vec<usize>,
    next: Vec<NodeId>,
    previous: Vec<NodeId>,
}

impl FreeGroups {
    fn new(colors: &[usize], classes: usize) -> Self {
        let mut state = Self {
            head: vec![NONE; classes],
            count: vec![0; classes],
            next: vec![NONE; colors.len()],
            previous: vec![NONE; colors.len()],
        };
        for (node, &color) in colors.iter().enumerate().rev() {
            let following = state.head[color];
            state.next[node] = following;
            state.set_free(node, color, true);
        }
        state
    }

    fn set_free(&mut self, node: NodeId, color: usize, free: bool) {
        let (previous, next) = (self.previous[node], self.next[node]);
        let (from, to) = if free {
            self.count[color] += 1;
            (node, node)
        } else {
            self.count[color] -= 1;
            (previous, next)
        };
        if previous == NONE {
            self.head[color] = to;
        } else {
            self.next[previous] = to;
        }
        if next != NONE {
            self.previous[next] = from;
        }
    }
}

struct CutScratch {
    mapped_count: Vec<usize>,
    groups: Vec<[usize; 2]>,
    mapped_touched: Vec<NodeId>,
    group_touched: Vec<usize>,
}

impl CutScratch {
    // The VF2++ cuts compare mapped adjacencies and the number of
    // neighbors/edges in each label and frontier class of unmapped vertices.
    fn check(
        &mut self,
        rows: [&[Neighbor]; 2],
        map: &[Vec<NodeId>; 2],
        color: &[Vec<usize>; 2],
        frontier: &[Vec<[usize; 2]>; 2],
        extra_vertices: usize,
    ) -> bool {
        for node in self.mapped_touched.drain(..) {
            self.mapped_count[node] = 0;
        }
        for class in self.group_touched.drain(..) {
            self.groups[class] = [0; 2];
        }
        if rows[1].len() < rows[0].len() || rows[1].len() - rows[0].len() > extra_vertices {
            return false;
        }
        for &Neighbor { node, count } in rows[0] {
            let mapped = map[0][node];
            if mapped == NONE {
                let class = 4 * color[0][node] + flags(frontier[0][node]);
                let entry = &mut self.groups[class];
                if entry[0] == 0 {
                    self.group_touched.push(class);
                }
                entry[0] += 1;
                entry[1] += count;
            } else {
                self.mapped_touched.push(mapped);
                self.mapped_count[mapped] = count;
            }
        }
        let mut missing_mapped = self.mapped_touched.len();
        for &Neighbor { node, count } in rows[1] {
            if map[1][node] == NONE {
                let class = 4 * color[1][node] + flags(frontier[1][node]);
                let entry = &mut self.groups[class];
                if extra_vertices == 0 && (entry[0] == 0 || entry[1] < count) {
                    return false;
                }
                entry[0] = entry[0].saturating_sub(1);
                entry[1] = entry[1].saturating_sub(count);
            } else {
                if self.mapped_count[node] != count {
                    return false;
                }
                missing_mapped -= 1;
            }
        }
        missing_mapped == 0
            && self
                .group_touched
                .iter()
                .all(|&class| self.groups[class] == [0; 2])
    }
}

/// AHU subtree classification for two unordered rooted undirected trees.
/// Returns false if either input fails the tree test. Sorting child classes
/// gives O(n log n) expected time with hash interning and O(n) memory.
pub fn rooted_tree_is_isomorphic<A, B>(
    a: &Graph<A>,
    root_a: usize,
    b: &Graph<B>,
    root_b: usize,
) -> Result<bool, crate::GraphError> {
    for kind in [a.kind(), b.kind()] {
        if kind == crate::GraphKind::Directed {
            return Err(crate::GraphError::RequiresUndirected);
        }
    }
    a.validate_node(root_a)?;
    b.validate_node(root_b)?;
    if a.node_count() != b.node_count()
        || a.edge_count() != a.node_count() - 1
        || b.edge_count() != a.edge_count()
    {
        return Ok(false);
    }
    fn classify<W>(
        g: &Graph<W>,
        root: usize,
        classes: &mut HashMap<Vec<usize>, usize>,
    ) -> Option<usize> {
        let n = g.node_count();
        let mut parent = vec![NONE; n];
        parent[root] = root;
        let mut order = vec![root];
        let mut head = 0;
        while let Some(&u) = order.get(head) {
            head += 1;
            for a in g.out(u) {
                if parent[a.node] == NONE {
                    parent[a.node] = u;
                    order.push(a.node);
                }
            }
        }
        if order.len() != n {
            return None;
        }
        let mut label = vec![0; n];
        for u in order.into_iter().rev() {
            let mut children: Vec<_> = g
                .out(u)
                .iter()
                .filter(|a| parent[a.node] == u)
                .map(|a| label[a.node])
                .collect();
            children.sort_unstable();
            let next = classes.len();
            label[u] = *classes.entry(children).or_insert(next);
        }
        Some(label[root])
    }
    let mut classes = HashMap::new();
    Ok(classify(a, root_a, &mut classes)
        .is_some_and(|label| classify(b, root_b, &mut classes) == Some(label)))
}
