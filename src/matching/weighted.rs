//! Weighted primal-dual blossom matching on dense vertex and edge IDs.
//!
//! Blossom expansion and lifting follow the reference algorithms in NetworkX
//! matching.py (Joris van Rantwijk and the NetworkX developers).
//! See LICENSE.networkx for the applicable BSD-3-Clause notice.

use super::{GeneralMatching, VACANT, maximum_general_matching_gabow};
use crate::connectivity::DisjointSet;
use crate::{EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

/// An exact integer objective value and a matching in the original graph.
#[derive(Clone, Debug)]
pub struct WeightedMatching {
    pub matching: GeneralMatching,
    pub total_weight: i128,
}

/// Maximum-weight matching in an undirected graph. Set `max_cardinality` to
/// maximize the number of edges first, then break ties by total weight.
/// Without this option the empty matching is eligible and zero-weight edges
/// are omitted. Self-loops are ignored. Parallel edges retain original IDs.
///
/// Uses a sparse primal-dual blossom implementation with O(V³ + E) worst-case
/// time after coalescing parallel edges. Components are solved independently.
pub fn maximum_weight_matching(
    graph: &Graph<i64>,
    max_cardinality: bool,
) -> Result<WeightedMatching, GraphError> {
    weighted_matching(graph, false, max_cardinality)
}

/// Minimum cost among all maximum-cardinality matchings.
pub fn minimum_cost_matching(graph: &Graph<i64>) -> Result<WeightedMatching, GraphError> {
    weighted_matching(graph, true, true)
}

/// Minimum-cost matching that covers every vertex, or `None` if infeasible.
pub fn minimum_cost_perfect_matching(
    graph: &Graph<i64>,
) -> Result<Option<WeightedMatching>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    if graph.node_count() % 2 != 0 {
        return Ok(None);
    }
    let matching = minimum_cost_matching(graph)?;
    Ok((2 * matching.matching.len() == graph.node_count()).then_some(matching))
}

#[derive(Clone, Copy)]
struct Score {
    value: i128,
    original: EdgeId,
}

fn weighted_matching(
    graph: &Graph<i64>,
    minimize: bool,
    max_cardinality: bool,
) -> Result<WeightedMatching, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    let mut sets = DisjointSet::new(n);
    let mut best = vec![NONE; n];
    let mut neighbors = Vec::new();
    let mut selected = Vec::new();
    for node in 0..n {
        for arc in graph.out(node) {
            if arc.node <= node || (!max_cardinality && graph.edges()[arc.edge].weight <= 0) {
                continue;
            }
            let previous = best[arc.node];
            if previous == NONE {
                neighbors.push(arc.node);
                best[arc.node] = arc.edge;
            } else {
                let old_weight = graph.edges()[previous].weight;
                let new_weight = graph.edges()[arc.edge].weight;
                if (minimize && new_weight < old_weight) || (!minimize && new_weight > old_weight) {
                    best[arc.node] = arc.edge;
                }
            }
        }
        for neighbor in neighbors.drain(..) {
            sets.union(node, neighbor);
            selected.push(std::mem::replace(&mut best[neighbor], NONE));
        }
    }

    let mut components = vec![Vec::new(); n];
    for id in selected {
        components[sets.find(graph.edges()[id].from)].push(id);
    }
    let mut result = WeightedMatching {
        matching: GeneralMatching {
            mates: vec![VACANT; n],
            len: 0,
        },
        total_weight: 0,
    };
    let mut local = best;
    let mut vertices = Vec::new();
    for edges in components.into_iter().filter(|edges| !edges.is_empty()) {
        if edges.len() == 1 {
            add_result(graph, &mut result, edges[0]);
            continue;
        }
        vertices.clear();
        for &id in &edges {
            let edge = &graph.edges()[id];
            for node in [edge.from, edge.to] {
                if local[node] == NONE {
                    local[node] = vertices.len();
                    vertices.push(node);
                }
            }
        }
        let first = graph.edges()[edges[0]].weight;
        let mut uniform = true;
        let mut build = Graph::builder(vertices.len(), GraphKind::Undirected);
        for id in edges {
            let edge = &graph.edges()[id];
            uniform &= edge.weight == first;
            let value = i128::from(edge.weight);
            let score = Score {
                original: id,
                value: if minimize { -value } else { value },
            };
            build
                .add_edge(local[edge.from], local[edge.to], score)
                .unwrap();
        }
        let component = build.build();
        if uniform {
            let matching = maximum_general_matching_gabow(&component).unwrap();
            for id in matching.edges() {
                add_result(graph, &mut result, component.edges()[id].weight.original);
            }
        } else {
            let mates = WeightedSearch::new(&component, max_cardinality).solve();
            for (node, end) in mates.into_iter().enumerate() {
                if end != NONE && node < endpoint(&component, end ^ 1) {
                    add_result(
                        graph,
                        &mut result,
                        component.edges()[end / 2].weight.original,
                    );
                }
            }
        }
    }
    Ok(result)
}

fn add_result(graph: &Graph<i64>, result: &mut WeightedMatching, id: EdgeId) {
    let edge = &graph.edges()[id];
    result.matching.pair(edge.from, edge.to, id);
    result.matching.len += 1;
    result.total_weight += i128::from(edge.weight);
}

// End 2e names the source of edge e and 2e+1 names its target. An oriented
// edge is represented by its starting end, so XOR 1 reverses its orientation.
fn endpoint(graph: &Graph<Score>, end: usize) -> NodeId {
    let edge = &graph.edges()[end / 2];
    if end & 1 == 0 { edge.from } else { edge.to }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
enum Label {
    #[default]
    Free,
    Outer,
    Inner,
    Visited,
}

struct Block {
    parent: usize,
    base: NodeId,
    label: Label,
    entry: usize,
    best: EdgeId,
    dual: i128,
}

impl Block {
    fn new(base: NodeId, dual: i128) -> Self {
        Self {
            parent: NONE,
            base,
            label: Label::Free,
            entry: NONE,
            best: NONE,
            dual,
        }
    }
}

#[derive(Default)]
struct Cycle {
    children: Vec<usize>,
    exits: Vec<usize>,
    /// Least-slack edges to distinct neighboring outer blocks in this stage.
    cache: Vec<EdgeId>,
    cached: bool,
}

impl Cycle {
    fn step(&self, index: usize, forward: bool) -> usize {
        if forward {
            (index + 1) % self.children.len()
        } else {
            index.checked_sub(1).unwrap_or(self.children.len() - 1)
        }
    }

    fn exit(&self, index: usize, forward: bool) -> usize {
        if forward {
            self.exits[index]
        } else {
            self.exits[self.step(index, false)] ^ 1
        }
    }
}

enum Edit {
    Start(usize, NodeId),
    Walk {
        block: usize,
        cursor: usize,
        start: usize,
        vertex: NodeId,
    },
    Pair(usize),
}

#[derive(Clone, Copy)]
enum Event {
    Stop,
    Tighten(EdgeId),
    Expand(usize),
}

struct WeightedSearch<'a> {
    graph: &'a Graph<Score>,
    n: usize,
    cardinality_first: bool,
    blocks: Vec<Block>,
    cycles: Vec<Cycle>,
    top: Vec<usize>,
    mate: Vec<usize>,
    allowed: Vec<bool>,
    queue: Vec<NodeId>,
    free_blocks: Vec<usize>,
    // Work buffers for nested blossoms, contraction, and iterative lifting.
    leaves: Vec<NodeId>,
    visit: Vec<usize>,
    trail: Vec<usize>,
    expansion: Vec<usize>,
    edits: Vec<Edit>,
    best_to: Vec<EdgeId>,
    touched: Vec<usize>,
}

impl<'a> WeightedSearch<'a> {
    fn new(graph: &'a Graph<Score>, cardinality_first: bool) -> Self {
        let n = graph.node_count();
        let maximum = graph
            .edges()
            .iter()
            .map(|edge| edge.weight.value)
            .max()
            .unwrap_or(0)
            .max(0);
        let blocks = (0..n)
            .map(|node| Block::new(node, maximum))
            .chain((0..n).map(|_| Block::new(NONE, 0)))
            .collect();
        Self {
            graph,
            n,
            cardinality_first,
            blocks,
            cycles: (0..n).map(|_| Cycle::default()).collect(),
            top: (0..n).collect(),
            mate: vec![NONE; n],
            allowed: vec![false; graph.edge_count()],
            queue: Vec::new(),
            free_blocks: (n..2 * n).collect(),
            leaves: Vec::new(),
            visit: Vec::new(),
            trail: Vec::new(),
            expansion: Vec::new(),
            edits: Vec::new(),
            best_to: vec![NONE; 2 * n],
            touched: Vec::new(),
        }
    }

    fn solve(mut self) -> Vec<usize> {
        while self.augment_stage() {
            for block in self.n..2 * self.n {
                let state = &self.blocks[block];
                if state.base != NONE
                    && state.parent == NONE
                    && state.label == Label::Outer
                    && state.dual == 0
                {
                    self.expand(block, true);
                }
            }
        }
        self.mate
    }

    fn source(&self, end: usize) -> NodeId {
        endpoint(self.graph, end)
    }

    fn target(&self, end: usize) -> NodeId {
        self.source(end ^ 1)
    }

    fn slack(&self, id: EdgeId) -> i128 {
        let edge = &self.graph.edges()[id];
        self.blocks[edge.from].dual + self.blocks[edge.to].dual - 2 * edge.weight.value
    }

    fn improves(&self, previous: EdgeId, slack: i128) -> bool {
        previous == NONE || slack < self.slack(previous)
    }

    fn list_leaves(&mut self, block: usize) {
        self.leaves.clear();
        self.visit.push(block);
        while let Some(item) = self.visit.pop() {
            if item < self.n {
                self.leaves.push(item);
            } else {
                self.visit
                    .extend(self.cycles[item - self.n].children.iter().rev());
            }
        }
    }

    fn label(&mut self, vertex: NodeId, label: Label, entry: usize) {
        let block = self.top[vertex];
        for id in [vertex, block] {
            self.blocks[id].label = label;
            self.blocks[id].entry = entry;
            self.blocks[id].best = NONE;
        }
    }

    fn outer(&mut self, vertex: NodeId, entry: usize) {
        self.label(vertex, Label::Outer, entry);
        self.list_leaves(self.top[vertex]);
        self.queue.extend_from_slice(&self.leaves);
    }

    fn inner(&mut self, vertex: NodeId, entry: usize) {
        self.label(vertex, Label::Inner, entry);
        let base = self.blocks[self.top[vertex]].base;
        let matched = self.mate[base];
        self.outer(self.target(matched), matched);
    }

    fn augment_stage(&mut self) -> bool {
        for block in &mut self.blocks {
            block.label = Label::Free;
            block.entry = NONE;
            block.best = NONE;
        }
        for cycle in &mut self.cycles {
            cycle.cache.clear();
            cycle.cached = false;
        }
        self.allowed.fill(false);
        self.queue.clear();
        for node in 0..self.n {
            if self.mate[node] == NONE && self.blocks[self.top[node]].label == Label::Free {
                self.outer(node, NONE);
            }
        }
        // Every augmenting search alternates tight-edge scans and dual updates.
        loop {
            while let Some(node) = self.queue.pop() {
                if self.scan(node) {
                    return true;
                }
            }
            let Some((delta, event)) = self.next_event() else {
                return false;
            };
            if matches!(event, Event::Stop) {
                return false;
            }
            if delta != 0 {
                self.adjust_duals(delta);
            }
            match event {
                Event::Tighten(id) => {
                    self.allowed[id] = true;
                    let edge = &self.graph.edges()[id];
                    let node = if self.blocks[self.top[edge.from]].label == Label::Outer {
                        edge.from
                    } else {
                        edge.to
                    };
                    self.queue.push(node);
                }
                Event::Expand(block) => self.expand(block, false),
                Event::Stop => unreachable!(),
            }
        }
    }

    fn scan(&mut self, node: NodeId) -> bool {
        for &arc in self.graph.out(node) {
            let other = arc.node;
            let left = self.top[node];
            let right = self.top[other];
            if left == right {
                continue;
            }
            let reduced = if self.allowed[arc.edge] {
                0
            } else {
                self.slack(arc.edge)
            };
            if reduced <= 0 {
                self.allowed[arc.edge] = true;
                let edge = &self.graph.edges()[arc.edge];
                let end = 2 * arc.edge + usize::from(edge.to == node);
                match self.blocks[right].label {
                    Label::Free => self.inner(other, end),
                    Label::Outer => {
                        if let Some(base) = self.common_base(node, other) {
                            self.contract(base, end);
                        } else {
                            self.augment(end);
                            return true;
                        }
                    }
                    Label::Inner => {
                        if self.blocks[other].label == Label::Free {
                            self.blocks[other].label = Label::Inner;
                            self.blocks[other].entry = end;
                        }
                    }
                    Label::Visited => unreachable!(),
                }
            } else if self.blocks[right].label == Label::Outer {
                if self.improves(self.blocks[left].best, reduced) {
                    self.blocks[left].best = arc.edge;
                }
            } else if self.blocks[other].label == Label::Free
                && self.improves(self.blocks[other].best, reduced)
            {
                self.blocks[other].best = arc.edge;
            }
        }
        false
    }

    fn common_base(&mut self, mut a: NodeId, mut b: NodeId) -> Option<NodeId> {
        let mut base = None;
        while a != NONE {
            let block = self.top[a];
            if self.blocks[block].label == Label::Visited {
                base = Some(self.blocks[block].base);
                break;
            }
            self.blocks[block].label = Label::Visited;
            self.trail.push(block);
            let previous = self.blocks[block].entry;
            a = if previous == NONE {
                NONE
            } else {
                let inner = self.top[self.source(previous)];
                self.source(self.blocks[inner].entry)
            };
            if b != NONE {
                std::mem::swap(&mut a, &mut b);
            }
        }
        for item in self.trail.drain(..) {
            self.blocks[item].label = Label::Outer;
        }
        base
    }

    fn contract(&mut self, base: NodeId, join: usize) {
        let block = self.free_blocks.pop().unwrap();
        let base_child = self.top[base];
        let mut left = self.top[self.source(join)];
        let mut right = self.top[self.target(join)];
        let mut cycle = std::mem::take(&mut self.cycles[block - self.n]);
        while left != base_child {
            self.blocks[left].parent = block;
            cycle.children.push(left);
            let entry = self.blocks[left].entry;
            cycle.exits.push(entry);
            left = self.top[self.source(entry)];
        }
        self.blocks[base_child].parent = block;
        cycle.children.push(base_child);
        cycle.children.reverse();
        cycle.exits.reverse();
        cycle.exits.push(join);
        while right != base_child {
            self.blocks[right].parent = block;
            cycle.children.push(right);
            let entry = self.blocks[right].entry;
            cycle.exits.push(entry ^ 1);
            right = self.top[self.source(entry)];
        }
        self.blocks[block] = Block::new(base, 0);
        self.blocks[block].label = Label::Outer;
        self.blocks[block].entry = self.blocks[base_child].entry;
        self.cycles[block - self.n] = cycle;
        self.list_leaves(block);
        for &node in &self.leaves {
            if self.blocks[self.top[node]].label == Label::Inner {
                self.queue.push(node);
            }
            self.top[node] = block;
        }
        self.refresh_boundary(block);
    }

    fn refresh_boundary(&mut self, block: usize) {
        for index in 0..self.cycles[block - self.n].children.len() {
            let child = self.cycles[block - self.n].children[index];
            self.blocks[child].best = NONE;
            if child >= self.n && self.cycles[child - self.n].cached {
                let cache = std::mem::take(&mut self.cycles[child - self.n].cache);
                self.cycles[child - self.n].cached = false;
                for id in cache {
                    self.boundary_edge(block, id);
                }
            } else {
                self.list_leaves(child);
                for leaf in 0..self.leaves.len() {
                    for arc in self.graph.out(self.leaves[leaf]) {
                        self.boundary_edge(block, arc.edge);
                    }
                }
            }
        }
        let mut best = NONE;
        let mut slack = i128::MAX;
        let cache = &mut self.cycles[block - self.n];
        for neighbor in self.touched.drain(..) {
            cache
                .cache
                .push(std::mem::replace(&mut self.best_to[neighbor], NONE));
        }
        cache.cached = true;
        for &id in &self.cycles[block - self.n].cache {
            let cost = self.slack(id);
            if cost < slack {
                best = id;
                slack = cost;
            }
        }
        self.blocks[block].best = best;
    }

    fn boundary_edge(&mut self, block: usize, id: EdgeId) {
        let edge = &self.graph.edges()[id];
        let a = self.top[edge.from];
        let b = self.top[edge.to];
        let neighbor = if b == block { a } else { b };
        if neighbor == block || self.blocks[neighbor].label != Label::Outer {
            return;
        }
        let previous = self.best_to[neighbor];
        if self.improves(previous, self.slack(id)) {
            if previous == NONE {
                self.touched.push(neighbor);
            }
            self.best_to[neighbor] = id;
        }
    }

    fn next_event(&self) -> Option<(i128, Event)> {
        let mut best = None;
        if !self.cardinality_first {
            let minimum = self.blocks[..self.n]
                .iter()
                .map(|state| state.dual)
                .min()
                .unwrap();
            best = Some((minimum, Event::Stop));
        }
        for node in 0..self.n {
            let id = self.blocks[node].best;
            if id != NONE && self.blocks[self.top[node]].label == Label::Free {
                earlier(&mut best, self.slack(id), Event::Tighten(id));
            }
        }
        for (block, state) in self.blocks.iter().enumerate() {
            if state.parent != NONE || state.base == NONE {
                continue;
            }
            if state.label == Label::Outer && state.best != NONE {
                let slack = self.slack(state.best);
                debug_assert_eq!(slack % 2, 0);
                earlier(&mut best, slack / 2, Event::Tighten(state.best));
            } else if state.label == Label::Inner && block >= self.n {
                earlier(&mut best, state.dual, Event::Expand(block));
            }
        }
        best
    }

    fn adjust_duals(&mut self, delta: i128) {
        for node in 0..self.n {
            match self.blocks[self.top[node]].label {
                Label::Outer => self.blocks[node].dual -= delta,
                Label::Inner => self.blocks[node].dual += delta,
                _ => {}
            }
        }
        for state in &mut self.blocks[self.n..] {
            if state.parent == NONE && state.base != NONE {
                match state.label {
                    Label::Outer => state.dual += delta,
                    Label::Inner => state.dual -= delta,
                    _ => {}
                }
            }
        }
    }

    fn expand(&mut self, block: usize, end_stage: bool) {
        self.expansion.push(block);
        while let Some(block) = self.expansion.pop() {
            for index in 0..self.cycles[block - self.n].children.len() {
                let child = self.cycles[block - self.n].children[index];
                self.blocks[child].parent = NONE;
                if end_stage && child >= self.n && self.blocks[child].dual == 0 {
                    self.expansion.push(child);
                } else {
                    self.list_leaves(child);
                    for &node in &self.leaves {
                        self.top[node] = child;
                    }
                }
            }
            if !end_stage {
                self.relabel_expansion(block);
            }
            let cycle = &mut self.cycles[block - self.n];
            cycle.children.clear();
            cycle.exits.clear();
            cycle.cache.clear();
            cycle.cached = false;
            self.blocks[block] = Block::new(NONE, 0);
            self.free_blocks.push(block);
        }
    }

    fn relabel_expansion(&mut self, block: usize) {
        let slot = block - self.n;
        let mut entry = self.blocks[block].entry;
        let entry_child = self.top[self.target(entry)];
        let mut index = self.cycles[slot]
            .children
            .iter()
            .position(|&child| child == entry_child)
            .unwrap();
        let forward = index % 2 == 1;
        while index != 0 {
            let exit = self.cycles[slot].exit(index, forward);
            let current = self.target(entry);
            let next = self.target(exit);
            self.blocks[current].label = Label::Free;
            self.blocks[next].label = Label::Free;
            self.inner(current, entry);
            self.allowed[exit / 2] = true;
            index = self.cycles[slot].step(index, forward);
            entry = self.cycles[slot].exit(index, forward);
            self.allowed[entry / 2] = true;
            index = self.cycles[slot].step(index, forward);
        }
        // The base's mate is already outer. Label the base without crossing it again.
        self.label(self.target(entry), Label::Inner, entry);
        index = self.cycles[slot].step(index, forward);
        while self.cycles[slot].children[index] != entry_child {
            let child = self.cycles[slot].children[index];
            if self.blocks[child].label != Label::Outer {
                self.list_leaves(child);
                if let Some(&node) = self
                    .leaves
                    .iter()
                    .find(|&&node| self.blocks[node].label == Label::Inner)
                {
                    let entry = self.blocks[node].entry;
                    self.blocks[node].label = Label::Free;
                    let matched = self.target(self.mate[self.blocks[child].base]);
                    self.blocks[matched].label = Label::Free;
                    self.inner(node, entry);
                }
            }
            index = self.cycles[slot].step(index, forward);
        }
    }

    fn augment(&mut self, join: usize) {
        for mut current in [join, join ^ 1] {
            loop {
                let vertex = self.source(current);
                let outer = self.top[vertex];
                self.lift(outer, vertex);
                self.mate[vertex] = current;
                let previous = self.blocks[outer].entry;
                if previous == NONE {
                    break;
                }
                let inner = self.top[self.source(previous)];
                current = self.blocks[inner].entry;
                let target = self.target(current);
                self.lift(inner, target);
                self.mate[target] = current ^ 1;
            }
        }
    }

    fn lift(&mut self, block: usize, vertex: NodeId) {
        self.edits.push(Edit::Start(block, vertex));
        while let Some(edit) = self.edits.pop() {
            match edit {
                Edit::Start(block, vertex) => {
                    if block < self.n {
                        continue;
                    }
                    let mut child = vertex;
                    while self.blocks[child].parent != block {
                        child = self.blocks[child].parent;
                    }
                    let start = self.cycles[block - self.n]
                        .children
                        .iter()
                        .position(|&item| item == child)
                        .unwrap();
                    self.edits.push(Edit::Walk {
                        block,
                        cursor: start,
                        start,
                        vertex,
                    });
                    self.edits.push(Edit::Start(child, vertex));
                }
                Edit::Walk {
                    block,
                    cursor,
                    start,
                    vertex,
                } => {
                    let cycle = &mut self.cycles[block - self.n];
                    if cursor == 0 {
                        cycle.children.rotate_left(start);
                        cycle.exits.rotate_left(start);
                        self.blocks[block].base = vertex;
                    } else {
                        let forward = start % 2 == 1;
                        let first = cycle.step(cursor, forward);
                        let second = cycle.step(first, forward);
                        let edge = cycle.exit(first, forward);
                        let (a, b) = (cycle.children[first], cycle.children[second]);
                        self.edits.push(Edit::Walk {
                            block,
                            cursor: second,
                            start,
                            vertex,
                        });
                        self.edits.push(Edit::Pair(edge));
                        self.edits.push(Edit::Start(b, self.target(edge)));
                        self.edits.push(Edit::Start(a, self.source(edge)));
                    }
                }
                Edit::Pair(edge) => {
                    let (a, b) = (self.source(edge), self.target(edge));
                    self.mate[a] = edge;
                    self.mate[b] = edge ^ 1;
                }
            }
        }
    }
}

fn earlier(best: &mut Option<(i128, Event)>, value: i128, event: Event) {
    if best.is_none_or(|(previous, _)| value < previous) {
        *best = Some((value, event));
    }
}
