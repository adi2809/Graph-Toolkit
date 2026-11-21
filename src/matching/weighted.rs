//! Edmonds' primal-dual algorithm on the quotient by the current odd sets.
//! Each ring stores only its cyclic edges. Internal matches are materialized
//! when the ring opens, so augmentation touches only the quotient forest.

use super::{GeneralMatching, VACANT, maximum_general_matching_gabow};
use crate::connectivity::DisjointSet;
use crate::queue::IndexedHeap;
use crate::{EdgeId, Graph, GraphError, GraphKind, NONE};

#[derive(Clone, Debug)]
pub struct WeightedMatching {
    pub matching: GeneralMatching,
    pub total_weight: i128,
}

/// Maximizes total weight, optionally prioritizing cardinality. Pure weight
/// optimization excludes nonpositive edges. Returned IDs belong to the input.
pub fn maximum_weight_matching(
    graph: &Graph<i64>,
    max_cardinality: bool,
) -> Result<WeightedMatching, GraphError> {
    optimize(graph, max_cardinality, false)
}

/// Minimizes total weight among maximum-cardinality matchings.
pub fn minimum_cost_matching(graph: &Graph<i64>) -> Result<WeightedMatching, GraphError> {
    optimize(graph, true, true)
}

pub fn minimum_cost_perfect_matching(
    graph: &Graph<i64>,
) -> Result<Option<WeightedMatching>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    if graph.node_count() % 2 == 1 {
        return Ok(None);
    }
    let result = minimum_cost_matching(graph)?;
    Ok((result.matching.len() * 2 == graph.node_count()).then_some(result))
}

fn optimize(
    graph: &Graph<i64>,
    cardinality: bool,
    minimize: bool,
) -> Result<WeightedMatching, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    let eligible = |e: &crate::Edge<i64>| e.from != e.to && (cardinality || e.weight > 0);
    let mut sets = DisjointSet::new(n);
    let mut groups = vec![Vec::new(); n];
    // Keep a single optimal representative of each unordered pair. The first
    // edge wins ties, independently of hash iteration order.
    let mut representatives = Vec::<EdgeId>::new();
    let mut positions = std::collections::HashMap::new();
    for (id, e) in graph
        .edges()
        .iter()
        .enumerate()
        .filter(|(_, e)| eligible(e))
    {
        let pair = (e.from.min(e.to), e.from.max(e.to));
        match positions.entry(pair) {
            std::collections::hash_map::Entry::Vacant(slot) => {
                slot.insert(representatives.len());
                representatives.push(id);
                sets.union(e.from, e.to);
            }
            std::collections::hash_map::Entry::Occupied(slot) => {
                let best = &mut representatives[*slot.get()];
                let old = graph.edges()[*best].weight;
                if (minimize && e.weight < old) || (!minimize && e.weight > old) {
                    *best = id;
                }
            }
        }
    }
    drop(positions);
    for id in representatives {
        groups[sets.find(graph.edges()[id].from)].push(id);
    }
    let mut local = vec![NONE; n];
    let mut result = WeightedMatching {
        matching: GeneralMatching {
            mates: vec![VACANT; n],
            len: 0,
        },
        total_weight: 0,
    };
    for ids in groups.into_iter().filter(|g| !g.is_empty()) {
        let mut count = 0;
        let mut arcs = Vec::with_capacity(ids.len());
        for &id in &ids {
            let e = graph.edges()[id];
            for v in [e.from, e.to] {
                if local[v] == NONE {
                    local[v] = count;
                    count += 1;
                }
            }
            let w = i128::from(e.weight);
            arcs.push(Link {
                u: local[e.from],
                v: local[e.to],
                profit: if minimize { -w } else { w },
            });
        }
        let selected = if arcs.iter().all(|a| a.profit == arcs[0].profit) {
            let mut g = Graph::builder(count, GraphKind::Undirected);
            for a in &arcs {
                g.add_edge(a.u, a.v, ())?;
            }
            maximum_general_matching_gabow(&g.build())?
                .edges()
                .collect()
        } else {
            Quotient::new(count, &arcs).solve(cardinality)
        };
        for i in selected {
            let id = ids[i];
            let e = graph.edges()[id];
            result.matching.pair(e.from, e.to, id);
            result.matching.len += 1;
            result.total_weight += i128::from(e.weight);
        }
    }
    Ok(result)
}

#[derive(Clone, Copy)]
struct Link {
    u: usize,
    v: usize,
    profit: i128,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Forest {
    Outside,
    Even,
    Odd(EdgeId),
}

struct Part {
    live: bool, // currently a vertex of the quotient
    forest: Forest,
    mate: Option<EdgeId>,
    // (child, edge leading to the next child), starting at the original base.
    ring: Vec<(usize, EdgeId)>,
    price: i128,
}

impl Part {
    fn singleton() -> Self {
        Self {
            live: true,
            forest: Forest::Outside,
            mate: None,
            ring: vec![],
            price: 0,
        }
    }
}

enum Event {
    End,
    Edge(usize, EdgeId),
    Open(usize),
}

struct Quotient<'a> {
    edges: &'a [Link],
    owner: Vec<usize>,
    potential: Vec<i128>, // doubled vertex dual
    parts: Vec<Part>,
    recycled: Vec<usize>,
    incident: Vec<Vec<EdgeId>>,
    // Predicted absolute time at which an inter-part edge becomes tight.
    // Obsolete entries are checked against their current quotient endpoints.
    agenda: IndexedHeap<(i128, bool)>,
    clock: i128,
}

impl<'a> Quotient<'a> {
    fn new(n: usize, edges: &'a [Link]) -> Self {
        let initial = edges.iter().map(|e| e.profit).max().unwrap().max(0);
        let mut incident = vec![Vec::new(); n];
        for (id, e) in edges.iter().enumerate() {
            incident[e.u].push(id);
            incident[e.v].push(id);
        }
        Self {
            edges,
            incident,

            agenda: IndexedHeap::new(edges.len()),
            clock: 0,
            owner: (0..n).collect(),
            potential: vec![initial; n],
            parts: (0..n).map(|_| Part::singleton()).collect(),
            recycled: vec![],
        }
    }

    fn other(&self, id: EdgeId, part: usize) -> usize {
        let e = self.edges[id];
        if self.owner[e.u] == part {
            self.owner[e.v]
        } else {
            self.owner[e.u]
        }
    }

    fn solve(mut self, cardinality: bool) -> Vec<EdgeId> {
        'stages: loop {
            // Disjoint tight edges between exposed parts can be taken together.
            for (id, e) in self.edges.iter().enumerate() {
                let (a, b) = (self.owner[e.u], self.owner[e.v]);
                if a != b
                    && self.parts[a].mate.is_none()
                    && self.parts[b].mate.is_none()
                    && self.potential[e.u] + self.potential[e.v] == 2 * e.profit
                {
                    self.parts[a].mate = Some(id);
                    self.parts[b].mate = Some(id);
                }
            }
            for p in self.parts.iter_mut().filter(|p| p.live) {
                p.forest = if p.mate.is_none() {
                    Forest::Even
                } else {
                    Forest::Outside
                };
            }
            if !self.parts.iter().any(|p| p.live && p.mate.is_none()) {
                break;
            }
            self.clock = 0;
            let initial = (0..self.edges.len())
                .filter_map(|id| self.prediction(id).map(|key| (key, id)))
                .collect();
            self.agenda.rebuild(initial);
            loop {
                let (step, event) = self.next_event(cardinality);
                if matches!(event, Event::End) {
                    break 'stages;
                }
                if step != 0 {
                    self.advance(step);
                }
                match event {
                    Event::Open(p) => self.open(p, true),
                    Event::Edge(a, id) => {
                        let b = self.other(id, a);
                        if self.parts[b].forest == Forest::Outside {
                            self.parts[b].forest = Forest::Odd(id);
                            let partner = self.other(self.parts[b].mate.unwrap(), b);
                            self.parts[partner].forest = Forest::Even;
                            self.reschedule(partner);
                        } else if self.join(a, b, id) {
                            break;
                        }
                    }
                    Event::End => unreachable!(),
                }
            }
        }
        let mut pending: Vec<_> = (self.owner.len()..self.parts.len())
            .filter(|&p| self.parts[p].live)
            .collect();
        while let Some(p) = pending.pop() {
            pending.extend(
                self.parts[p]
                    .ring
                    .iter()
                    .map(|x| x.0)
                    .filter(|&c| c >= self.owner.len()),
            );
            self.open(p, false);
        }
        (0..self.owner.len())
            .filter_map(|v| self.parts[v].mate.filter(|&id| self.edges[id].u == v))
            .collect()
    }

    // Only edges between quotient vertices participate. Descendant odd-set
    // prices cancel from their reduced cost. Ties can be resolved arbitrarily.
    fn edge_event(&self, id: EdgeId) -> Option<(i128, usize)> {
        let e = self.edges[id];
        let (mut a, mut b) = (self.owner[e.u], self.owner[e.v]);
        if a == b {
            return None;
        }
        if self.parts[a].forest != Forest::Even {
            std::mem::swap(&mut a, &mut b);
        }
        if self.parts[a].forest != Forest::Even {
            return None;
        }
        let slack = self.potential[e.u] + self.potential[e.v] - 2 * e.profit;
        let delay = match self.parts[b].forest {
            Forest::Odd(_) => return None,
            Forest::Outside => slack,
            Forest::Even => {
                debug_assert_eq!(slack % 2, 0);
                slack / 2
            }
        };
        Some((self.clock + delay, a))
    }

    fn prediction(&self, id: EdgeId) -> Option<(i128, bool)> {
        let (time, a) = self.edge_event(id)?;
        let b = self.other(id, a);
        // An edge directly joining two exposed roots immediately augments.
        let needs_search = self.parts[a].mate.is_some() || self.parts[b].mate.is_some();
        Some((time, needs_search))
    }

    fn reschedule(&mut self, part: usize) {
        let mut stack = vec![part];
        while let Some(p) = stack.pop() {
            if p >= self.owner.len() {
                stack.extend(self.parts[p].ring.iter().map(|x| x.0));
                continue;
            }
            for &id in &self.incident[p] {
                if let Some(event) = self.prediction(id) {
                    self.agenda.offer(id, event);
                }
            }
        }
    }

    fn next_event(&mut self, cardinality: bool) -> (i128, Event) {
        let mut delta = i128::MAX;
        if !cardinality {
            for (v, &p) in self.owner.iter().enumerate() {
                if self.parts[p].forest == Forest::Even {
                    delta = delta.min(self.potential[v]);
                }
            }
        }
        let mut event = Event::End;
        for (p, part) in self.parts.iter().enumerate().skip(self.owner.len()) {
            if part.live && matches!(part.forest, Forest::Odd(_)) && part.price < delta {
                delta = part.price;
                event = Event::Open(p);
            }
        }
        while let Some(((time, _), id)) = self.agenda.peek() {
            if time - self.clock >= delta {
                break;
            }
            self.agenda.pop();
            if let Some((current, a)) = self.edge_event(id) {
                if current == time {
                    return (time - self.clock, Event::Edge(a, id));
                }
                self.agenda.decrease(id, self.prediction(id).unwrap());
            }
        }
        (delta, event)
    }

    fn advance(&mut self, delta: i128) {
        self.clock += delta;
        for (y, &p) in self.potential.iter_mut().zip(&self.owner) {
            match self.parts[p].forest {
                Forest::Even => *y -= delta,
                Forest::Odd(_) => *y += delta,
                _ => {}
            }
        }
        for p in self
            .parts
            .iter_mut()
            .skip(self.owner.len())
            .filter(|p| p.live)
        {
            match p.forest {
                Forest::Even => p.price += delta,
                Forest::Odd(_) => p.price -= delta,
                _ => {}
            }
        }
    }

    // Trace the alternating path to its exposed root in the quotient.
    fn ancestors(&self, mut v: usize) -> Vec<usize> {
        let mut path = vec![v];
        while let Some(id) = self.parts[v].mate {
            let odd = self.other(id, v);
            let Forest::Odd(incoming) = self.parts[odd].forest else {
                unreachable!()
            };
            v = self.other(incoming, odd);
            path.extend([odd, v]);
        }
        path
    }

    fn forest_edge(&self, child: usize) -> EdgeId {
        match self.parts[child].forest {
            Forest::Odd(e) => e,
            _ => self.parts[child].mate.unwrap(),
        }
    }

    fn join(&mut self, a: usize, b: usize, bridge: EdgeId) -> bool {
        let mut left = self.ancestors(a);
        let mut right = self.ancestors(b);
        if left.last() != right.last() {
            for path in [&left, &right] {
                self.parts[path[0]].mate = Some(bridge);
                for pair in path[1..].chunks_exact(2) {
                    let edge = self.forest_edge(pair[0]);
                    for &p in pair {
                        self.parts[p].mate = Some(edge);
                    }
                }
            }
            return true;
        }
        let mut base = NONE;
        while left.last().is_some() && left.last() == right.last() {
            base = left.pop().unwrap();
            right.pop();
        }
        let mut ring = Vec::new();
        let mut previous = base;
        for &p in left.iter().rev() {
            ring.push((previous, self.forest_edge(p)));
            previous = p;
        }
        ring.push((previous, bridge));
        for &p in &right {
            ring.push((p, self.forest_edge(p)));
        }
        let id = self.recycled.pop().unwrap_or(self.parts.len());
        // Reassign only leaves of the contracted subtrees.
        for &(child, _) in &ring {
            self.parts[child].live = false;
            self.assign(child, id);
        }
        let changed: Vec<_> = ring
            .iter()
            .map(|x| x.0)
            .filter(|&c| self.parts[c].forest != Forest::Even)
            .collect();
        let part = Part {
            live: true,
            forest: Forest::Even,
            mate: self.parts[base].mate,
            ring,
            price: 0,
        };
        if id == self.parts.len() {
            self.parts.push(part);
        } else {
            self.parts[id] = part;
        }
        for child in changed {
            self.reschedule(child);
        }
        false
    }

    fn assign(&mut self, child: usize, owner: usize) {
        let mut stack = vec![child];
        while let Some(p) = stack.pop() {
            if p < self.owner.len() {
                self.owner[p] = owner;
            } else {
                stack.extend(self.parts[p].ring.iter().map(|x| x.0));
            }
        }
    }

    fn open(&mut self, p: usize, restore_forest: bool) {
        let matched = self.parts[p].mate;
        let incoming = match self.parts[p].forest {
            Forest::Odd(e) => Some(e),
            _ => None,
        };
        let interior = |id: usize| {
            let e = self.edges[id];
            if self.owner[e.u] == p { e.u } else { e.v }
        };
        let exit_vertex = matched.map(interior);
        let entry_vertex = incoming.map(interior);
        let mut ring = std::mem::take(&mut self.parts[p].ring);
        for &(child, _) in &ring {
            self.assign(child, child);
            self.parts[child].live = true;
            self.parts[child].forest = Forest::Outside;
        }
        if let Some(v) = exit_vertex {
            let position = ring.iter().position(|x| x.0 == self.owner[v]).unwrap();
            ring.rotate_left(position);
        }
        self.parts[ring[0].0].mate = matched;
        for pair in ring[1..].chunks_exact(2) {
            for &(child, _) in pair {
                self.parts[child].mate = Some(pair[0].1);
            }
        }
        if restore_forest {
            // Enter at the child hit by the forest edge, leave through the
            // base. Exactly one of the two arcs of an odd ring has even length.
            let v = entry_vertex.unwrap();
            let mut i = ring.iter().position(|x| x.0 == self.owner[v]).unwrap();
            let forward = i % 2 == 1;
            let mut label = Forest::Odd(incoming.unwrap());
            loop {
                self.parts[ring[i].0].forest = label;
                if i == 0 {
                    break;
                }
                let next = if forward { (i + 1) % ring.len() } else { i - 1 };
                let edge = ring[if forward { i } else { next }].1;
                label = match label {
                    Forest::Odd(_) => Forest::Even,
                    _ => Forest::Odd(edge),
                };
                i = next;
            }
        }
        self.parts[p].live = false;
        self.recycled.push(p);
        if restore_forest {
            for (child, _) in ring {
                if !matches!(self.parts[child].forest, Forest::Odd(_)) {
                    self.reschedule(child);
                }
            }
        }
    }
}
