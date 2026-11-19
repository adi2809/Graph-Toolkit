use std::collections::VecDeque;

use super::{MaxFlow, Residual, residual_at, residual_neighbors, validate_terminals};
use crate::{Graph, GraphError, NONE, NodeId};

/// Edmonds–Karp maximum flow, using shortest residual augmenting paths.
/// Runs in O(V E²) time and O(V + E) space.
pub fn max_flow_edmonds_karp(
    graph: &Graph<u64>,
    source: NodeId,
    sink: NodeId,
) -> Result<MaxFlow, GraphError> {
    validate_terminals(graph, source, sink)?;
    let mut state = Network::new(graph);
    let n = graph.node_count();
    let mut parent = vec![None; n];
    let mut visited = vec![0; n];
    let mut queue = Vec::new();
    let mut path = Vec::new();
    let mut phase = 0;
    loop {
        phase += 1;
        queue.clear();
        queue.push(source);
        visited[source] = phase;
        let mut head = 0;
        while let Some(&node) = queue.get(head) {
            head += 1;
            if node == sink {
                break;
            }
            for arc in residual_neighbors(graph, node) {
                if visited[arc.to] != phase && state.capacity(arc) != 0 {
                    visited[arc.to] = phase;
                    parent[arc.to] = Some((node, arc));
                    queue.push(arc.to);
                }
            }
        }
        if visited[sink] != phase {
            break;
        }
        let mut node = sink;
        path.clear();
        while node != source {
            let (from, arc) = parent[node].unwrap();
            path.push(arc);
            node = from;
        }
        state.augment(&path);
    }
    Ok(state.finish(source))
}

/// Dinic maximum flow with level graphs and iterative blocking-flow search.
/// Runs in O(V² E) time and O(V + E) space.
pub fn max_flow_dinic(
    graph: &Graph<u64>,
    source: NodeId,
    sink: NodeId,
) -> Result<MaxFlow, GraphError> {
    validate_terminals(graph, source, sink)?;
    let mut state = Network::new(graph);
    let mut level = vec![NONE; graph.node_count()];
    let mut cursor = vec![0; level.len()];
    let mut queue = Vec::new();
    let mut path: Vec<(NodeId, Residual)> = Vec::new();
    loop {
        state.levels(source, &mut level, &mut queue);
        if level[sink] == NONE {
            return Ok(state.finish(source));
        }
        cursor.fill(0);
        let mut node = source;
        loop {
            if node == sink {
                let amount = path
                    .iter()
                    .map(|&(_, arc)| state.capacity(arc))
                    .min()
                    .unwrap();
                let mut first_saturated = None;
                for (position, &(parent, arc)) in path.iter().enumerate() {
                    state.push(arc, amount);
                    if first_saturated.is_none() && state.capacity(arc) == 0 {
                        first_saturated = Some((position, parent));
                    }
                }
                state.value += u128::from(amount);
                let (position, parent) = first_saturated.unwrap();
                node = parent;
                cursor[node] += 1;
                path.truncate(position);
                continue;
            }

            let Some(arc) = residual_at(graph, node, cursor[node]) else {
                level[node] = NONE;
                let Some((parent, _)) = path.pop() else { break };
                node = parent;
                cursor[node] += 1;
                continue;
            };
            if level[arc.to] == level[node] + 1 && state.capacity(arc) != 0 {
                path.push((node, arc));
                node = arc.to;
            } else {
                cursor[node] += 1;
            }
        }
    }
}

/// Boykov–Kolmogorov maximum flow. Reuses source/sink search trees across
/// augmentations. Designed for networks with short cuts, including grids.
pub fn max_flow_boykov_kolmogorov(
    graph: &Graph<u64>,
    source: NodeId,
    sink: NodeId,
) -> Result<MaxFlow, GraphError> {
    validate_terminals(graph, source, sink)?;
    let n = graph.node_count();
    let mut search = Trees {
        residual: Network::new(graph),
        source,
        sink,
        side: vec![0; n],
        parent: vec![None; n],
        queued: vec![false; n],
        active: VecDeque::new(),
        orphan: VecDeque::new(),
        timestamp: vec![0; n],
        epoch: 0,
        distance: vec![0; n],
        ancestors: Vec::new(),
    };
    search.side[source] = 1;
    search.side[sink] = 2;
    search.activate(source);
    search.activate(sink);
    while let Some((tail, crossing)) = search.grow() {
        search.augment(tail, crossing);
        search.adopt();
    }
    Ok(search.residual.finish(source))
}

struct Network<'g> {
    graph: &'g Graph<u64>,
    capacity: Vec<[u64; 2]>,
    value: u128,
}

impl<'g> Network<'g> {
    fn new(graph: &'g Graph<u64>) -> Self {
        Self {
            graph,
            capacity: graph.edges().iter().map(|edge| [edge.weight, 0]).collect(),
            value: 0,
        }
    }

    fn capacity(&self, arc: Residual) -> u64 {
        self.capacity[arc.edge][arc.reversed]
    }

    fn push(&mut self, arc: Residual, amount: u64) {
        self.capacity[arc.edge][arc.reversed] -= amount;
        self.capacity[arc.edge][arc.reversed ^ 1] += amount;
    }

    fn augment(&mut self, path: &[Residual]) {
        let amount = path.iter().map(|&arc| self.capacity(arc)).min().unwrap();
        for &arc in path {
            self.push(arc, amount);
        }
        self.value += u128::from(amount);
    }

    fn levels(&self, source: NodeId, level: &mut [usize], queue: &mut Vec<NodeId>) {
        level.fill(NONE);
        level[source] = 0;
        queue.clear();
        queue.push(source);
        let mut head = 0;
        while let Some(&node) = queue.get(head) {
            head += 1;
            for arc in residual_neighbors(self.graph, node) {
                if level[arc.to] == NONE && self.capacity(arc) > 0 {
                    level[arc.to] = level[node] + 1;
                    queue.push(arc.to);
                }
            }
        }
    }

    fn finish(self, source: NodeId) -> MaxFlow {
        let mut level = vec![NONE; self.graph.node_count()];
        self.levels(source, &mut level, &mut Vec::new());
        MaxFlow {
            value: self.value,
            edge_flows: self.capacity.into_iter().map(|pair| pair[1]).collect(),
            source_side: level.into_iter().map(|rank| rank != NONE).collect(),
        }
    }
}

struct Trees<'g> {
    residual: Network<'g>,
    source: NodeId,
    sink: NodeId,
    side: Vec<u8>,
    // (parent vertex, residual arc oriented along a source-to-sink path)
    parent: Vec<Option<(NodeId, Residual)>>,
    queued: Vec<bool>,
    active: VecDeque<NodeId>,
    orphan: VecDeque<NodeId>,
    timestamp: Vec<usize>,
    epoch: usize,
    distance: Vec<usize>,
    ancestors: Vec<NodeId>,
}

impl Trees<'_> {
    fn activate(&mut self, node: NodeId) {
        if !self.queued[node] {
            self.queued[node] = true;
            self.active.push_back(node);
        }
    }

    fn grow(&mut self) -> Option<(NodeId, Residual)> {
        while let Some(node) = self.active.pop_front() {
            self.queued[node] = false;
            let side = self.side[node];
            if side == 0 {
                continue;
            }
            for arc in residual_neighbors(self.residual.graph, node) {
                let neighbor = arc.to;
                let oriented = if side == 1 {
                    arc
                } else {
                    Residual {
                        to: node,
                        reversed: arc.reversed ^ 1,
                        ..arc
                    }
                };
                if self.residual.capacity(oriented) == 0 {
                    continue;
                }
                if self.side[neighbor] == 0 {
                    self.side[neighbor] = side;
                    self.parent[neighbor] = Some((node, oriented));
                    self.activate(neighbor);
                } else if self.side[neighbor] != side {
                    self.activate(node);
                    return Some((if side == 1 { node } else { neighbor }, oriented));
                }
            }
        }
        None
    }

    fn augment(&mut self, start: NodeId, crossing: Residual) {
        let endpoints = [(start, self.source), (crossing.to, self.sink)];
        let mut amount = self.residual.capacity(crossing);
        for (mut node, root) in endpoints {
            while node != root {
                let (parent, arc) = self.parent[node].unwrap();
                amount = amount.min(self.residual.capacity(arc));
                node = parent;
            }
        }
        self.residual.push(crossing, amount);
        self.residual.value += u128::from(amount);
        for (mut node, root) in endpoints {
            while node != root {
                let (parent, arc) = self.parent[node].unwrap();
                self.residual.push(arc, amount);
                if self.residual.capacity(arc) == 0 {
                    self.orphan.push_back(node);
                    self.parent[node] = None;
                }
                node = parent;
            }
        }
    }

    fn adopt(&mut self) {
        if self.orphan.is_empty() {
            return;
        }
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.timestamp.fill(0);
            self.epoch = 1;
        }
        for root in [self.source, self.sink] {
            self.timestamp[root] = self.epoch;
            self.distance[root] = 0;
        }
        while let Some(node) = self.orphan.pop_front() {
            let side = self.side[node];
            let mut best = None;
            let mut shortest = NONE;
            for arc in residual_neighbors(self.residual.graph, node) {
                if self.side[arc.to] != side || arc.to == node {
                    continue;
                }
                let oriented = if side == 1 {
                    Residual {
                        to: node,
                        reversed: arc.reversed ^ 1,
                        ..arc
                    }
                } else {
                    arc
                };
                if self.residual.capacity(oriented) == 0 {
                    continue;
                }
                if let Some(distance) = self.distance_to_root(arc.to) {
                    if distance < shortest {
                        shortest = distance;
                        best = Some((arc.to, oriented));
                    }
                }
            }
            if let Some(parent) = best {
                self.parent[node] = Some(parent);
                self.timestamp[node] = self.epoch;
                self.distance[node] = shortest + 1;
                continue;
            }
            self.side[node] = 0;
            for arc in residual_neighbors(self.residual.graph, node) {
                // Either search may now grow through the released vertex.
                if self.side[arc.to] != 0 {
                    self.activate(arc.to);
                }
                if self.parent[arc.to].is_some_and(|(parent, _)| parent == node) {
                    self.parent[arc.to] = None;
                    self.orphan.push_back(arc.to);
                }
            }
        }
    }

    fn distance_to_root(&mut self, vertex: NodeId) -> Option<usize> {
        self.ancestors.clear();
        let mut v = vertex;
        while self.timestamp[v] != self.epoch {
            let (parent, _) = self.parent[v]?;
            self.ancestors.push(v);
            v = parent;
        }
        let mut depth = self.distance[v];
        while let Some(node) = self.ancestors.pop() {
            depth += 1;
            self.timestamp[node] = self.epoch;
            self.distance[node] = depth;
        }
        Some(depth)
    }
}
