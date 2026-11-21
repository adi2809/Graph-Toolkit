//! Exact integer maximum flow, minimum-cost flow, and minimum cuts.

use std::cmp::Reverse;
use std::collections::HashMap;

use crate::queue::IndexedHeap;
use crate::{EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

mod augment;
mod core_cycles;
pub use core_cycles::{CoreCycleFlow, EmbeddedCycleHierarchy, max_flow_core_cycles};
mod capacity_cost;
mod cut_tree;
mod mean;
mod ratio;
pub use capacity_cost::{
    minimum_cost_flow_capacity_scaling, minimum_cost_flow_minimum_mean_cycle,
    minimum_cost_flow_with_lower_bounds,
};
pub use mean::minimum_mean_cycle;
pub use ratio::{
    CycleRatioWeight, CycleStep, RatioCycle, minimum_ratio_cycle, minimum_undirected_ratio_cycle,
};
mod min_cost;
mod network_simplex;
mod transport;
pub use augment::{
    TreeCycleBasis, TreeCycleFlow, max_flow_boykov_kolmogorov, max_flow_capacity_scaling,
    max_flow_dinic, max_flow_dynamic_trees, max_flow_edmonds_karp, max_flow_ford_fulkerson,
    max_flow_network_simplex, max_flow_st_planar, max_flow_tree_cycle, recover_dyadic_max_flow,
};
pub use cut_tree::{GomoryHuTree, gomory_hu_tree};
pub use min_cost::{
    ArcCost, CostFlow, CostFlowError, MinCostMaxFlow, flow_cost, min_cost_max_flow,
    minimum_cost_flow,
};

pub use network_simplex::network_simplex;
pub use transport::{
    min_cost_max_flow_successive_shortest_path, minimum_cost_flow_cycle_canceling,
    minimum_cost_flow_successive_shortest_path,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MaxFlow {
    pub value: u128,
    /// The flow for each logical directed edge, indexed by its original edge ID.
    pub edge_flows: Vec<u64>,
    /// Membership in the source-side residual reachability set of a minimum cut.
    pub source_side: Vec<bool>,
}

/// Highest-label push–relabel with current arcs, global relabeling and gap detection.
/// Uses the input graph's forward and reverse CSR without materializing residual adjacency.
pub fn max_flow(graph: &Graph<u64>, source: NodeId, sink: NodeId) -> Result<MaxFlow, GraphError> {
    validate_terminals(graph, source, sink)?;
    Ok(Preflow::new(graph, source, sink).solve())
}

fn validate_terminals<W>(graph: &Graph<W>, source: NodeId, sink: NodeId) -> Result<(), GraphError> {
    graph.require(GraphKind::Directed)?;
    graph.validate_node(source)?;
    graph.validate_node(sink)?;
    if source == sink {
        return Err(GraphError::IdenticalTerminals);
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct Residual {
    to: NodeId,
    edge: EdgeId,
    reversed: usize,
}

struct Preflow<'a> {
    graph: &'a Graph<u64>,
    source: NodeId,
    sink: NodeId,
    residual: Vec<[u64; 2]>,
    excess: Vec<u128>,
    height: Vec<usize>,
    current: Vec<usize>,
    height_count: Vec<usize>,
    // Intrusive buckets. Each active vertex has at most one bucket entry.
    head: Vec<NodeId>,
    next: Vec<NodeId>,
    highest: usize,
    work: usize,
    queue: Vec<NodeId>,
}

impl<'a> Preflow<'a> {
    fn new(graph: &'a Graph<u64>, source: NodeId, sink: NodeId) -> Self {
        let n = graph.node_count();
        let mut state = Self {
            graph,
            source,
            sink,
            residual: graph.edges().iter().map(|edge| [edge.weight, 0]).collect(),
            excess: vec![0; n],
            height: vec![0; n],
            current: vec![0; n],
            height_count: vec![0; n],
            head: vec![NONE; 2 * n + 1],
            next: vec![NONE; n],
            highest: 0,
            work: 0,
            queue: Vec::new(),
        };
        for arc in graph.out(source).iter().filter(|arc| arc.node != source) {
            let capacity = state.residual[arc.edge][0];
            state.residual[arc.edge] = [0, capacity];
            state.excess[arc.node] += u128::from(capacity);
        }
        state
    }

    fn solve(mut self) -> MaxFlow {
        let budget = self
            .graph
            .node_count()
            .saturating_add(self.graph.edge_count().saturating_mul(6));
        self.global_relabel();
        while let Some(node) = self.pop_active() {
            self.discharge(node);
            if self.work >= budget {
                self.global_relabel();
            }
        }
        let source_side = self.source_side();
        MaxFlow {
            value: self.excess[self.sink],
            edge_flows: self.residual.into_iter().map(|pair| pair[1]).collect(),
            source_side,
        }
    }

    fn discharge(&mut self, node: NodeId) {
        while self.excess[node] > 0 {
            self.work = self.work.saturating_add(1);
            if let Some(arc) = residual_at(self.graph, node, self.current[node]) {
                let capacity = self.residual[arc.edge][arc.reversed];
                if capacity > 0 && self.height[node] == self.height[arc.to] + 1 {
                    let amount = self.excess[node].min(u128::from(capacity)) as u64;
                    self.residual[arc.edge][arc.reversed] -= amount;
                    self.residual[arc.edge][arc.reversed ^ 1] += amount;
                    self.excess[node] -= u128::from(amount);
                    if arc.to != self.source {
                        let was_inactive = self.excess[arc.to] == 0;
                        self.excess[arc.to] += u128::from(amount);
                        if was_inactive && arc.to != self.sink {
                            self.activate(arc.to);
                        }
                    }
                } else {
                    self.current[node] += 1;
                }
            } else if self.relabel(node) {
                // Gap handling rebuilt all buckets, including the current vertex.
                return;
            }
        }
    }

    fn relabel(&mut self, node: NodeId) -> bool {
        let mut best = usize::MAX;
        for arc in residual_neighbors(self.graph, node) {
            if arc.to != node && self.residual[arc.edge][arc.reversed] != 0 {
                best = best.min(self.height[arc.to]);
            }
        }
        self.work = self
            .work
            .saturating_add(self.graph.out(node).len() + self.graph.inc(node).len());
        debug_assert!(best < 2 * self.graph.node_count());
        let old = self.height[node];
        self.height[node] = best + 1;
        self.current[node] = 0;
        if let Some(count) = self.height_count.get_mut(best + 1) {
            *count += 1;
        }
        if let Some(count) = self.height_count.get_mut(old) {
            *count -= 1;
            if *count == 0 {
                self.apply_gap(old);
                return true;
            }
        }
        false
    }

    fn apply_gap(&mut self, gap: usize) {
        let n = self.graph.node_count();
        for (node, height) in self.height.iter_mut().enumerate() {
            if gap < *height && *height < n {
                self.height_count[*height] -= 1;
                *height = n + 1;
                self.current[node] = 0;
            }
        }
        self.rebuild_active();
        self.work = self.work.saturating_add(n);
    }

    fn global_relabel(&mut self) {
        let n = self.graph.node_count();
        self.height.fill(2 * n);
        self.height[self.source] = n;
        self.height[self.sink] = 0;
        self.reverse_bfs(self.sink);
        // Vertices cut off from the sink can still return excess to the source.
        self.reverse_bfs(self.source);
        self.current.fill(0);
        self.height_count.fill(0);
        for &height in &self.height {
            if let Some(count) = self.height_count.get_mut(height) {
                *count += 1;
            }
        }
        self.rebuild_active();
        self.work = 0;
    }

    fn reverse_bfs(&mut self, start: NodeId) {
        let untouched = 2 * self.graph.node_count();
        self.queue.clear();
        self.queue.push(start);
        let mut head = 0;
        while let Some(&node) = self.queue.get(head) {
            head += 1;
            for arc in residual_neighbors(self.graph, node) {
                if self.height[arc.to] == untouched && self.residual[arc.edge][arc.reversed ^ 1] > 0
                {
                    self.height[arc.to] = self.height[node] + 1;
                    self.queue.push(arc.to);
                }
            }
        }
    }

    fn source_side(&mut self) -> Vec<bool> {
        let mut reached = vec![false; self.graph.node_count()];
        self.queue.clear();
        self.queue.push(self.source);
        reached[self.source] = true;
        let mut head = 0;
        while let Some(&node) = self.queue.get(head) {
            head += 1;
            for arc in residual_neighbors(self.graph, node) {
                if !reached[arc.to] && self.residual[arc.edge][arc.reversed] > 0 {
                    reached[arc.to] = true;
                    self.queue.push(arc.to);
                }
            }
        }
        reached
    }

    fn rebuild_active(&mut self) {
        self.head.fill(NONE);
        self.highest = 0;
        for node in 0..self.graph.node_count() {
            if node != self.source && node != self.sink && self.excess[node] != 0 {
                self.activate(node);
            }
        }
    }

    fn activate(&mut self, node: NodeId) {
        let height = self.height[node];
        self.next[node] = self.head[height];
        self.head[height] = node;
        self.highest = self.highest.max(height);
    }

    fn pop_active(&mut self) -> Option<NodeId> {
        while self.head[self.highest] == NONE {
            self.highest = self.highest.checked_sub(1)?;
        }
        let node = self.head[self.highest];
        self.head[self.highest] = self.next[node];
        Some(node)
    }
}

fn residual_neighbors<W>(graph: &Graph<W>, node: NodeId) -> impl Iterator<Item = Residual> + '_ {
    graph
        .out(node)
        .iter()
        .map(|arc| Residual {
            to: arc.node,
            edge: arc.edge,
            reversed: 0,
        })
        .chain(graph.inc(node).iter().map(|arc| Residual {
            to: arc.node,
            edge: arc.edge,
            reversed: 1,
        }))
}

fn residual_at<W>(graph: &Graph<W>, node: NodeId, index: usize) -> Option<Residual> {
    let outgoing = graph.out(node);
    if let Some(arc) = outgoing.get(index) {
        return Some(Residual {
            to: arc.node,
            edge: arc.edge,
            reversed: 0,
        });
    }
    let arc = graph.inc(node).get(index - outgoing.len())?;
    Some(Residual {
        to: arc.node,
        edge: arc.edge,
        reversed: 1,
    })
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct GlobalMinCut {
    pub value: u128,
    /// Membership in one nonempty side of the cut.
    pub side: Vec<bool>,
}

/// Stoer–Wagner global cut for undirected graphs using sparse contractions and
/// an indexed max-heap. Graphs with fewer than two vertices yield `None`.
pub fn global_min_cut(graph: &Graph<u64>) -> Result<Option<GlobalMinCut>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    if n < 2 {
        return Ok(None);
    }
    let mut adjacency = vec![HashMap::<NodeId, u128>::new(); n];
    for edge in graph
        .edges()
        .iter()
        .filter(|edge| edge.from != edge.to && edge.weight != 0)
    {
        let weight = u128::from(edge.weight);
        *adjacency[edge.from].entry(edge.to).or_default() += weight;
        *adjacency[edge.to].entry(edge.from).or_default() += weight;
    }

    let mut best = u128::MAX;
    let mut best_group = NONE;
    let mut active: Vec<_> = (0..n).collect();
    let mut root = active.clone();
    let mut merges = Vec::with_capacity(n - 1);
    let mut visited = vec![0; n];
    let mut weight = vec![0; n];
    let mut heap = IndexedHeap::new(n);
    for phase in 1..n {
        for &node in &active {
            weight[node] = 0;
            heap.decrease(node, Reverse(0));
        }
        let mut previous = NONE;
        for remaining in (0..active.len()).rev() {
            let (Reverse(cut), node) = heap.pop().unwrap();
            if remaining == 0 {
                if cut < best {
                    best = cut;
                    best_group = root[node];
                }
                if best != 0 {
                    merges.push((root[previous], root[node]));
                    root[previous] = n + merges.len() - 1;
                    contract(&mut adjacency, previous, node);
                    let index = active.iter().position(|&item| item == node).unwrap();
                    active.swap_remove(index);
                }
                break;
            }

            visited[node] = phase;
            for (&neighbor, &capacity) in &adjacency[node] {
                if visited[neighbor] != phase {
                    weight[neighbor] += capacity;
                    heap.decrease(neighbor, Reverse(weight[neighbor]));
                }
            }
            previous = node;
        }
        if best == 0 {
            break;
        }
    }

    let mut side = vec![false; n];
    let mut stack = vec![best_group];
    while let Some(group) = stack.pop() {
        if let Some(&(a, b)) = group.checked_sub(n).and_then(|index| merges.get(index)) {
            stack.extend([a, b]);
        } else {
            side[group] = true;
        }
    }
    Ok(Some(GlobalMinCut { value: best, side }))
}

fn contract(adjacency: &mut [HashMap<NodeId, u128>], keep: NodeId, remove: NodeId) {
    adjacency[keep].remove(&remove);
    for (neighbor, capacity) in std::mem::take(&mut adjacency[remove]) {
        if neighbor != keep {
            adjacency[neighbor].remove(&remove);
            *adjacency[neighbor].entry(keep).or_default() += capacity;
            *adjacency[keep].entry(neighbor).or_default() += capacity;
        }
    }
}
