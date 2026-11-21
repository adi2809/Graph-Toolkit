use std::collections::HashMap;
use std::error::Error;
use std::fmt;

use crate::graph::validate_node;
use crate::{GraphError, NodeId};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConnectivityEvent {
    Add(NodeId, NodeId),
    Remove(NodeId, NodeId),
    Connected(NodeId, NodeId),
    ComponentCount,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConnectivityAnswer {
    Connected(bool),
    ComponentCount(usize),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DynamicConnectivityError {
    Graph(GraphError),
    InactiveEdge {
        event: usize,
        from: NodeId,
        to: NodeId,
    },
}

impl fmt::Display for DynamicConnectivityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Graph(error) => error.fmt(f),
            Self::InactiveEdge { event, from, to } => {
                write!(f, "event {event} removes inactive edge ({from}, {to})")
            }
        }
    }
}
impl Error for DynamicConnectivityError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Graph(error) => Some(error),
            Self::InactiveEdge { .. } => None,
        }
    }
}

/// Offline undirected connectivity with insertions, deletions, and queries.
/// The graph starts empty. An edge can be inserted multiple times and stays
/// active until its final copy is removed. Query answers follow input order.
///
/// A segment tree groups edge lifetimes. A rollback union–find processes
/// every group in O(T log T log V) time and O(V + T log T) space.
pub fn offline_dynamic_connectivity(
    n: usize,
    events: &[ConnectivityEvent],
) -> Result<Vec<ConnectivityAnswer>, DynamicConnectivityError> {
    use ConnectivityEvent::{Add, ComponentCount, Connected, Remove};
    let base = events
        .len()
        .checked_next_power_of_two()
        .and_then(|base| base.checked_mul(2))
        .ok_or(DynamicConnectivityError::Graph(
            GraphError::GraphSizeOverflow,
        ))?
        / 2;
    let mut tree = vec![Vec::new(); base * 2];
    let mut active = HashMap::<(usize, usize), (usize, usize)>::new();
    for (time, event) in events.iter().copied().enumerate() {
        let (a, b) = match event {
            Add(a, b) | Remove(a, b) | Connected(a, b) => (a, b),
            ComponentCount => continue,
        };
        for node in [a, b] {
            validate_node(node, n).map_err(DynamicConnectivityError::Graph)?;
        }
        let edge = (a.min(b), a.max(b));
        match event {
            Add(..) => {
                let state = active.entry(edge).or_insert((time, 0));
                state.1 += 1;
            }
            Remove(..) => {
                let Some(state) = active.get_mut(&edge) else {
                    return Err(DynamicConnectivityError::InactiveEdge {
                        event: time,
                        from: a,
                        to: b,
                    });
                };
                state.1 -= 1;
                if state.1 == 0 {
                    let (start, _) = active.remove(&edge).unwrap();
                    add_lifetime(&mut tree, start + base, time + base, edge);
                }
            }
            _ => {}
        }
    }
    for (edge, (start, _)) in active {
        add_lifetime(&mut tree, start + base, events.len() + base, edge);
    }
    let mut answers = Vec::new();
    let mut sets = RollbackUnionFind::new(n);
    process(&tree, events, base, 1, &mut sets, &mut answers);
    Ok(answers)
}

fn add_lifetime(
    tree: &mut [Vec<(usize, usize)>],
    mut left: usize,
    mut right: usize,
    edge: (usize, usize),
) {
    while left < right {
        if left % 2 == 1 {
            tree[left].push(edge);
            left += 1;
        }
        if right % 2 == 1 {
            right -= 1;
            tree[right].push(edge);
        }
        left /= 2;
        right /= 2;
    }
}

fn process(
    tree: &[Vec<(usize, usize)>],
    events: &[ConnectivityEvent],
    base: usize,
    node: usize,
    sets: &mut RollbackUnionFind,
    answers: &mut Vec<ConnectivityAnswer>,
) {
    let checkpoint = sets.checkpoint();
    for &(a, b) in &tree[node] {
        sets.union(a, b);
    }
    if node < base {
        process(tree, events, base, 2 * node, sets, answers);
        process(tree, events, base, 2 * node + 1, sets, answers);
    } else {
        match events.get(node - base) {
            Some(ConnectivityEvent::Connected(a, b)) => {
                answers.push(ConnectivityAnswer::Connected(
                    sets.find(*a) == sets.find(*b),
                ));
            }
            Some(ConnectivityEvent::ComponentCount) => {
                answers.push(ConnectivityAnswer::ComponentCount(sets.count));
            }
            _ => {}
        }
    }
    sets.rollback(checkpoint);
}

/// Union by size with O(log n) finds and reversible merges. A checkpoint
/// denotes a prefix of this instance's current merge history.
pub struct RollbackUnionFind {
    parent: Vec<isize>,
    history: Vec<(usize, usize, isize)>,
    count: usize,
}

impl RollbackUnionFind {
    pub fn new(n: usize) -> Self {
        Self {
            parent: vec![-1; n],
            history: Vec::new(),
            count: n,
        }
    }
    pub fn component_count(&self) -> usize {
        self.count
    }
    pub fn checkpoint(&self) -> usize {
        self.history.len()
    }
    pub fn find(&self, mut node: usize) -> usize {
        while self.parent[node] >= 0 {
            node = self.parent[node] as usize;
        }
        node
    }

    pub fn union(&mut self, a: usize, b: usize) {
        let (mut a, mut b) = (self.find(a), self.find(b));
        if a == b {
            return;
        }
        if self.parent[a] > self.parent[b] {
            std::mem::swap(&mut a, &mut b);
        }
        self.history.push((a, b, self.parent[b]));
        self.parent[a] += self.parent[b];
        self.parent[b] = a as isize;
        self.count -= 1;
    }

    /// Undoes all merges after a checkpoint on the current history. A depth
    /// greater than the current depth leaves the sets unchanged.
    pub fn rollback(&mut self, checkpoint: usize) {
        while self.history.len() > checkpoint {
            let (a, b, size) = self.history.pop().unwrap();
            self.parent[a] -= size;
            self.parent[b] = size;
            self.count += 1;
        }
    }
}
