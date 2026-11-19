//! Primal network simplex with block pricing and a strongly feasible basis.

use std::ops::{Add, Neg, Sub};

use super::{ArcCost, CostFlow, CostFlowError, flow_cost, min_cost::validate_supply};
use crate::{Graph, NONE, NodeId};

/// Exact minimum-cost flow with finite capacities. Positive supplies send
/// flow. A nonzero artificial flow at termination indicates infeasibility.
///
/// The primal network simplex uses cyclic block pricing and the last blocking
/// arc rule for degenerate pivots. Self-loops are optimized independently.
/// Like other primal simplex variants it can take exponentially many pivots.
pub fn network_simplex(
    graph: &Graph<ArcCost>,
    supply: &[i64],
) -> Result<Option<CostFlow>, CostFlowError> {
    if !validate_supply(graph, supply)? {
        return Ok(None);
    }
    let n = graph.node_count();

    let mut edge_flows = vec![0; graph.edge_count()];
    let mut original = Vec::new();
    let mut arcs = Vec::new();
    for (id, edge) in graph.edges().iter().enumerate() {
        if edge.from == edge.to {
            if edge.weight.unit_cost < 0 {
                edge_flows[id] = edge.weight.capacity;
            }
        } else if edge.weight.capacity > 0 {
            original.push(id);
            arcs.push(BasisArc {
                from: edge.from,
                to: edge.to,
                capacity: edge.weight.capacity.into(),
                cost: Price(0, edge.weight.unit_cost.into()),
                flow: 0,
                basic: false,
            });
        }
    }
    let real = arcs.len();

    let mut tree = Tree {
        parent: vec![n; n + 1],
        edge: (real..real + n).chain([NONE]).collect(),
        depth: vec![1; n + 1],
        price: vec![Price::ZERO; n + 1],
        adjacency: (0..=n).map(|_| Vec::new()).collect(),
        position: vec![[0, 0]; real + n],
    };
    tree.parent[n] = NONE;
    tree.depth[n] = 0;
    for (node, &supply) in supply.iter().enumerate() {
        // At zero supply, orient the artificial arc away from the root.
        // The root must be able to send positive residual flow to every vertex.
        let (from, to, sign) = if supply > 0 {
            (node, n, 1)
        } else {
            (n, node, -1)
        };
        tree.price[node] = Price(sign, 0);
        let id = arcs.len();
        arcs.push(BasisArc {
            from,
            to,
            cost: Price(1, 0),
            capacity: u128::MAX,
            flow: supply.unsigned_abs().into(),
            basic: true,
        });
        tree.add_arc(id, arcs[id]);
    }

    let mut pricing = Pricing::new(arcs.len());
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut stack = Vec::new();
    while let Some((entry, start)) = pricing.entering(&arcs, &tree.price) {
        let end = arcs[entry].other(start);
        tree.cycle(entry, start, end, &mut left, &mut right);
        let (leave, amount) = left
            .iter()
            .rev()
            .map(|&(id, tail)| (id, arcs[id].residual(tail)))
            .min_by_key(|&(_, residual)| residual)
            .unwrap();

        for &(id, tail) in &left {
            let arc = &mut arcs[id];
            if tail == arc.from {
                arc.flow += amount;
            } else {
                arc.flow -= amount;
            }
        }

        if leave != entry {
            arcs[leave].basic = false;
            arcs[entry].basic = true;
            let old = arcs[leave];
            let detached = if tree.edge[old.from] == leave {
                old.from
            } else {
                old.to
            };
            tree.remove_arc(leave, &arcs);

            // Exactly one entry endpoint is in the detached subtree. Reroot
            // by scanning only that subtree, updating depths and potentials.
            let mut ancestor = start;
            while ancestor != detached && tree.depth[ancestor] > tree.depth[detached] {
                ancestor = tree.parent[ancestor];
            }
            let (parent, child) = if ancestor == detached {
                (end, start)
            } else {
                (start, end)
            };
            tree.add_arc(entry, arcs[entry]);
            stack.push((parent, child, entry));
            while let Some((parent, node, id)) = stack.pop() {
                tree.parent[node] = parent;
                tree.edge[node] = id;
                tree.depth[node] = tree.depth[parent] + 1;
                let arc = arcs[id];
                let delta = if arc.from == node {
                    arc.cost
                } else {
                    -arc.cost
                };
                tree.price[node] = tree.price[parent] + delta;
                for &neighbor in &tree.adjacency[node] {
                    if neighbor != id {
                        stack.push((node, arcs[neighbor].other(node), neighbor));
                    }
                }
            }
        }
    }

    if arcs[real..].iter().any(|arc| arc.flow != 0) {
        return Ok(None);
    }
    for (&id, arc) in original.iter().zip(arcs) {
        edge_flows[id] = arc.flow as u64;
    }
    Ok(Some(CostFlow {
        total_cost: flow_cost(graph, &edge_flows)?,
        edge_flows,
    }))
}

// Lexicographic prices avoid a synthetic big-M cost and its overflow risks.
// A simple tree path has at most V arcs, each with an i64 real cost.
#[derive(Clone, Copy, Eq, Ord, PartialEq, PartialOrd)]
struct Price(i128, i128);

impl Price {
    const ZERO: Self = Self(0, 0);
}

impl Add for Price {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0, self.1 + rhs.1)
    }
}

impl Sub for Price {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0, self.1 - rhs.1)
    }
}

impl Neg for Price {
    type Output = Self;

    fn neg(self) -> Self {
        Self(-self.0, -self.1)
    }
}

#[derive(Clone, Copy)]
struct BasisArc {
    from: NodeId,
    to: NodeId,
    capacity: u128,
    cost: Price,
    flow: u128,
    basic: bool,
}

impl BasisArc {
    fn other(self, node: NodeId) -> NodeId {
        if node == self.from {
            self.to
        } else {
            self.from
        }
    }

    fn residual(self, tail: NodeId) -> u128 {
        if tail == self.from {
            self.capacity - self.flow
        } else {
            self.flow
        }
    }
}

struct Tree {
    parent: Vec<NodeId>,
    edge: Vec<usize>,
    depth: Vec<usize>,
    price: Vec<Price>,
    adjacency: Vec<Vec<usize>>,
    position: Vec<[usize; 2]>,
}

impl Tree {
    fn add_arc(&mut self, id: usize, arc: BasisArc) {
        for (side, node) in [arc.from, arc.to].into_iter().enumerate() {
            self.position[id][side] = self.adjacency[node].len();
            self.adjacency[node].push(id);
        }
    }

    fn remove_arc(&mut self, id: usize, arcs: &[BasisArc]) {
        for (side, node) in [arcs[id].from, arcs[id].to].into_iter().enumerate() {
            let position = self.position[id][side];
            let adjacent = &mut self.adjacency[node];
            debug_assert_eq!(adjacent[position], id);
            adjacent.swap_remove(position);
            if let Some(&moved) = adjacent.get(position) {
                self.position[moved][usize::from(arcs[moved].to == node)] = position;
            }
        }
    }

    fn cycle(
        &self,
        entry: usize,
        start: NodeId,
        end: NodeId,
        left: &mut Vec<(usize, NodeId)>,
        right: &mut Vec<(usize, NodeId)>,
    ) {
        left.clear();
        right.clear();
        let (mut a, mut b) = (start, end);
        while a != b {
            if self.depth[a] >= self.depth[b] {
                left.push((self.edge[a], self.parent[a]));
                a = self.parent[a];
            } else {
                right.push((self.edge[b], b));
                b = self.parent[b];
            }
        }
        // Follow the cycle from the apex, through entry, and back to the
        // apex so the final blocker preserves strong feasibility.
        left.reverse();
        left.push((entry, start));
        left.append(right);
    }
}

struct Pricing {
    cursor: usize,
    block: usize,
}

impl Pricing {
    fn new(count: usize) -> Self {
        let root = count.isqrt();
        Self {
            cursor: 0,
            block: root + usize::from(root * root != count),
        }
    }

    fn entering(&mut self, arcs: &[BasisArc], prices: &[Price]) -> Option<(usize, NodeId)> {
        let mut remaining = arcs.len();
        while remaining > 0 {
            let length = remaining.min(self.block);
            let mut best = None;
            let mut best_cost = Price::ZERO;
            for _ in 0..length {
                let id = self.cursor;
                self.cursor += 1;
                if self.cursor == arcs.len() {
                    self.cursor = 0;
                }
                let arc = arcs[id];
                if arc.basic {
                    continue;
                }
                let reduced = arc.cost - prices[arc.from] + prices[arc.to];
                let (reduced, tail) = if arc.flow == 0 {
                    (reduced, arc.from)
                } else {
                    debug_assert_eq!(arc.flow, arc.capacity);
                    (-reduced, arc.to)
                };
                if reduced < best_cost {
                    best = Some((id, tail));
                    best_cost = reduced;
                }
            }
            if best.is_some() {
                return best;
            }
            remaining -= length;
        }
        None
    }
}
