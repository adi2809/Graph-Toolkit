//! Undirected all-pairs minimum cuts using a Gomory–Hu tree.

use super::max_flow;
use crate::{Graph, GraphError, GraphKind, NONE, NodeId};

#[derive(Clone, Copy, Debug)]
struct Jump {
    parent: NodeId,
    bottleneck: u128,
}

/// A cut-equivalent tree with O(log V) minimum-cut value queries.
#[derive(Clone, Debug)]
pub struct GomoryHuTree {
    tree: Graph<u128>,
    depth: Vec<usize>,
    jumps: Vec<Vec<Jump>>,
}

impl GomoryHuTree {
    /// Tree edge weights are minimum-cut values in the original graph.
    pub fn tree(&self) -> &Graph<u128> {
        &self.tree
    }

    /// Returns the minimum capacity of an undirected cut separating the terminals.
    pub fn minimum_cut_value(&self, mut a: NodeId, mut b: NodeId) -> Result<u128, GraphError> {
        self.tree.validate_node(a)?;
        self.tree.validate_node(b)?;
        if a == b {
            return Err(GraphError::IdenticalTerminals);
        }
        if self.depth[a] < self.depth[b] {
            std::mem::swap(&mut a, &mut b);
        }
        let difference = self.depth[a] - self.depth[b];
        let mut value = u128::MAX;
        for (power, row) in self.jumps.iter().enumerate() {
            if difference & (1 << power) != 0 {
                value = value.min(row[a].bottleneck);
                a = row[a].parent;
            }
        }
        if a == b {
            return Ok(value);
        }
        for row in self.jumps.iter().rev() {
            if row[a].parent != row[b].parent {
                value = value.min(row[a].bottleneck).min(row[b].bottleneck);
                a = row[a].parent;
                b = row[b].parent;
            }
        }
        Ok(value
            .min(self.jumps[0][a].bottleneck)
            .min(self.jumps[0][b].bottleneck))
    }

    fn new(tree: Graph<u128>) -> Self {
        let n = tree.node_count();
        let mut depth = vec![NONE; n];
        let mut jumps = Vec::new();
        if n > 0 {
            let mut row = vec![
                Jump {
                    parent: 0,
                    bottleneck: u128::MAX
                };
                n
            ];
            let mut queue = vec![0];
            depth[0] = 0;
            let mut head = 0;
            while let Some(&node) = queue.get(head) {
                head += 1;
                for arc in tree.out(node) {
                    if depth[arc.node] == NONE {
                        depth[arc.node] = depth[node] + 1;
                        row[arc.node] = Jump {
                            parent: node,
                            bottleneck: tree.edges()[arc.edge].weight,
                        };
                        queue.push(arc.node);
                    }
                }
            }
            jumps.push(row);
            for level in 1..=n.ilog2() {
                let previous = &jumps[level as usize - 1];
                let doubled = previous
                    .iter()
                    .map(|&jump| {
                        let next = previous[jump.parent];
                        Jump {
                            parent: next.parent,
                            bottleneck: jump.bottleneck.min(next.bottleneck),
                        }
                    })
                    .collect();
                jumps.push(doubled);
            }
        }
        Self { tree, depth, jumps }
    }
}

/// Builds a Gomory–Hu tree using Gusfield's algorithm with V − 1 maximum-flow
/// calls on nonempty inputs. Disconnected components are joined by zero cuts.
/// Self-loops are ignored. Parallel capacities contribute separately.
pub fn gomory_hu_tree(graph: &Graph<u64>) -> Result<GomoryHuTree, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    let mut projection = Graph::builder(n, GraphKind::Directed);
    for edge in graph.edges() {
        if edge.from != edge.to && edge.weight > 0 {
            projection.add_edge(edge.from, edge.to, edge.weight)?;
            projection.add_edge(edge.to, edge.from, edge.weight)?;
        }
    }
    let projection = projection.build();
    let mut parent = vec![0; n];
    let mut weights = vec![0; n];
    for source in 1..n {
        let sink = parent[source];
        let flow = max_flow(&projection, source, sink)?;
        for (node, parent) in parent.iter_mut().enumerate().skip(1) {
            if node != source && *parent == sink && flow.source_side[node] {
                *parent = source;
            }
        }
        if sink != 0 && flow.source_side[parent[sink]] {
            parent[source] = parent[sink];
            parent[sink] = source;
            weights[source] = weights[sink];
            weights[sink] = flow.value;
        } else {
            weights[source] = flow.value;
        }
    }
    let mut builder = Graph::builder(n, GraphKind::Undirected);
    for node in 1..n {
        builder.add_edge(node, parent[node], weights[node])?;
    }
    Ok(GomoryHuTree::new(builder.build()))
}
