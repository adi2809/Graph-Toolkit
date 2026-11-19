use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};

use super::{INF, NO_LINK, Path, trace_path};
use crate::queue::RadixHeap;
use crate::{EdgeId, Graph, GraphError, NodeId};

/// Returns up to k lowest-cost simple paths, in nondecreasing cost order.
/// Yen's spur search uses indexed radix Dijkstra and a prefix trie of
/// previously emitted paths. Parallel arcs yield distinct edge-ID paths.
pub fn k_shortest_simple_paths(
    graph: &Graph<u64>,
    source: NodeId,
    target: NodeId,
    k: usize,
) -> Result<Vec<Path>, GraphError> {
    graph.validate_node(source)?;
    graph.validate_node(target)?;
    if k == 0 {
        return Ok(Vec::new());
    }

    let mut search = SpurSearch::new(graph);
    let mut blocked_nodes = vec![false; graph.node_count()];
    let mut blocked_edges = vec![false; graph.edge_count()];
    let Some(first) = search.find(source, target, &blocked_nodes, &blocked_edges) else {
        return Ok(Vec::new());
    };
    let mut seen = HashSet::from([first.edges.clone()]);
    let mut result = vec![first];
    let mut trie = vec![HashMap::<EdgeId, usize>::new()];
    type Candidate = Reverse<(u128, Vec<EdgeId>, Vec<NodeId>)>;
    let mut candidates = BinaryHeap::<Candidate>::new();

    while result.len() < k {
        let path = result.last().unwrap();
        let mut prefix = 0;
        for &edge in &path.edges {
            let new = trie.len();
            prefix = *trie[prefix].entry(edge).or_insert(new);
            if prefix == new {
                trie.push(HashMap::new());
            }
        }

        let mut prefix = 0;
        let mut prefix_cost = 0;
        for (position, &edge) in path.edges.iter().enumerate() {
            let spur_node = path.nodes[position];
            for &excluded in trie[prefix].keys() {
                blocked_edges[excluded] = true;
            }
            if let Some(spur) = search.find(spur_node, target, &blocked_nodes, &blocked_edges) {
                let edges: Vec<_> = path.edges[..position]
                    .iter()
                    .copied()
                    .chain(spur.edges)
                    .collect();
                if seen.insert(edges.clone()) {
                    let nodes = path.nodes[..position]
                        .iter()
                        .copied()
                        .chain(spur.nodes)
                        .collect();
                    candidates.push(Reverse((prefix_cost + spur.cost, edges, nodes)));
                }
            }
            for &excluded in trie[prefix].keys() {
                blocked_edges[excluded] = false;
            }
            blocked_nodes[spur_node] = true;
            prefix_cost += u128::from(graph.edges()[edge].weight);
            prefix = trie[prefix][&edge];
        }
        for &node in &path.nodes[..path.edges.len()] {
            blocked_nodes[node] = false;
        }

        let Some(Reverse((cost, edges, nodes))) = candidates.pop() else {
            break;
        };
        result.push(Path { cost, nodes, edges });
    }

    Ok(result)
}

struct SpurSearch<'g> {
    graph: &'g Graph<u64>,
    distance: Vec<u128>,
    parent: Vec<(NodeId, EdgeId)>,
    queue: RadixHeap,
    touched: Vec<NodeId>,
}

impl<'g> SpurSearch<'g> {
    fn new(graph: &'g Graph<u64>) -> Self {
        let n = graph.node_count();
        Self {
            graph,
            distance: vec![INF; n],
            parent: vec![NO_LINK; n],
            queue: RadixHeap::new(n),
            touched: Vec::new(),
        }
    }

    fn find(
        &mut self,
        source: NodeId,
        target: NodeId,
        nodes: &[bool],
        edges: &[bool],
    ) -> Option<Path> {
        self.queue.reset();
        for node in self.touched.drain(..) {
            self.distance[node] = INF;
        }
        self.touched.push(source);
        self.distance[source] = 0;
        self.parent[source] = NO_LINK;
        self.queue.decrease(source, INF, 0);
        while let Some((cost, node)) = self.queue.pop(&self.distance) {
            if node == target {
                let (nodes, edges) = trace_path(&self.parent, target);
                return Some(Path { cost, nodes, edges });
            }
            for arc in self.graph.out(node) {
                if nodes[arc.node] || edges[arc.edge] {
                    continue;
                }
                let candidate = cost + u128::from(self.graph.edges()[arc.edge].weight);
                let old = self.distance[arc.node];
                if candidate < old {
                    if old == INF {
                        self.touched.push(arc.node);
                    }
                    self.queue.decrease(arc.node, old, candidate);
                    self.distance[arc.node] = candidate;
                    self.parent[arc.node] = (node, arc.edge);
                }
            }
        }
        None
    }
}
