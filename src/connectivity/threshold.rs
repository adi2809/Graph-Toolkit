use super::{DisjointSet, kruskal};
use crate::{Graph, GraphError, NONE, NodeId, graph::validate_node, structures::RangeMinimumIndex};

/// Connectivity under a maximum allowed edge weight in a static undirected
/// multigraph. Kruskal merges form a binary reconstruction forest. The lowest
/// common ancestor of two leaves is the merge that first connected them.
/// Construction takes O(V + E log(E+1)) time. Queries take O(1).
pub struct ThresholdConnectivityIndex<W> {
    first: Vec<usize>,
    merge_weight: Vec<W>,
    euler: RangeMinimumIndex<(usize, usize)>,
}

impl<W: Copy + Ord> ThresholdConnectivityIndex<W> {
    pub fn new(graph: &Graph<W>) -> Result<Self, GraphError> {
        let (forest, _) = kruskal(graph)?;
        let n = graph.node_count();
        let mut sets = DisjointSet::new(n);
        let mut root: Vec<_> = graph.nodes().collect();
        let mut children = Vec::with_capacity(forest.len());
        let mut weights = Vec::with_capacity(forest.len());
        for id in forest {
            let edge = graph.edges()[id];
            let a = sets.find(edge.from);
            let b = sets.find(edge.to);
            children.push((root[a], root[b]));
            weights.push(edge.weight);
            root[sets.union(a, b).unwrap()] = n + children.len() - 1;
        }
        // The virtual root joins different connected components above all
        // real merges. It also gives an edgeless graph a well-defined tour.
        let virtual_root = n + children.len();
        let mut euler = vec![(0, virtual_root)];
        let mut first = vec![NONE; n];
        let mut stack = Vec::new();
        for v in graph.nodes() {
            if sets.find(v) != v {
                continue;
            }
            stack.push((root[v], 1, false));
            while let Some((node, depth, returning)) = stack.pop() {
                euler.push((depth, node));
                if returning {
                    continue;
                }
                if node < n {
                    first[node] = euler.len() - 1;
                } else {
                    let (left, right) = children[node - n];
                    stack.extend([
                        (node, depth, true),
                        (right, depth + 1, false),
                        (node, depth, true),
                        (left, depth + 1, false),
                    ]);
                }
            }
            euler.push((0, virtual_root));
        }
        Ok(Self {
            first,
            merge_weight: weights,
            euler: RangeMinimumIndex::new(euler),
        })
    }

    /// Includes the zero-edge path from a vertex to itself for every threshold.
    pub fn connected_at(&self, a: NodeId, b: NodeId, threshold: W) -> Result<bool, GraphError> {
        let n = self.first.len();
        validate_node(a, n)?;
        validate_node(b, n)?;
        if a == b {
            return Ok(true);
        }
        Ok(self.merge_weight(a, b).is_some_and(|w| w <= threshold))
    }

    // Distinct, validated leaves. The virtual root has no merge weight.
    pub(super) fn merge_weight(&self, a: usize, b: usize) -> Option<W> {
        let (i, j) = (self.first[a], self.first[b]);
        let &(_, ancestor) = self.euler.minimum(i.min(j)..i.max(j) + 1).unwrap();
        self.merge_weight.get(ancestor - self.first.len()).copied()
    }
}
