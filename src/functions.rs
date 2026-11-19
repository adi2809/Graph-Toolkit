//! Inspection, attributed payloads, and read-only subgraphs for CSR graphs.

use crate::{Edge, EdgeId, Graph, GraphBuilder, GraphError, GraphKind, NONE, NodeId};

impl<W> Graph<W> {
    pub fn nodes(&self) -> std::ops::Range<NodeId> {
        0..self.node_count()
    }

    /// Whether there are zero logical edges.
    pub fn is_empty(&self) -> bool {
        self.edge_count() == 0
    }

    pub fn is_directed(&self) -> bool {
        self.kind() == GraphKind::Directed
    }

    /// Graphs have immutable structure after construction.
    pub const fn is_frozen(&self) -> bool {
        true
    }

    /// Neighbor IDs in forward adjacency order, including parallel incidences.
    pub fn neighbors(
        &self,
        node: NodeId,
    ) -> Result<impl ExactSizeIterator<Item = NodeId> + '_, GraphError> {
        Ok(self.outgoing(node)?.iter().map(|arc| arc.node))
    }

    /// Forward and backward neighbors on a digraph, preserving multiplicity.
    pub fn all_neighbors(
        &self,
        node: NodeId,
    ) -> Result<impl Iterator<Item = NodeId> + '_, GraphError> {
        let reverse = if self.is_directed() {
            self.incoming(node)?
        } else {
            &[]
        };
        Ok(self
            .outgoing(node)?
            .iter()
            .chain(reverse)
            .map(|arc| arc.node))
    }

    /// Total degree. Each self-loop contributes two.
    pub fn degree(&self, node: NodeId) -> Result<usize, GraphError> {
        let outgoing = self.outgoing(node)?;
        Ok(outgoing.len()
            + if self.is_directed() {
                self.inc(node).len()
            } else {
                outgoing.iter().filter(|arc| arc.node == node).count()
            })
    }

    /// Vertices without incoming or outgoing incidences.
    pub fn isolates(&self) -> impl Iterator<Item = NodeId> + '_ {
        self.nodes()
            .filter(|&node| self.out(node).is_empty() && self.inc(node).is_empty())
    }

    pub fn degree_histogram(&self) -> Vec<usize> {
        let mut histogram = Vec::new();
        for node in self.nodes() {
            let degree = self.degree(node).unwrap();
            if degree >= histogram.len() {
                histogram.resize(degree + 1, 0);
            }
            histogram[degree] += 1;
        }
        histogram
    }

    /// Logical edges divided by the number of possible simple non-loop edges.
    /// Loops and parallel edges can make the ratio exceed one.
    pub fn density(&self) -> f64 {
        let n = self.node_count();
        if n < 2 {
            return 0.0;
        }
        let divisor = n as f64 * (n - 1) as f64;
        let numerator = self.edge_count() as f64;
        if self.is_directed() {
            numerator / divisor
        } else {
            2.0 * numerator / divisor
        }
    }

    pub fn selfloop_edges(&self) -> impl Iterator<Item = (EdgeId, &Edge<W>)> {
        self.edges()
            .iter()
            .enumerate()
            .filter(|(_, edge)| edge.from == edge.to)
    }

    pub fn number_of_selfloops(&self) -> usize {
        self.selfloop_edges().count()
    }

    pub fn nodes_with_selfloops(&self) -> impl Iterator<Item = NodeId> + '_ {
        let mut seen = vec![false; self.node_count()];
        self.selfloop_edges().filter_map(move |(_, edge)| {
            (!std::mem::replace(&mut seen[edge.from], true)).then_some(edge.from)
        })
    }

    pub fn non_neighbors(
        &self,
        node: NodeId,
    ) -> Result<impl Iterator<Item = NodeId> + '_, GraphError> {
        let mut excluded = vec![false; self.node_count()];
        for arc in self.outgoing(node)? {
            excluded[arc.node] = true;
        }
        excluded[node] = true;
        Ok(self.nodes().filter(move |&other| !excluded[other]))
    }

    /// Common neighbors of two vertices in an undirected graph. Results are
    /// distinct and exclude the endpoints.
    pub fn common_neighbors(&self, a: NodeId, b: NodeId) -> Result<Vec<NodeId>, GraphError> {
        self.require(GraphKind::Undirected)?;
        let (mut first, mut second) = (self.outgoing(a)?, self.outgoing(b)?);
        if first.len() > second.len() {
            std::mem::swap(&mut first, &mut second);
        }
        let mut candidates: std::collections::HashSet<_> =
            first.iter().map(|arc| arc.node).collect();
        candidates.remove(&a);
        candidates.remove(&b);
        Ok(second
            .iter()
            .filter_map(|arc| candidates.remove(&arc.node).then_some(arc.node))
            .collect())
    }

    /// Distinct missing non-loop edges in row-major order. For an undirected
    /// graph only pairs with first vertex less than second are produced.
    pub fn non_edges(&self) -> impl Iterator<Item = (NodeId, NodeId)> + '_ {
        let mut mark = vec![NONE; self.node_count()];
        let (mut from, mut to) = (0, 0);
        let mut loaded = false;
        std::iter::from_fn(move || {
            while from < self.node_count() {
                if !loaded {
                    for arc in self.out(from) {
                        mark[arc.node] = from;
                    }
                    to = if self.is_directed() { 0 } else { from + 1 };
                    loaded = true;
                }
                while to < self.node_count() {
                    let candidate = to;
                    to += 1;
                    if candidate != from && mark[candidate] != from {
                        return Some((from, candidate));
                    }
                }
                from += 1;
                loaded = false;
            }
            None
        })
    }

    /// Checks whether every consecutive pair has an edge in the given
    /// direction. Vertices may repeat. An empty walk is valid.
    pub fn is_path(&self, nodes: &[NodeId]) -> bool {
        nodes.iter().all(|&node| node < self.node_count())
            && nodes
                .windows(2)
                .all(|pair| self.out(pair[0]).iter().any(|arc| arc.node == pair[1]))
    }

    pub fn create_empty_copy(&self) -> Self {
        GraphBuilder::new(self.node_count(), self.kind()).build()
    }

    /// Clones a directed graph, or emits both directed arcs of each undirected
    /// edge, preserving a single logical arc for a self-loop.
    pub fn to_directed(&self) -> Self
    where
        W: Clone,
    {
        if self.is_directed() {
            return self.clone();
        }
        let mut builder = GraphBuilder::new(self.node_count(), GraphKind::Directed);
        for edge in self.edges() {
            builder
                .add_edge(edge.from, edge.to, edge.weight.clone())
                .unwrap();
            if edge.from != edge.to {
                builder
                    .add_edge(edge.to, edge.from, edge.weight.clone())
                    .unwrap();
            }
        }
        builder.build()
    }

    /// Each directed logical edge becomes one undirected edge. Opposite arcs
    /// become parallel edges and retain their distinct payloads.
    pub fn to_undirected(&self) -> Self
    where
        W: Clone,
    {
        if !self.is_directed() {
            return self.clone();
        }
        let mut builder = GraphBuilder::with_edge_capacity(
            self.node_count(),
            GraphKind::Undirected,
            self.edge_count(),
        );
        for edge in self.edges() {
            builder
                .add_edge(edge.from, edge.to, edge.weight.clone())
                .unwrap();
        }
        builder.build()
    }
}

impl Graph<i64> {
    pub fn is_negatively_weighted(&self) -> bool {
        self.edges().iter().any(|edge| edge.weight < 0)
    }

    /// Walk cost, choosing the cheapest parallel edge at each step.
    pub fn path_weight(&self, nodes: &[NodeId]) -> Option<i128> {
        if nodes.iter().any(|&node| node >= self.node_count()) {
            return None;
        }
        nodes.windows(2).try_fold(0, |total, pair| {
            self.out(pair[0])
                .iter()
                .filter(|arc| arc.node == pair[1])
                .map(|arc| i128::from(self.edges()[arc.edge].weight))
                .min()
                .map(|cost| total + cost)
        })
    }

    pub fn weighted_degree(&self, node: NodeId) -> Result<i128, GraphError> {
        let mut total = 0;
        for arc in self.outgoing(node)? {
            let cost = i128::from(self.edges()[arc.edge].weight);
            total += cost;
            if !self.is_directed() && arc.node == node {
                total += cost;
            }
        }
        if self.is_directed() {
            total += self
                .inc(node)
                .iter()
                .map(|arc| i128::from(self.edges()[arc.edge].weight))
                .sum::<i128>();
        }
        Ok(total)
    }
}

/// A stable-ID view over an immutable graph. State consists of inclusion masks.
#[derive(Clone, Debug)]
pub struct GraphView<'g, W> {
    graph: &'g Graph<W>,
    nodes: Vec<bool>,
    edges: Option<Vec<bool>>,
}

/// An owned CSR graph and its original vertex and edge IDs, in new ID order.
#[derive(Clone, Debug)]
pub struct MaterializedGraph<W> {
    pub graph: Graph<W>,
    pub original_nodes: Vec<NodeId>,
    pub original_edges: Vec<EdgeId>,
}

pub fn induced_subgraph<'g, W>(
    graph: &'g Graph<W>,
    selected: impl IntoIterator<Item = NodeId>,
) -> Result<GraphView<'g, W>, GraphError> {
    let mut nodes = vec![false; graph.node_count()];
    for node in selected {
        graph.validate_node(node)?;
        nodes[node] = true;
    }
    Ok(GraphView {
        graph,
        nodes,
        edges: None,
    })
}

pub use induced_subgraph as subgraph;

/// Whether the graph has an edge and every payload exposes a weight through
/// the supplied accessor. Built-in weighted algorithms take typed payloads.
pub fn is_weighted<W, T>(graph: &Graph<W>, weight: impl Fn(&W) -> Option<T>) -> bool {
    !graph.is_empty()
        && graph
            .edges()
            .iter()
            .all(|edge| weight(&edge.weight).is_some())
}

pub fn restricted_view<'g, W>(
    graph: &'g Graph<W>,
    hidden_nodes: impl IntoIterator<Item = NodeId>,
    hidden_edges: impl IntoIterator<Item = EdgeId>,
) -> Result<GraphView<'g, W>, GraphError> {
    let mut nodes = vec![true; graph.node_count()];
    let mut edges = vec![true; graph.edge_count()];
    for node in hidden_nodes {
        graph.validate_node(node)?;
        nodes[node] = false;
    }
    for id in hidden_edges {
        mask_edge(&mut edges, id, false)?;
    }
    Ok(GraphView {
        graph,
        nodes,
        edges: Some(edges),
    })
}

pub fn edge_subgraph<'g, W>(
    graph: &'g Graph<W>,
    selected: impl IntoIterator<Item = EdgeId>,
) -> Result<GraphView<'g, W>, GraphError> {
    let mut nodes = vec![false; graph.node_count()];
    let mut edges = vec![false; graph.edge_count()];
    for id in selected {
        mask_edge(&mut edges, id, true)?;
        let edge = &graph.edges()[id];
        nodes[edge.from] = true;
        nodes[edge.to] = true;
    }
    Ok(GraphView {
        graph,
        nodes,
        edges: Some(edges),
    })
}

fn mask_edge(mask: &mut [bool], edge: EdgeId, value: bool) -> Result<(), GraphError> {
    let edge_count = mask.len();
    *mask
        .get_mut(edge)
        .ok_or(GraphError::InvalidEdge { edge, edge_count })? = value;
    Ok(())
}

impl<W> GraphView<'_, W> {
    pub fn contains_node(&self, node: NodeId) -> bool {
        self.nodes.get(node).copied().unwrap_or(false)
    }

    pub fn nodes(&self) -> impl Iterator<Item = NodeId> + '_ {
        self.graph.nodes().filter(|&node| self.nodes[node])
    }

    pub fn edges(&self) -> impl Iterator<Item = (EdgeId, &Edge<W>)> {
        self.graph.edges().iter().enumerate().filter(|&(id, edge)| {
            self.nodes[edge.from]
                && self.nodes[edge.to]
                && self.edges.as_ref().is_none_or(|mask| mask[id])
        })
    }

    /// Forward incidences of a visible vertex, with original IDs.
    pub fn neighbors(&self, node: NodeId) -> Option<impl Iterator<Item = crate::Arc> + '_> {
        self.contains_node(node).then(|| {
            self.graph.out(node).iter().copied().filter(move |arc| {
                self.nodes[arc.node] && self.edges.as_ref().is_none_or(|mask| mask[arc.edge])
            })
        })
    }

    /// Builds an owned subgraph with dense IDs and both original-ID maps.
    pub fn materialize(&self) -> MaterializedGraph<W>
    where
        W: Clone,
    {
        let original_nodes: Vec<_> = self.nodes().collect();
        let mut index = vec![NONE; self.graph.node_count()];
        for (local, &original) in original_nodes.iter().enumerate() {
            index[original] = local;
        }
        let mut builder = GraphBuilder::new(original_nodes.len(), self.graph.kind());
        let mut original_edges = Vec::new();
        for (id, edge) in self.edges() {
            builder
                .add_edge(index[edge.from], index[edge.to], edge.weight.clone())
                .unwrap();
            original_edges.push(id);
        }
        MaterializedGraph {
            graph: builder.build(),
            original_nodes,
            original_edges,
        }
    }
}

/// Dense typed attribute storage keyed by vertex or by logical edge ID.
/// Separate maps can hold independent attribute types for the same graph.
#[derive(Clone, Debug)]
pub struct AttributeMap<T> {
    values: Vec<Option<T>>,
    domain: AttributeDomain,
}

#[derive(Clone, Copy, Debug)]
enum AttributeDomain {
    Nodes,
    Edges,
}

impl<T> AttributeMap<T> {
    pub fn for_nodes<W>(graph: &Graph<W>) -> Self {
        Self {
            values: std::iter::repeat_with(|| None)
                .take(graph.node_count())
                .collect(),
            domain: AttributeDomain::Nodes,
        }
    }

    pub fn for_edges<W>(graph: &Graph<W>) -> Self {
        Self {
            values: std::iter::repeat_with(|| None)
                .take(graph.edge_count())
                .collect(),
            domain: AttributeDomain::Edges,
        }
    }

    pub fn set(&mut self, id: usize, value: T) -> Result<Option<T>, GraphError> {
        let len = self.values.len();
        let Some(slot) = self.values.get_mut(id) else {
            return Err(match self.domain {
                AttributeDomain::Nodes => GraphError::InvalidNode {
                    node: id,
                    node_count: len,
                },
                AttributeDomain::Edges => GraphError::InvalidEdge {
                    edge: id,
                    edge_count: len,
                },
            });
        };
        Ok(slot.replace(value))
    }

    pub fn get(&self, id: usize) -> Option<&T> {
        self.values.get(id)?.as_ref()
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &T)> {
        self.values
            .iter()
            .enumerate()
            .filter_map(|(id, slot)| Some((id, slot.as_ref()?)))
    }
}
