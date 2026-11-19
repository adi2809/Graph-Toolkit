use std::error::Error;
use std::fmt;

/// A dense, zero-based vertex identifier.
pub type NodeId = usize;
/// The zero-based insertion index of a logical edge.
pub type EdgeId = usize;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GraphKind {
    Directed,
    Undirected,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Edge<W> {
    pub from: NodeId,
    pub to: NodeId,
    pub weight: W,
}

/// One incidence in an adjacency list. For incoming incidences, `node` is the origin.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct Arc {
    pub node: NodeId,
    pub edge: EdgeId,
}

#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum GraphError {
    InvalidNode { node: NodeId, node_count: usize },
    InvalidEdge { edge: EdgeId, edge_count: usize },
    GraphSizeOverflow,
    InvalidModelParameters,
    RequiresDirected,
    RequiresUndirected,
    IdenticalTerminals,
    InvalidBipartition { edge: EdgeId },
    NonBinaryWeight { edge: EdgeId, weight: u64 },
    WeightOverflow,
    Cycle,
}

impl fmt::Display for GraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::InvalidNode { node, node_count } => {
                write!(f, "vertex {node} is outside 0..{node_count}")
            }
            Self::InvalidEdge { edge, edge_count } => {
                write!(f, "edge {edge} is outside 0..{edge_count}")
            }
            Self::GraphSizeOverflow => f.write_str("the requested graph size exceeds index limits"),
            Self::InvalidModelParameters => f.write_str("invalid graph model parameters"),
            Self::RequiresDirected => f.write_str("the algorithm requires a directed graph"),
            Self::RequiresUndirected => f.write_str("the algorithm requires an undirected graph"),
            Self::IdenticalTerminals => f.write_str("source and sink must be distinct"),
            Self::InvalidBipartition { edge } => {
                write!(
                    f,
                    "edge {edge} has both endpoints in the same bipartite partition"
                )
            }
            Self::NonBinaryWeight { edge, weight } => {
                write!(f, "edge {edge} has weight {weight}, expected 0 or 1")
            }
            Self::WeightOverflow => f.write_str("the total edge weight does not fit i128"),
            Self::Cycle => f.write_str("the graph contains a directed cycle"),
        }
    }
}

impl Error for GraphError {}

/// Mutable construction phase for a graph with a known vertex count.
#[derive(Clone, Debug)]
pub struct GraphBuilder<W> {
    kind: GraphKind,
    node_count: usize,
    edges: Vec<Edge<W>>,
}

impl<W> GraphBuilder<W> {
    pub fn new(node_count: usize, kind: GraphKind) -> Self {
        Self::with_edge_capacity(node_count, kind, 0)
    }

    pub fn with_edge_capacity(node_count: usize, kind: GraphKind, capacity: usize) -> Self {
        Self {
            kind,
            node_count,
            edges: Vec::with_capacity(capacity),
        }
    }

    /// Preserves parallel edges and self-loops. Failed insertions leave the builder unchanged.
    pub fn add_edge(&mut self, from: NodeId, to: NodeId, weight: W) -> Result<EdgeId, GraphError> {
        validate_node(from, self.node_count)?;
        validate_node(to, self.node_count)?;
        let id = self.edges.len();
        self.edges.push(Edge { from, to, weight });
        Ok(id)
    }

    /// Appends consecutive edges. Every vertex is validated before insertion.
    pub fn add_path(&mut self, nodes: &[NodeId], weight: W) -> Result<(), GraphError>
    where
        W: Clone,
    {
        self.validate_vertices(nodes)?;
        for pair in nodes.windows(2) {
            self.edges.push(Edge {
                from: pair[0],
                to: pair[1],
                weight: weight.clone(),
            });
        }
        Ok(())
    }

    /// Appends a path and its closing edge. A singleton creates a self-loop.
    pub fn add_cycle(&mut self, nodes: &[NodeId], weight: W) -> Result<(), GraphError>
    where
        W: Clone,
    {
        self.add_path(nodes, weight.clone())?;
        if let (Some(&first), Some(&last)) = (nodes.first(), nodes.last()) {
            self.edges.push(Edge {
                from: last,
                to: first,
                weight,
            });
        }
        Ok(())
    }

    /// Appends one edge from the first vertex to each remaining vertex.
    pub fn add_star(&mut self, nodes: &[NodeId], weight: W) -> Result<(), GraphError>
    where
        W: Clone,
    {
        self.validate_vertices(nodes)?;
        if let Some((&center, leaves)) = nodes.split_first() {
            self.edges.extend(leaves.iter().map(|&to| Edge {
                from: center,
                to,
                weight: weight.clone(),
            }));
        }
        Ok(())
    }

    fn validate_vertices(&self, nodes: &[NodeId]) -> Result<(), GraphError> {
        for &node in nodes {
            validate_node(node, self.node_count)?;
        }
        Ok(())
    }

    /// Builds contiguous adjacency lists in insertion order in O(V + E) time.
    pub fn build(self) -> Graph<W> {
        let outgoing = Csr::build(self.node_count, &self.edges, self.kind, false);
        let incoming = (self.kind == GraphKind::Directed)
            .then(|| Csr::build(self.node_count, &self.edges, self.kind, true));
        Graph {
            kind: self.kind,
            edges: self.edges,
            outgoing,
            incoming,
        }
    }
}

/// Immutable compressed-sparse-row multigraph. Edge payloads are stored once.
/// Directed graphs precompute reverse adjacency for queries against incoming arcs.
#[derive(Clone, Debug)]
pub struct Graph<W> {
    kind: GraphKind,
    edges: Vec<Edge<W>>,
    outgoing: Csr,
    incoming: Option<Csr>,
}

impl<W> Graph<W> {
    pub fn builder(node_count: usize, kind: GraphKind) -> GraphBuilder<W> {
        GraphBuilder::new(node_count, kind)
    }

    pub fn kind(&self) -> GraphKind {
        self.kind
    }

    pub fn node_count(&self) -> usize {
        self.outgoing.offsets.len() - 1
    }

    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Logical edges in insertion order. The slice index is the edge ID.
    pub fn edges(&self) -> &[Edge<W>] {
        &self.edges
    }

    pub fn edge(&self, id: EdgeId) -> Option<&Edge<W>> {
        self.edges.get(id)
    }

    /// Undirected self-loops appear once.
    pub fn outgoing(&self, node: NodeId) -> Result<&[Arc], GraphError> {
        self.validate_node(node)?;
        Ok(self.out(node))
    }

    /// For an undirected graph this returns the same slice as outgoing adjacency.
    pub fn incoming(&self, node: NodeId) -> Result<&[Arc], GraphError> {
        self.validate_node(node)?;
        Ok(self.inc(node))
    }

    pub(crate) fn out(&self, node: NodeId) -> &[Arc] {
        self.outgoing.at(node)
    }

    pub(crate) fn inc(&self, node: NodeId) -> &[Arc] {
        self.incoming.as_ref().unwrap_or(&self.outgoing).at(node)
    }

    pub(crate) fn validate_node(&self, node: NodeId) -> Result<(), GraphError> {
        validate_node(node, self.node_count())
    }

    pub(crate) fn require(&self, kind: GraphKind) -> Result<(), GraphError> {
        if self.kind == kind {
            return Ok(());
        }
        Err(match kind {
            GraphKind::Directed => GraphError::RequiresDirected,
            GraphKind::Undirected => GraphError::RequiresUndirected,
        })
    }
}

pub(crate) fn validate_node(node: NodeId, node_count: usize) -> Result<(), GraphError> {
    if node < node_count {
        Ok(())
    } else {
        Err(GraphError::InvalidNode { node, node_count })
    }
}

#[derive(Clone, Debug)]
struct Csr {
    offsets: Vec<usize>,
    arcs: Vec<Arc>,
}

impl Csr {
    fn build<W>(n: usize, edges: &[Edge<W>], kind: GraphKind, reverse: bool) -> Self {
        let mut offsets = vec![0; n + 1];
        let incidences = edges.iter().enumerate().flat_map(|(edge, item)| {
            let (from, to) = if reverse {
                (item.to, item.from)
            } else {
                (item.from, item.to)
            };
            let second = (kind == GraphKind::Undirected && from != to)
                .then_some((to, Arc { node: from, edge }));
            [(from, Arc { node: to, edge })].into_iter().chain(second)
        });
        for (from, _) in incidences.clone() {
            offsets[from + 1] += 1;
        }
        for node in 0..n {
            offsets[node + 1] += offsets[node];
        }
        let mut arcs = vec![Arc { node: 0, edge: 0 }; offsets[n]];
        let mut cursor = offsets[..n].to_vec();
        for (from, arc) in incidences {
            arcs[cursor[from]] = arc;
            cursor[from] += 1;
        }
        Self { offsets, arcs }
    }

    fn at(&self, node: NodeId) -> &[Arc] {
        &self.arcs[self.offsets[node]..self.offsets[node + 1]]
    }
}
