use super::{
    StrongComponents, WeakComponents, strongly_connected_components, undirected_low_links,
    weakly_connected_components,
};
use crate::flow::{global_min_cut, max_flow};
use crate::{EdgeId, Graph, GraphBuilder, GraphError, GraphKind, NodeId};

pub fn connected_components<W>(graph: &Graph<W>) -> Result<WeakComponents, GraphError> {
    graph.require(GraphKind::Undirected)?;
    Ok(weakly_connected_components(graph))
}

/// Whether a nonempty undirected graph is connected.
pub fn is_connected<W>(graph: &Graph<W>) -> Result<bool, GraphError> {
    Ok(connected_components(graph)?.count == 1)
}

pub fn is_strongly_connected<W>(graph: &Graph<W>) -> bool {
    strongly_connected_components(graph).count == 1
}

pub fn is_weakly_connected<W>(graph: &Graph<W>) -> bool {
    weakly_connected_components(graph).count == 1
}

/// Test reachability with bidirectional BFS, expanding the frontier with the
/// smaller total adjacency. Traverses incoming arcs from the destination.
pub fn st_connected<W>(
    graph: &Graph<W>,
    source: NodeId,
    target: NodeId,
) -> Result<bool, GraphError> {
    graph.validate_node(source)?;
    graph.validate_node(target)?;
    if source == target {
        return Ok(true);
    }
    let mut sides = vec![0; graph.node_count()];
    sides[source] = 1;
    sides[target] = 2;
    let mut forward = vec![source];
    let mut backward = vec![target];
    let mut scratch = Vec::new();
    let mut incoming = graph.inc(target).len();
    let mut outgoing = graph.out(source).len();
    while !forward.is_empty() && !backward.is_empty() {
        let reverse = incoming < outgoing;
        let (side, frontier, work) = if reverse {
            (2, &mut backward, &mut incoming)
        } else {
            (1, &mut forward, &mut outgoing)
        };
        scratch.clear();
        *work = 0;
        for &node in frontier.iter() {
            let arcs = if reverse {
                graph.inc(node)
            } else {
                graph.out(node)
            };
            for arc in arcs {
                match sides[arc.node] {
                    0 => {
                        sides[arc.node] = side;
                        scratch.push(arc.node);
                        *work += if reverse {
                            graph.inc(arc.node)
                        } else {
                            graph.out(arc.node)
                        }
                        .len();
                    }
                    mark if mark != side => return Ok(true),
                    _ => {}
                }
            }
        }
        std::mem::swap(frontier, &mut scratch);
    }
    Ok(false)
}

#[derive(Clone, Debug)]
pub struct Condensation {
    pub graph: Graph<()>,
    pub components: StrongComponents,
}

/// Contracts strongly connected components and deduplicates inter-component arcs.
pub fn condensation<W>(graph: &Graph<W>) -> Condensation {
    let components = strongly_connected_components(graph);
    let mut members = vec![Vec::new(); components.count];
    for (node, &component) in components.component_of.iter().enumerate() {
        members[component].push(node);
    }
    let mut builder = GraphBuilder::new(components.count, GraphKind::Directed);
    let mut seen = vec![crate::NONE; components.count];
    for (from, nodes) in members.iter().enumerate() {
        for &node in nodes {
            for arc in graph.out(node) {
                let to = components.component_of[arc.node];
                if from != to && seen[to] != from {
                    seen[to] = from;
                    builder.add_edge(from, to, ()).unwrap();
                }
            }
        }
    }
    Condensation {
        graph: builder.build(),
        components,
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BiconnectedComponents {
    /// Each block is given as a list of logical edge IDs. A bridge or self-loop
    /// is a singleton block. Isolated vertices have no edge block.
    pub edge_blocks: Vec<Vec<EdgeId>>,
    pub articulation_points: Vec<NodeId>,
}

pub fn biconnected_components<W>(graph: &Graph<W>) -> Result<BiconnectedComponents, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let (cuts, edge_blocks) = undirected_low_links(graph, true);
    Ok(BiconnectedComponents {
        edge_blocks,
        articulation_points: cuts.articulation_points,
    })
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EdgeCut {
    /// Edge IDs whose removal separates source and sink.
    pub edges: Vec<EdgeId>,
    pub source_side: Vec<bool>,
}

/// Local edge connectivity with a minimum-cut certificate. Parallel edges
/// count independently. Directed cuts separate outgoing reachability.
pub fn minimum_st_edge_cut<W>(
    graph: &Graph<W>,
    source: NodeId,
    sink: NodeId,
) -> Result<EdgeCut, GraphError> {
    graph.validate_node(source)?;
    graph.validate_node(sink)?;
    if source == sink {
        return Err(GraphError::IdenticalTerminals);
    }
    let undirected = graph.kind() == GraphKind::Undirected;
    let mut builder = GraphBuilder::new(graph.node_count(), GraphKind::Directed);
    for edge in graph.edges() {
        if edge.from == edge.to {
            continue;
        }
        builder.add_edge(edge.from, edge.to, 1).unwrap();
        if undirected {
            builder.add_edge(edge.to, edge.from, 1).unwrap();
        }
    }
    let flow = max_flow(&builder.build(), source, sink)?;
    let edges = graph
        .edges()
        .iter()
        .enumerate()
        .filter_map(|(id, e)| {
            let (from, to) = (flow.source_side[e.from], flow.source_side[e.to]);
            (from != to && (undirected || from)).then_some(id)
        })
        .collect();
    Ok(EdgeCut {
        edges,
        source_side: flow.source_side,
    })
}

/// Global edge connectivity of an undirected graph. Returns zero for
/// disconnected graphs and graphs with fewer than two vertices.
pub fn edge_connectivity<W>(graph: &Graph<W>) -> Result<usize, GraphError> {
    graph.require(GraphKind::Undirected)?;
    if graph.node_count() < 2 || !is_connected(graph)? {
        return Ok(0);
    }
    let cuts = undirected_low_links(graph, false).0;
    if !cuts.bridges.is_empty() {
        return Ok(1);
    }
    let mut minimum_degree = usize::MAX;
    for node in 0..graph.node_count() {
        minimum_degree =
            minimum_degree.min(graph.out(node).iter().filter(|e| e.node != node).count());
    }
    if minimum_degree == 2 {
        return Ok(2);
    }
    let mut builder = GraphBuilder::new(graph.node_count(), GraphKind::Undirected);
    for edge in graph.edges().iter().filter(|edge| edge.from != edge.to) {
        builder.add_edge(edge.from, edge.to, 1).unwrap();
    }
    Ok(global_min_cut(&builder.build())?.unwrap().value as usize)
}
