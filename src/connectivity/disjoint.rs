use super::EdgeWalk;
use crate::flow::max_flow;
use crate::{EdgeId, Graph, GraphBuilder, GraphError, GraphKind, NONE, NodeId};

/// Maximum family of edge-disjoint s–t paths. Each logical edge can appear
/// on at most one path, including when traversing an undirected graph.
/// Parallel edges are distinct.
pub fn edge_disjoint_paths<W>(
    graph: &Graph<W>,
    source: NodeId,
    sink: NodeId,
) -> Result<Vec<EdgeWalk>, GraphError> {
    terminals(graph, source, sink)?;
    let (network, mapping) = unit_network(graph, false, 1);
    disjoint_flow(
        &network,
        &mapping,
        source,
        sink,
        false,
        !graph.is_directed(),
    )
}

/// Maximum family of internally vertex-disjoint s–t paths via vertex
/// splitting. Distinct parallel direct edges can each form a path.
pub fn node_disjoint_paths<W>(
    graph: &Graph<W>,
    source: NodeId,
    sink: NodeId,
) -> Result<Vec<EdgeWalk>, GraphError> {
    terminals(graph, source, sink)?;
    let (network, mapping) = unit_network(graph, true, 1);
    disjoint_flow(
        &network,
        &mapping,
        2 * source + 1,
        2 * sink,
        true,
        !graph.is_directed(),
    )
}

/// A minimum set of internal vertices whose removal makes the sink
/// unreachable from the source. None means the terminals have a direct arc.
/// A disconnected pair has an empty cut.
pub fn minimum_st_vertex_cut<W>(
    graph: &Graph<W>,
    source: NodeId,
    sink: NodeId,
) -> Result<Option<Vec<NodeId>>, GraphError> {
    terminals(graph, source, sink)?;
    if graph.out(source).iter().any(|arc| arc.node == sink) {
        return Ok(None);
    }
    let n = graph.node_count();
    let (network, _) = unit_network(graph, true, n as u64);
    let cut = max_flow(&network, 2 * source + 1, 2 * sink)?;
    Ok(Some(
        graph
            .nodes()
            .filter(|&node| {
                node != source
                    && node != sink
                    && cut.source_side[2 * node]
                    && !cut.source_side[2 * node + 1]
            })
            .collect(),
    ))
}

fn terminals<W>(graph: &Graph<W>, source: NodeId, sink: NodeId) -> Result<(), GraphError> {
    graph.validate_node(source)?;
    graph.validate_node(sink)?;
    if source == sink {
        return Err(GraphError::IdenticalTerminals);
    }
    Ok(())
}

fn unit_network<W>(graph: &Graph<W>, split: bool, capacity: u64) -> (Graph<u64>, Vec<EdgeId>) {
    let n = graph.node_count();
    let mut builder = GraphBuilder::new(if split { 2 * n } else { n }, GraphKind::Directed);
    let mut mapping = Vec::new();
    for (id, edge) in graph.edges().iter().enumerate() {
        if edge.from == edge.to {
            continue;
        }
        let mut add = |a, b| {
            let (from, to) = if split { (2 * a + 1, 2 * b) } else { (a, b) };
            builder.add_edge(from, to, capacity).unwrap();
            mapping.push(id);
        };
        add(edge.from, edge.to);
        if !graph.is_directed() {
            add(edge.to, edge.from);
        }
    }
    if split {
        for node in graph.nodes() {
            builder.add_edge(2 * node, 2 * node + 1, 1).unwrap();
            mapping.push(NONE);
        }
    }
    (builder.build(), mapping)
}

fn disjoint_flow(
    network: &Graph<u64>,
    mapping: &[EdgeId],
    source: NodeId,
    sink: NodeId,
    split: bool,
    undirected: bool,
) -> Result<Vec<EdgeWalk>, GraphError> {
    let mut flow = max_flow(network, source, sink)?;
    if undirected && !split {
        // Paired arcs represent one logical undirected capacity. Opposing
        // flows cancel without changing balances or the total value.
        for pair in flow.edge_flows.chunks_exact_mut(2) {
            if pair[0] != 0 && pair[1] != 0 {
                pair.fill(0);
            }
        }
    }
    let original_node = |node| if split { node / 2 } else { node };
    let mut result = Vec::with_capacity(flow.value as usize);
    let mut cursor = vec![0; network.node_count()];
    let mut position = vec![NONE; network.node_count()];
    let mut nodes = vec![source];
    let mut edges: Vec<EdgeId> = Vec::new();
    position[source] = 0;
    while result.len() < flow.value as usize {
        let node = *nodes.last().unwrap();
        if node == sink {
            let mut path = EdgeWalk {
                nodes: vec![original_node(source)],
                edges: Vec::new(),
            };
            for &edge in &edges {
                if mapping[edge] != NONE {
                    path.edges.push(mapping[edge]);
                    path.nodes.push(original_node(network.edges()[edge].to));
                }
            }
            for v in nodes.drain(1..) {
                position[v] = NONE;
            }
            edges.clear();
            result.push(path);
            continue;
        }
        let arcs = network.out(node);
        let next = &mut cursor[node];
        while flow.edge_flows[arcs[*next].edge] == 0 {
            *next += 1;
        }
        let arc = arcs[*next];
        *next += 1;
        flow.edge_flows[arc.edge] = 0;
        let existing = position[arc.node];
        if existing == NONE {
            position[arc.node] = nodes.len();
            nodes.push(arc.node);
            edges.push(arc.edge);
        } else {
            // Delete the cycle from the unfinished walk and keep its prefix.
            for v in nodes.drain(existing + 1..) {
                position[v] = NONE;
            }
            edges.truncate(existing);
        }
    }
    Ok(result)
}
