use crate::{EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EdgeWalk {
    pub nodes: Vec<NodeId>,
    pub edges: Vec<EdgeId>,
}

/// Eulerian path or circuit using every logical edge exactly once. A caller
/// can request a starting vertex. On graphs with edges, isolated vertices
/// are ignored. A singleton graph with no edges has a zero-length path.
pub fn eulerian_path<W>(
    graph: &Graph<W>,
    source: Option<NodeId>,
) -> Result<Option<EdgeWalk>, GraphError> {
    eulerian_walk(graph, source, false)
}

/// Eulerian circuit, using Hierholzer's iterative O(V + E) algorithm.
pub fn eulerian_circuit<W>(
    graph: &Graph<W>,
    source: Option<NodeId>,
) -> Result<Option<EdgeWalk>, GraphError> {
    eulerian_walk(graph, source, true)
}

pub fn is_eulerian<W>(graph: &Graph<W>) -> bool {
    eulerian_walk(graph, None, true).unwrap().is_some()
}

pub fn has_eulerian_path<W>(graph: &Graph<W>) -> bool {
    eulerian_walk(graph, None, false).unwrap().is_some()
}

fn eulerian_walk<W>(
    graph: &Graph<W>,
    requested: Option<NodeId>,
    circuit: bool,
) -> Result<Option<EdgeWalk>, GraphError> {
    if let Some(start) = requested {
        graph.validate_node(start)?;
    }
    let n = graph.node_count();
    if graph.edge_count() == 0 {
        return Ok((n == 1).then(|| EdgeWalk {
            nodes: vec![0],
            edges: Vec::new(),
        }));
    }

    let mut required = Vec::new();
    let mut end_count = 0;
    let mut any = NONE;
    for node in graph.nodes() {
        let outgoing = graph.out(node).len();
        if outgoing > 0 && any == NONE {
            any = node;
        }
        if graph.is_directed() {
            match outgoing as i128 - graph.inc(node).len() as i128 {
                0 => {}
                1 => required.push(node),
                -1 => end_count += 1,
                _ => return Ok(None),
            }
        } else if graph.degree(node).unwrap() % 2 != 0 {
            required.push(node);
        }
    }
    let valid_open = if graph.is_directed() {
        required.len() == 1 && end_count == 1
    } else {
        required.len() == 2
    };
    if !(required.is_empty() && end_count == 0 || !circuit && valid_open) {
        return Ok(None);
    }
    let start = requested.unwrap_or_else(|| required.first().copied().unwrap_or(any));
    if (!required.is_empty() && !required.contains(&start)) || graph.out(start).is_empty() {
        return Ok(None);
    }

    let mut cursor = vec![0; n];
    let mut used = vec![false; graph.edge_count()];
    let mut stack = vec![(start, NONE)];
    let mut nodes = Vec::with_capacity(used.len() + 1);
    let mut edges = Vec::with_capacity(used.len());
    while let Some(&(node, _)) = stack.last() {
        if let Some(arc) = graph.out(node).get(cursor[node]) {
            cursor[node] += 1;
            if !std::mem::replace(&mut used[arc.edge], true) {
                stack.push((arc.node, arc.edge));
            }
        } else {
            let (node, edge) = stack.pop().unwrap();
            nodes.push(node);
            if edge != NONE {
                edges.push(edge);
            }
        }
    }
    if edges.len() != graph.edge_count() {
        return Ok(None);
    }
    nodes.reverse();
    edges.reverse();
    Ok(Some(EdgeWalk { nodes, edges }))
}

/// Fundamental cycle basis of an undirected multigraph. Self-loops and
/// pairs of parallel edges give cycles of length one and two. Runs in
/// O(V + E + total basis length) time.
pub fn cycle_basis<W>(graph: &Graph<W>) -> Result<Vec<EdgeWalk>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let n = graph.node_count();
    let mut parent = vec![NONE; n];
    let mut pred = vec![NONE; n];
    let mut depth = vec![0; n];
    let mut tree = vec![false; graph.edge_count()];
    let mut stack = Vec::new();
    for root in 0..n {
        if parent[root] != NONE {
            continue;
        }
        parent[root] = root;
        stack.push(root);
        while let Some(from) = stack.pop() {
            for arc in graph.out(from) {
                if parent[arc.node] == NONE {
                    tree[arc.edge] = true;
                    parent[arc.node] = from;
                    pred[arc.node] = arc.edge;
                    depth[arc.node] = depth[from] + 1;
                    stack.push(arc.node);
                }
            }
        }
    }

    let mut cycles = Vec::new();
    let mut right_nodes = Vec::new();
    let mut right_edges = Vec::new();
    for (id, edge) in graph.edges().iter().enumerate() {
        if tree[id] {
            continue;
        }
        let (mut a, mut b) = (edge.from, edge.to);
        let mut nodes = vec![a];
        let mut edges = Vec::new();
        right_nodes.clear();
        right_edges.clear();
        while a != b {
            if depth[a] >= depth[b] {
                edges.push(pred[a]);
                a = parent[a];
                nodes.push(a);
            } else {
                right_nodes.push(b);
                right_edges.push(pred[b]);
                b = parent[b];
            }
        }
        nodes.extend(right_nodes.iter().rev());
        edges.extend(right_edges.iter().rev());
        nodes.push(edge.from);
        edges.push(id);
        cycles.push(EdgeWalk { nodes, edges });
    }
    Ok(cycles)
}
