//! Deterministic multiscale ball growing for low-average-stretch trees.
use crate::{Graph, GraphError, GraphKind, NONE};

#[derive(Clone, Copy)]
struct Edge {
    u: usize,
    v: usize,
    id: usize,
    class: usize,
}

// Grow BFS balls in the quotient graph of currently eligible length classes.
// Continue while some class has more than a 1/x fraction of its internal
// edge count on the boundary. Accepted layers grow that class geometrically.
// A rejected next layer is discovered, but its adjacency is never scanned.
fn partition(
    n: usize,
    edges: &[Edge],
    phase: usize,
    classes: usize,
    x: u128,
    tree: &mut Vec<usize>,
) -> (Vec<usize>, usize) {
    let mut b = Graph::builder(n, GraphKind::Undirected);
    for (i, e) in edges.iter().enumerate() {
        if e.class <= phase {
            b.add_edge(e.u, e.v, i).unwrap();
        }
    }
    let active = b.build();
    let mut owner = vec![NONE; n];
    let mut depth = vec![0; n];
    let mut discovered = vec![NONE; n];
    let mut incoming = vec![NONE; n];
    let mut layer = Vec::new();
    let mut next = Vec::new();
    let mut inside = vec![0u128; classes];
    let mut crossing = vec![0u128; classes];
    let mut count = 0;
    for root in 0..n {
        if owner[root] != NONE {
            continue;
        }
        owner[root] = count;
        let mut d = 0;
        layer.clear();
        layer.push(root);
        inside.fill(0);
        loop {
            crossing.fill(0);
            next.clear();
            for &u in &layer {
                for a in active.out(u) {
                    let e = edges[active.edges()[a.edge].weight];
                    let v = a.node;
                    if owner[v] == NONE {
                        crossing[e.class] += 1;
                        if discovered[v] != root {
                            discovered[v] = root;
                            incoming[v] = e.id;
                            next.push(v);
                        }
                    } else if owner[v] == count && depth[v] == d && u < v {
                        inside[e.class] += 1;
                    }
                }
            }
            if crossing
                .iter()
                .zip(&inside)
                .all(|(&boundary, &internal)| boundary * x <= internal)
            {
                break;
            }
            d += 1;
            for &v in &next {
                owner[v] = count;
                depth[v] = d;
                tree.push(incoming[v]);
            }
            for (sum, &added) in inside.iter_mut().zip(&crossing) {
                *sum += added;
            }
            std::mem::swap(&mut layer, &mut next);
        }
        count += 1;
    }
    (owner, count)
}

/// Static AKPW-style spanning forest for undirected positive integer edge
/// lengths. Returns original edge IDs. The construction groups comparable
/// lengths and contracts low-radius clusters while bounding the fraction of
/// eligible edges left between clusters in each class. The objective is low
/// average tree-distance / edge-length stretch, including parallel edges.
/// Zero lengths are rejected. Disconnected components are processed together.
pub fn low_stretch_spanning_forest(graph: &Graph<u64>) -> Result<Vec<usize>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    if graph.edges().iter().any(|e| e.weight == 0) {
        return Err(GraphError::InvalidModelParameters);
    }
    Ok(build(
        graph.node_count(),
        graph
            .edges()
            .iter()
            .enumerate()
            .map(|(id, e)| (e.from, e.to, e.weight, id)),
    ))
}

pub(crate) fn build(
    n: usize,
    input: impl ExactSizeIterator<Item = (usize, usize, u64, usize)>,
) -> Vec<usize> {
    // Using max(n,m) also handles multigraphs with unbounded multiplicities.
    let log = (n.max(input.len()).max(3) as f64).ln();
    let x = (log * log.ln().max(1.)).sqrt().exp().ceil().max(2.) as u128;
    let rho = (3. * log / (x as f64).ln()).ceil();
    let base = x * (9. * rho * log).ceil() as u128;
    let mut ends = vec![base];
    while *ends.last().unwrap() <= u128::from(u64::MAX) {
        ends.push(ends.last().unwrap() * base);
    }
    let mut edges: Vec<_> = input
        .filter(|(u, v, _, _)| u != v)
        .map(|(u, v, length, id)| Edge {
            u,
            v,
            id,
            class: ends.partition_point(|&end| end <= u128::from(length)),
        })
        .collect();
    let classes = ends.len();
    let mut tree = Vec::new();
    let mut nodes = n;
    let mut phase = 0;
    while !edges.is_empty() {
        let (part, count) = partition(nodes, &edges, phase, classes, x, &mut tree);
        edges.retain_mut(|e| {
            e.u = part[e.u];
            e.v = part[e.v];
            e.u != e.v
        });
        nodes = count;
        phase += 1;
    }
    tree
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn each_active_class_retains_only_a_small_boundary() {
        let n = 400;
        let mut edges = Vec::new();
        for u in 0..n {
            for (hop, class) in [(1, 0), (1, 1), (2, 0), (3, 1), (5, 2)] {
                if u + hop < n {
                    edges.push(Edge {
                        u,
                        v: u + hop,
                        id: edges.len(),
                        class,
                    });
                }
            }
        }
        let mut tree = Vec::new();
        let (part, count) = partition(n, &edges, 2, 3, 3, &mut tree);
        assert!(count > 1 && count < n);
        assert_eq!(tree.len() + count, n);
        let mut internal = [0; 3];
        let mut cut = [0; 3];
        for e in edges {
            if part[e.u] == part[e.v] {
                internal[e.class] += 1;
            } else {
                cut[e.class] += 1;
            }
        }
        for c in 0..3 {
            assert!(3 * cut[c] <= internal[c]);
        }
    }
}
