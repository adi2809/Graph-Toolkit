//! Exact triangle counts, clustering coefficients and k-truss decomposition.

use std::collections::VecDeque;

use crate::{EdgeId, Graph, GraphError, GraphKind};

/// Local triangle counts after projecting parallel edges to one simple edge
/// and ignoring self-loops. Uses degree-oriented triangle enumeration.
pub fn triangles<W>(graph: &Graph<W>) -> Result<Vec<u128>, GraphError> {
    let topology = SimpleTopology::new(graph)?;
    Ok(topology.count().1)
}

/// The ratio of closed wedges to all wedges in an undirected graph.
pub fn transitivity<W>(graph: &Graph<W>) -> Result<f64, GraphError> {
    let topology = SimpleTopology::new(graph)?;
    let (_, triangles) = topology.count();
    let wedges: u128 = topology.rows.iter().map(|row| choose_two(row.len())).sum();
    Ok(if wedges == 0 {
        0.0
    } else {
        triangles.into_iter().sum::<u128>() as f64 / wedges as f64
    })
}

/// Unweighted local clustering coefficients in vertex order.
pub fn clustering_coefficients<W>(graph: &Graph<W>) -> Result<Vec<f64>, GraphError> {
    let topology = SimpleTopology::new(graph)?;
    let (_, triangles) = topology.count();
    Ok(topology
        .rows
        .iter()
        .zip(triangles)
        .map(|(row, count)| {
            let possible = choose_two(row.len());
            if possible == 0 {
                0.0
            } else {
                count as f64 / possible as f64
            }
        })
        .collect())
}

/// Maximal undirected subgraph in which every edge belongs to at least k-2
/// triangles. Returns the first original edge ID of each retained simple pair.
pub fn k_truss<W>(graph: &Graph<W>, k: usize) -> Result<Vec<EdgeId>, GraphError> {
    if k < 2 {
        return Err(GraphError::InvalidModelParameters);
    }
    let topology = SimpleTopology::new(graph)?;
    if k == 2 {
        return Ok(topology.original());
    }
    let threshold = k - 2;
    let (mut support, _) = topology.count();
    let mut live = vec![true; topology.edges.len()];
    let mut queue: VecDeque<_> = support
        .iter()
        .enumerate()
        .filter_map(|(edge, &count)| (count < threshold).then_some(edge))
        .collect();
    while let Some(edge) = queue.pop_front() {
        live[edge] = false;
        if support[edge] == 0 {
            continue;
        }
        let (a, b, _) = topology.edges[edge];
        let (mut first, mut second) = (&topology.rows[a], &topology.rows[b]);
        if first.len() > second.len() {
            std::mem::swap(&mut first, &mut second);
        }
        let mut offset = 0;
        for &(node, left) in first {
            if !live[left] {
                continue;
            }
            offset += second[offset..].partition_point(|&(other, _)| other < node);
            let Some(&(other, right)) = second.get(offset) else {
                break;
            };
            if other != node || !live[right] {
                continue;
            }
            for remaining in [left, right] {
                if support[remaining] == threshold {
                    queue.push_back(remaining);
                }
                support[remaining] -= 1;
            }
        }
    }
    let mut result: Vec<_> = topology
        .edges
        .iter()
        .zip(live)
        .filter_map(|(edge, live)| live.then_some(edge.2))
        .collect();
    result.sort_unstable();
    Ok(result)
}

fn choose_two(n: usize) -> u128 {
    n as u128 * n.saturating_sub(1) as u128 / 2
}

struct SimpleTopology {
    edges: Vec<(usize, usize, EdgeId)>,
    rows: Vec<Vec<(usize, usize)>>,
}

impl SimpleTopology {
    fn new<W>(graph: &Graph<W>) -> Result<Self, GraphError> {
        graph.require(GraphKind::Undirected)?;
        let mut edges: Vec<_> = graph
            .edges()
            .iter()
            .enumerate()
            .filter_map(|(id, edge)| {
                (edge.from != edge.to).then_some((
                    edge.from.min(edge.to),
                    edge.from.max(edge.to),
                    id,
                ))
            })
            .collect();
        edges.sort_unstable();
        edges.dedup_by_key(|edge| (edge.0, edge.1));
        let mut rows = vec![Vec::new(); graph.node_count()];
        for (index, &(a, b, _)) in edges.iter().enumerate() {
            rows[a].push((b, index));
            rows[b].push((a, index));
        }
        Ok(Self { edges, rows })
    }

    fn count(&self) -> (Vec<usize>, Vec<u128>) {
        let n = self.rows.len();
        let mut oriented = vec![Vec::new(); n];
        for (index, &(mut a, mut b, _)) in self.edges.iter().enumerate() {
            if (self.rows[a].len(), a) > (self.rows[b].len(), b) {
                std::mem::swap(&mut a, &mut b);
            }
            oriented[a].push((b, index));
        }
        let mut mark = vec![None; n];
        let mut edge_count = vec![0; self.edges.len()];
        let mut vertex_count = vec![0; n];
        for a in 0..n {
            for &(node, edge) in &oriented[a] {
                mark[node] = Some((a, edge));
            }
            for &(b, ab) in &oriented[a] {
                for &(c, bc) in &oriented[b] {
                    if let Some((owner, ac)) = mark[c] {
                        if owner == a {
                            for edge in [ab, bc, ac] {
                                edge_count[edge] += 1;
                            }
                            for node in [a, b, c] {
                                vertex_count[node] += 1;
                            }
                        }
                    }
                }
            }
        }
        (edge_count, vertex_count)
    }

    fn original(&self) -> Vec<EdgeId> {
        let mut result: Vec<_> = self.edges.iter().map(|edge| edge.2).collect();
        result.sort_unstable();
        result
    }
}

/// Estimates the mean local clustering coefficient, counting degree-zero and
/// degree-one vertices as zero. Each trial samples one uniform vertex and two
/// distinct neighbors. After simple projection, a trial takes O(log V) time.
pub fn sampled_average_clustering<W>(
    graph: &Graph<W>,
    samples: usize,
    seed: u64,
) -> Result<f64, GraphError> {
    graph.require(GraphKind::Undirected)?;
    if samples == 0 {
        return Err(GraphError::InvalidModelParameters);
    }
    let topology = SimpleTopology::new(graph)?;
    if topology.rows.is_empty() {
        return Ok(0.0);
    }
    let mut rng = crate::random::Random(seed);
    let mut closed = 0;
    for _ in 0..samples {
        let row = &topology.rows[rng.below(topology.rows.len())];
        if row.len() < 2 {
            continue;
        }
        let i = rng.below(row.len());
        let j = rng.below(row.len() - 1);
        let a = row[i].0;
        let b = row[j + usize::from(j >= i)].0;
        if topology.rows[a]
            .binary_search_by_key(&b, |entry| entry.0)
            .is_ok()
        {
            closed += 1usize;
        }
    }
    Ok(closed as f64 / samples as f64)
}
