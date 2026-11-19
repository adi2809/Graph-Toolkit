use crate::connectivity::{EdgeWalk, biconnected_components};
use crate::{EdgeId, Graph, GraphBuilder, GraphKind, NONE, NodeId};

use super::{PlanarEmbedding, embed_simple, is_simple_planar, simple_projection};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum KuratowskiKind {
    K5,
    K33,
}

/// A subdivision of K₅ or K₃,₃ in the underlying undirected graph.
/// Each path represents one edge of the obstruction and uses original
/// vertex and edge IDs. Path interiors are pairwise disjoint and contain
/// no branch vertices. For K₃,₃ the first and last three branch vertices
/// give the two sides of the bipartition.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KuratowskiWitness {
    pub kind: KuratowskiKind,
    pub branch_vertices: Vec<NodeId>,
    pub paths: Vec<EdgeWalk>,
}

impl KuratowskiWitness {
    pub fn edge_ids(&self) -> impl Iterator<Item = EdgeId> + '_ {
        self.paths
            .iter()
            .flat_map(|path| path.edges.iter().copied())
    }
}

/// Returns an embedding on planar inputs and a Kuratowski subdivision on
/// nonplanar inputs. Direction, loops, and parallel copies have the same
/// projection semantics as `check_planarity`.
pub fn check_planarity_with_witness<W>(
    graph: &Graph<W>,
) -> Result<PlanarEmbedding, KuratowskiWitness> {
    let simple = simple_projection(graph, true).unwrap();
    embed_simple(&simple).ok_or_else(|| extract(graph, &simple))
}

/// Finds a Kuratowski subdivision, or returns None for a planar graph.
/// Extraction isolates a nonplanar biconnected block, suppresses degree-two
/// chains, and minimizes the remaining edges by batched planarity probes.
/// For a compressed block with B vertices and C edges the extraction bound
/// is O(V + E + C(B + C)) time and O(V + E) working space.
pub fn kuratowski_witness<W>(graph: &Graph<W>) -> Option<KuratowskiWitness> {
    if graph.node_count() < 5 || graph.edge_count() < 9 {
        return None;
    }
    let simple = simple_projection(graph, true).unwrap();
    (!is_simple_planar(&simple)).then(|| extract(graph, &simple))
}

fn extract<W>(original: &Graph<W>, simple: &Graph<EdgeId>) -> KuratowskiWitness {
    let mut blocks = biconnected_components(simple).unwrap().edge_blocks;
    blocks.retain(|block| block.len() >= 9);
    // Prefer small blocks, which can often certify the entire graph.
    if let Some((i, _)) = blocks
        .iter()
        .enumerate()
        .min_by_key(|(_, block)| block.len())
    {
        blocks.swap(0, i);
    }
    let mut projection = Compact::new(simple.node_count());
    for edges in blocks {
        let block = projection.build(edges.into_iter().map(|id| {
            let edge = simple.edges()[id];
            (edge.from, edge.to, edge.weight)
        }));
        let kernel = Kernel::new(&block, &projection.nodes);
        let mut probe = Compact::new(kernel.nodes.len());
        if !kernel.test(&vec![true; kernel.chains.len()], &mut probe) {
            return kernel.witness(original, &kernel.minimize(&mut probe));
        }
    }
    unreachable!("a nonplanar simple graph has a nonplanar block")
}

/// Lazily remaps only the vertices touched by a probe.
struct Compact {
    local: Vec<usize>,
    nodes: Vec<NodeId>,
    edges: Vec<(usize, usize, usize)>,
}

impl Compact {
    fn new(n: usize) -> Self {
        Self {
            local: vec![NONE; n],
            nodes: Vec::new(),
            edges: Vec::new(),
        }
    }

    fn build(&mut self, edges: impl IntoIterator<Item = (usize, usize, usize)>) -> Graph<usize> {
        for node in self.nodes.drain(..) {
            self.local[node] = NONE;
        }
        self.edges.clear();
        for (u, v, data) in edges {
            for node in [u, v] {
                if self.local[node] == NONE {
                    self.local[node] = self.nodes.len();
                    self.nodes.push(node);
                }
            }
            self.edges.push((self.local[u], self.local[v], data));
        }
        let mut builder = GraphBuilder::with_edge_capacity(
            self.nodes.len(),
            GraphKind::Undirected,
            self.edges.len(),
        );
        for &(u, v, data) in &self.edges {
            builder.add_edge(u, v, data).unwrap();
        }
        builder.build()
    }
}

struct Chain {
    from: usize,
    to: usize,
    begin: usize,
    end: usize,
}

struct Kernel {
    nodes: Vec<NodeId>,
    chains: Vec<Chain>,
    edge_ids: Vec<EdgeId>,
}

impl Kernel {
    fn new(block: &Graph<EdgeId>, original_nodes: &[NodeId]) -> Self {
        let mut local = vec![NONE; block.node_count()];
        let mut nodes = Vec::new();
        for u in block.nodes().filter(|&u| block.out(u).len() != 2) {
            local[u] = nodes.len();
            nodes.push(original_nodes[u]);
        }
        let mut result = Self {
            nodes,
            chains: Vec::new(),
            edge_ids: Vec::new(),
        };
        let mut used = vec![false; block.edge_count()];
        let mut path = Vec::new();
        let mut last_source = vec![NONE; result.nodes.len()];
        // Trace maximal degree-two chains. Retaining one chain per parallel
        // pair suffices to preserve planarity, and the retained chains lift
        // to paths with disjoint interiors in the input graph.
        for start in block.nodes().filter(|&u| local[u] != NONE) {
            let from = local[start];
            for &arc in block.out(start) {
                if used[arc.edge] {
                    continue;
                }
                path.clear();
                let mut v = arc.node;
                let mut id = arc.edge;
                loop {
                    used[id] = true;
                    path.push(block.edges()[id].weight);
                    if local[v] != NONE {
                        break;
                    }
                    let row = block.out(v);
                    let next = row[usize::from(row[0].edge == id)];
                    v = next.node;
                    id = next.edge;
                }
                let to = local[v];
                if to == from || last_source[to] == from {
                    continue;
                }
                last_source[to] = from;
                let begin = result.edge_ids.len();
                result.edge_ids.extend_from_slice(&path);
                result.chains.push(Chain {
                    from,
                    to,
                    begin,
                    end: result.edge_ids.len(),
                });
            }
        }
        result
    }

    fn test(&self, alive: &[bool], probe: &mut Compact) -> bool {
        let graph = probe.build(
            self.chains
                .iter()
                .zip(alive)
                .enumerate()
                .filter_map(|(id, (chain, &keep))| keep.then_some((chain.from, chain.to, id))),
        );
        is_simple_planar(&graph)
    }

    fn minimize(&self, probe: &mut Compact) -> Vec<bool> {
        let mut alive = vec![true; self.chains.len()];
        let mut pending = vec![(0, alive.len())];
        while let Some((first, end)) = pending.pop() {
            alive[first..end].fill(false);
            if self.test(&alive, probe) {
                // A necessary singleton stays necessary after later deletions.
                alive[first..end].fill(true);
                if end - first > 1 {
                    let middle = first + (end - first) / 2;
                    pending.push((middle, end));
                    pending.push((first, middle));
                }
            }
        }
        alive
    }

    fn witness<W>(&self, original: &Graph<W>, alive: &[bool]) -> KuratowskiWitness {
        let mut builder = GraphBuilder::new(self.nodes.len(), GraphKind::Undirected);
        for (id, chain) in self.chains.iter().enumerate() {
            if alive[id] {
                builder.add_edge(chain.from, chain.to, id).unwrap();
            }
        }
        let graph = builder.build();
        let mut branch_vertices: Vec<_> = graph
            .nodes()
            .filter(|&u| graph.out(u).len() >= 3)
            .map(|u| self.nodes[u])
            .collect();
        branch_vertices.sort_unstable();
        let kind = if branch_vertices.len() == 5 {
            KuratowskiKind::K5
        } else {
            KuratowskiKind::K33
        };
        let mut used = vec![false; graph.edge_count()];
        let mut paths = Vec::new();
        for start in graph.nodes().filter(|&u| graph.out(u).len() >= 3) {
            for &first in graph.out(start) {
                if used[first.edge] {
                    continue;
                }
                let mut path = EdgeWalk {
                    nodes: vec![self.nodes[start]],
                    edges: Vec::new(),
                };
                let mut from = start;
                let mut next = first;
                loop {
                    used[next.edge] = true;
                    let chain = &self.chains[graph.edges()[next.edge].weight];
                    let ids = &self.edge_ids[chain.begin..chain.end];
                    if from == chain.from {
                        for &id in ids {
                            append(original, &mut path, id);
                        }
                    } else {
                        for &id in ids.iter().rev() {
                            append(original, &mut path, id);
                        }
                    }
                    from = next.node;
                    let row = graph.out(from);
                    if row.len() != 2 {
                        break;
                    }
                    next = row[usize::from(row[0].edge == next.edge)];
                }
                paths.push(path);
            }
        }
        if kind == KuratowskiKind::K33 {
            let first = branch_vertices[0];
            let right: Vec<_> = paths
                .iter()
                .filter_map(|path| {
                    let (u, v) = (path.nodes[0], *path.nodes.last().unwrap());
                    if u == first {
                        Some(v)
                    } else if v == first {
                        Some(u)
                    } else {
                        None
                    }
                })
                .collect();
            branch_vertices.sort_by_key(|v| right.contains(v));
            debug_assert_eq!(branch_vertices.len(), 6);
            debug_assert_eq!(paths.len(), 9);
        } else {
            debug_assert_eq!(paths.len(), 10);
        }
        KuratowskiWitness {
            kind,
            branch_vertices,
            paths,
        }
    }
}

fn append<W>(graph: &Graph<W>, walk: &mut EdgeWalk, id: EdgeId) {
    let edge = &graph.edges()[id];
    let from = *walk.nodes.last().unwrap();
    debug_assert!(from == edge.from || from == edge.to);
    walk.nodes.push(if edge.from == from {
        edge.to
    } else {
        edge.from
    });
    walk.edges.push(id);
}
