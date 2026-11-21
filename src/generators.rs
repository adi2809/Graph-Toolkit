//! Deterministic and seeded generators for dense-ID CSR graphs.
//! Topology generators emit unit payloads. Graph transforms can retain edge payloads.

mod algebraic;
mod geometric;
mod lattice;
mod named;
mod random;
pub use algebraic::{
    chordal_cycle_graph, hkn_harary_graph, kneser_graph, paley_graph, random_tournament,
};
pub use geometric::{GeometricGraph, geometric_edges, random_geometric_graph};
pub use lattice::{hexagonal_lattice_graph, triangular_lattice_graph};
pub use named::{
    bull_graph, chvatal_graph, diamond_graph, dorogovtsev_goltsev_mendes_graph, frucht_graph,
    heawood_graph, hoffman_singleton_graph, house_graph, house_x_graph, icosahedral_graph,
    krackhardt_kite_graph, lcf_graph, margulis_gabber_galil_graph, pappus_graph,
    sedgewick_maze_graph, shrikhande_graph, truncated_cube_graph, truncated_tetrahedron_graph,
    tutte_graph,
};
pub use random::{
    barabasi_albert_graph, bipartite_gnp_random_graph, configuration_model,
    connected_watts_strogatz_graph, degree_sequence_tree, dense_gnm_random_graph,
    directed_configuration_model, erdos_renyi_graph, expected_degree_graph, fast_gnp_random_graph,
    gn_graph, gnc_graph, gnm_random_graph, gnp_random_graph, gnr_graph, havel_hakimi_graph,
    newman_watts_strogatz_graph, planted_partition_graph, random_clustered_graph,
    random_intersection_graph, random_k_lift, random_labeled_tree, random_partition_graph,
    random_regular_bipartite_multigraph, random_regular_graph, stochastic_block_model,
    watts_strogatz_graph,
};

use crate::{Graph, GraphBuilder, GraphError, GraphKind};

pub fn empty_graph(n: usize, kind: GraphKind) -> Graph<()> {
    GraphBuilder::new(n, kind).build()
}

pub fn null_graph() -> Graph<()> {
    empty_graph(0, GraphKind::Undirected)
}

pub fn trivial_graph() -> Graph<()> {
    empty_graph(1, GraphKind::Undirected)
}

pub fn path_graph(n: usize, kind: GraphKind) -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(n, kind, n.saturating_sub(1));
    for node in 1..n {
        add(&mut builder, node - 1, node);
    }
    builder.build()
}

pub fn cycle_graph(n: usize, kind: GraphKind) -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(n, kind, n);
    for node in 1..n {
        add(&mut builder, node - 1, node);
    }
    if n > 0 && (n != 2 || kind == GraphKind::Directed) {
        add(&mut builder, n - 1, 0);
    }
    builder.build()
}

pub fn complete_graph(n: usize, kind: GraphKind) -> Result<Graph<()>, GraphError> {
    let count = candidate_count(n, kind)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, kind, count);
    for from in 0..n {
        let begin = if kind == GraphKind::Directed {
            0
        } else {
            from + 1
        };
        for to in begin..n {
            if to != from {
                add(&mut builder, from, to);
            }
        }
    }
    Ok(builder.build())
}

/// Vertices of each part form one contiguous block of IDs.
pub fn complete_multipartite_graph(sizes: &[usize]) -> Result<Graph<()>, GraphError> {
    let (mut n, mut edges) = (0usize, 0usize);
    for &size in sizes {
        edges = edges
            .checked_add(n.checked_mul(size).ok_or(GraphError::GraphSizeOverflow)?)
            .ok_or(GraphError::GraphSizeOverflow)?;
        n = n.checked_add(size).ok_or(GraphError::GraphSizeOverflow)?;
    }
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, edges);
    let mut first = 0;
    for &size in sizes {
        for to in first..first + size {
            for from in 0..first {
                add(&mut builder, from, to);
            }
        }
        first += size;
    }
    Ok(builder.build())
}

pub fn complete_bipartite_graph(left: usize, right: usize) -> Result<Graph<()>, GraphError> {
    complete_multipartite_graph(&[left, right])
}

pub fn turan_graph(n: usize, parts: usize) -> Result<Graph<()>, GraphError> {
    if parts == 0 || parts > n {
        return Err(GraphError::InvalidModelParameters);
    }
    let sizes: Vec<_> = (0..parts)
        .map(|i| n / parts + usize::from(i < n % parts))
        .collect();
    complete_multipartite_graph(&sizes)
}

pub fn star_graph(leaves: usize, kind: GraphKind) -> Result<Graph<()>, GraphError> {
    let n = leaves.checked_add(1).ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, kind, leaves);
    for to in 1..n {
        add(&mut builder, 0, to);
    }
    Ok(builder.build())
}

pub fn wheel_graph(n: usize) -> Graph<()> {
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    for to in 1..n {
        add(&mut builder, 0, to);
        if to > 1 {
            add(&mut builder, to - 1, to);
        }
    }
    if n > 3 {
        add(&mut builder, n - 1, 1);
    }
    builder.build()
}

pub fn ladder_graph(rungs: usize) -> Result<Graph<()>, GraphError> {
    ladder(rungs, false)
}

pub fn circular_ladder_graph(rungs: usize) -> Result<Graph<()>, GraphError> {
    ladder(rungs, true)
}

fn ladder(rungs: usize, circular: bool) -> Result<Graph<()>, GraphError> {
    let n = rungs.checked_mul(2).ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    for i in 0..rungs {
        add(&mut builder, i, i + rungs);
        if i > 0 {
            add(&mut builder, i - 1, i);
            add(&mut builder, i + rungs - 1, i + rungs);
        }
    }
    if circular && rungs > 2 {
        add(&mut builder, 0, rungs - 1);
        add(&mut builder, rungs, n - 1);
    }
    Ok(builder.build())
}

pub fn lollipop_graph(clique: usize, tail: usize) -> Result<Graph<()>, GraphError> {
    if clique < 2 {
        return Err(GraphError::InvalidModelParameters);
    }
    let n = clique
        .checked_add(tail)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    add_clique(&mut builder, 0..clique);
    for node in clique..n {
        add(&mut builder, node - 1, node);
    }
    Ok(builder.build())
}

pub fn barbell_graph(clique: usize, middle: usize) -> Result<Graph<()>, GraphError> {
    if clique < 2 {
        return Err(GraphError::InvalidModelParameters);
    }
    let n = clique
        .checked_mul(2)
        .and_then(|v| v.checked_add(middle))
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    add_clique(&mut builder, 0..clique);
    add_clique(&mut builder, clique + middle..n);
    for node in clique..=clique + middle {
        add(&mut builder, node - 1, node);
    }
    Ok(builder.build())
}

pub fn tadpole_graph(cycle: usize, tail: usize) -> Result<Graph<()>, GraphError> {
    if cycle < 3 {
        return Err(GraphError::InvalidModelParameters);
    }
    let n = cycle
        .checked_add(tail)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, n);
    for node in 1..n {
        add(&mut builder, node - 1, node);
    }
    add(&mut builder, cycle - 1, 0);
    Ok(builder.build())
}

fn add_clique(builder: &mut GraphBuilder<()>, nodes: std::ops::Range<usize>) {
    for a in nodes.clone() {
        for b in a + 1..nodes.end {
            add(builder, a, b);
        }
    }
}

/// Breadth-first labeling of a complete r-ary tree, truncated at n vertices.
pub fn full_rary_tree(branching: usize, n: usize) -> Result<Graph<()>, GraphError> {
    if branching == 0 && n > 1 {
        return Err(GraphError::InvalidModelParameters);
    }
    let mut builder =
        GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, n.saturating_sub(1));
    for node in 1..n {
        add(&mut builder, (node - 1) / branching, node);
    }
    Ok(builder.build())
}

pub fn balanced_tree(branching: usize, height: usize) -> Result<Graph<()>, GraphError> {
    let n = match branching {
        0 => 1,
        1 => height.checked_add(1).ok_or(GraphError::GraphSizeOverflow)?,
        _ => (0..height).try_fold(1usize, |n, _| {
            n.checked_mul(branching)
                .and_then(|v| v.checked_add(1))
                .ok_or(GraphError::GraphSizeOverflow)
        })?,
    };
    full_rary_tree(branching, n)
}

pub fn binomial_tree(order: u32) -> Result<Graph<()>, GraphError> {
    let n = 1usize
        .checked_shl(order)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, n - 1);
    for node in 1..n {
        let high_bit = 1usize << node.ilog2();
        add(&mut builder, node ^ high_bit, node);
    }
    Ok(builder.build())
}

/// Cartesian product of paths (or cycles along periodic axes). Axis zero
/// changes fastest in the dense vertex IDs.
pub fn grid_graph(dimensions: &[usize], periodic: &[bool]) -> Result<Graph<()>, GraphError> {
    if periodic.len() != dimensions.len() {
        return Err(GraphError::InvalidModelParameters);
    }
    if dimensions.contains(&0) || dimensions.is_empty() {
        return Ok(null_graph());
    }
    let mut strides = Vec::with_capacity(dimensions.len());
    let mut n = 1usize;
    for &size in dimensions {
        strides.push(n);
        n = n.checked_mul(size).ok_or(GraphError::GraphSizeOverflow)?;
    }
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    for node in 0..n {
        for ((&size, &stride), &wrap) in dimensions.iter().zip(&strides).zip(periodic) {
            let coordinate = node / stride % size;
            if coordinate + 1 < size {
                add(&mut builder, node, node + stride);
            } else if wrap && size > 2 {
                add(&mut builder, node, node - (size - 1) * stride);
            }
        }
    }
    Ok(builder.build())
}

pub fn grid_2d_graph(
    width: usize,
    height: usize,
    periodic: [bool; 2],
) -> Result<Graph<()>, GraphError> {
    grid_graph(&[width, height], &periodic)
}

pub fn hypercube_graph(dimensions: u32) -> Result<Graph<()>, GraphError> {
    let n = 1usize
        .checked_shl(dimensions)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let m = (n / 2)
        .checked_mul(dimensions as usize)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, m);
    for node in 0..n {
        let mut unused = (n - 1) ^ node;
        while unused != 0 {
            let bit = 1 << unused.trailing_zeros();
            unused &= unused - 1;
            add(&mut builder, node, node | bit);
        }
    }
    Ok(builder.build())
}

pub fn circulant_graph(n: usize, offsets: &[usize]) -> Graph<()> {
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    if n == 0 {
        return builder.build();
    }
    let mut offsets: Vec<_> = offsets
        .iter()
        .map(|&i| {
            let i = i % n;
            i.min(n - i)
        })
        .collect();
    offsets.sort_unstable();
    offsets.dedup();
    for step in offsets {
        let limit = if step != 0 && step * 2 == n { step } else { n };
        for node in 0..limit {
            add(&mut builder, node, advance(node, step, n));
        }
    }
    builder.build()
}

pub fn generalized_petersen_graph(n: usize, step: usize) -> Result<Graph<()>, GraphError> {
    if n < 3 || step == 0 || step >= n.div_ceil(2) {
        return Err(GraphError::InvalidModelParameters);
    }
    let nodes = n.checked_mul(2).ok_or(GraphError::GraphSizeOverflow)?;
    let edges = n.checked_mul(3).ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(nodes, GraphKind::Undirected, edges);
    for i in 0..n {
        add(&mut builder, i, advance(i, 1, n));
        add(&mut builder, i, i + n);
        add(&mut builder, i + n, advance(i, step, n) + n);
    }
    Ok(builder.build())
}

pub fn petersen_graph() -> Graph<()> {
    generalized_petersen_graph(5, 2).unwrap()
}
pub fn dodecahedral_graph() -> Graph<()> {
    generalized_petersen_graph(10, 2).unwrap()
}
pub fn desargues_graph() -> Graph<()> {
    generalized_petersen_graph(10, 3).unwrap()
}
pub fn moebius_kantor_graph() -> Graph<()> {
    generalized_petersen_graph(8, 3).unwrap()
}
pub fn cubical_graph() -> Graph<()> {
    hypercube_graph(3).unwrap()
}
pub fn tetrahedral_graph() -> Graph<()> {
    complete_graph(4, GraphKind::Undirected).unwrap()
}
pub fn octahedral_graph() -> Graph<()> {
    complete_multipartite_graph(&[2, 2, 2]).unwrap()
}

fn add(builder: &mut GraphBuilder<()>, from: usize, to: usize) {
    builder.add_edge(from, to, ()).unwrap();
}

fn advance(node: usize, offset: usize, n: usize) -> usize {
    if node >= n - offset {
        node - (n - offset)
    } else {
        node + offset
    }
}

fn candidate_count(n: usize, kind: GraphKind) -> Result<usize, GraphError> {
    let mut count = n as u128 * n.saturating_sub(1) as u128;
    if kind == GraphKind::Undirected {
        count /= 2;
    }
    usize::try_from(count).map_err(|_| GraphError::GraphSizeOverflow)
}

/// Disjoint union of equally sized cliques.
pub fn caveman_graph(cliques: usize, size: usize) -> Result<Graph<()>, GraphError> {
    let n = cliques
        .checked_mul(size)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let m = candidate_count(size, GraphKind::Undirected)?
        .checked_mul(cliques)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, m);
    for i in 0..cliques {
        add_clique(&mut builder, i * size..(i + 1) * size);
    }
    Ok(builder.build())
}

/// Clique interiors with one edge from each clique to the next clique in a ring.
pub fn ring_of_cliques(cliques: usize, size: usize) -> Result<Graph<()>, GraphError> {
    if cliques < 2 || size < 2 {
        return Err(GraphError::InvalidModelParameters);
    }
    let n = cliques
        .checked_mul(size)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    for i in 0..cliques {
        add_clique(&mut builder, i * size..(i + 1) * size);
        add(
            &mut builder,
            i * size + size - 1,
            ((i + 1) % cliques) * size,
        );
    }
    Ok(builder.build())
}

/// Union of cliques sharing exactly one center vertex.
pub fn windmill_graph(cliques: usize, size: usize) -> Result<Graph<()>, GraphError> {
    if size < 2 || cliques < 2 {
        return Err(GraphError::InvalidModelParameters);
    }
    let stride = size - 1;
    let n = cliques
        .checked_mul(stride)
        .and_then(|x| x.checked_add(1))
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    for i in 0..cliques {
        let members = 1 + i * stride..1 + (i + 1) * stride;
        for node in members.clone() {
            add(&mut builder, 0, node);
        }
        add_clique(&mut builder, members);
    }
    Ok(builder.build())
}

/// Vertices of the line graph correspond to original logical edge IDs.
/// Undirected edges are adjacent when they share an endpoint. Directed
/// edges are adjacent when they can be traversed consecutively.
pub fn line_graph<W>(graph: &Graph<W>) -> Graph<()> {
    let mut builder = GraphBuilder::new(graph.edge_count(), graph.kind());
    for node in graph.nodes() {
        let outgoing = graph.out(node);
        if graph.is_directed() {
            for a in graph.inc(node) {
                for b in outgoing {
                    add(&mut builder, a.edge, b.edge);
                }
            }
        } else {
            for (position, a) in outgoing.iter().enumerate() {
                for b in &outgoing[position + 1..] {
                    // A pair of parallel edges shares both endpoints.
                    if a.node != b.node || node <= a.node {
                        add(&mut builder, a.edge, b.edge);
                    }
                }
            }
        }
    }
    builder.build()
}

/// Hop-radius induced neighborhood using outgoing reachability.
pub fn ego_graph<W: Clone>(
    graph: &Graph<W>,
    center: usize,
    radius: usize,
    include_center: bool,
) -> Result<crate::functions::MaterializedGraph<W>, GraphError> {
    graph.validate_node(center)?;
    let mut seen = vec![false; graph.node_count()];
    seen[center] = true;
    let mut order = vec![center];
    let mut start = 0;
    for _ in 0..radius {
        let end = order.len();
        if start == end {
            break;
        }
        for i in start..end {
            for arc in graph.out(order[i]) {
                if !std::mem::replace(&mut seen[arc.node], true) {
                    order.push(arc.node);
                }
            }
        }
        start = end;
    }
    let selected = order
        .into_iter()
        .filter(|&node| include_center || node != center);
    Ok(crate::functions::induced_subgraph(graph, selected)?.materialize())
}
