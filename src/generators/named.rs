//! Named graphs and explicit algebraic constructions.
//!
//! Named graph IDs follow the constructions below. See docs/algorithm-survey.md.

use std::collections::HashSet;

use super::{add, cubical_graph, null_graph, tetrahedral_graph};
use crate::planarity::check_planarity;
use crate::{Graph, GraphBuilder, GraphError, GraphKind};

/// A cycle with chords encoded by LCF shifts. Repeated undirected pairs
/// are coalesced. A zero shift adds a self-loop.
pub fn lcf_graph(n: usize, shifts: &[isize], repeats: usize) -> Result<Graph<()>, GraphError> {
    if n == 0 {
        return Ok(null_graph());
    }
    let chords = repeats
        .checked_mul(shifts.len())
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    let mut seen = HashSet::new();
    let mut emit = |u, v| {
        let pair = if u <= v { (u, v) } else { (v, u) };
        if seen.insert(pair) {
            add(&mut builder, u, v);
        }
    };
    for u in 0..n {
        emit(u, (u + 1) % n);
    }
    for i in 0..chords {
        let u = i % n;
        let v = (u as i128 + shifts[i % shifts.len()] as i128).rem_euclid(n as i128) as usize;
        emit(u, v);
    }
    Ok(builder.build())
}

pub fn heawood_graph() -> Graph<()> {
    lcf_graph(14, &[5, -5], 7).unwrap()
}
pub fn pappus_graph() -> Graph<()> {
    lcf_graph(18, &[5, 7, -7, 7, -7, -5], 3).unwrap()
}

// Build small fixtures from paths in their defining drawings.
fn paths(n: usize, walks: &[&[usize]]) -> Graph<()> {
    let mut b = GraphBuilder::new(n, GraphKind::Undirected);
    for path in walks {
        for step in path.windows(2) {
            add(&mut b, step[0], step[1]);
        }
    }
    b.build()
}

pub fn bull_graph() -> Graph<()> {
    paths(5, &[&[3, 1, 0, 2, 4], &[1, 2]])
}
pub fn diamond_graph() -> Graph<()> {
    paths(4, &[&[0, 1, 3, 2, 0], &[1, 2]])
}
pub fn house_graph() -> Graph<()> {
    house(false)
}
pub fn house_x_graph() -> Graph<()> {
    house(true)
}
fn house(crossed: bool) -> Graph<()> {
    let mut walks = vec![&[0, 1, 3, 2, 0][..], &[2, 4, 3]];
    if crossed {
        walks.extend([&[0, 3][..], &[1, 2]]);
    }
    paths(5, &walks)
}

/// Hamiltonian-cycle numbering of the Frucht graph.
pub fn frucht_graph() -> Graph<()> {
    lcf_graph(12, &[-5, -2, -4, 2, 5, -2, 2, 5, -2, -5, 4, 2], 1).unwrap()
}

pub fn icosahedral_graph() -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(12, GraphKind::Undirected, 30);
    for i in 0..5 {
        let next = (i + 1) % 5;
        for (u, v) in [
            (10, i),
            (11, i + 5),
            (i, next),
            (i + 5, next + 5),
            (i, i + 5),
            (next, i + 5),
        ] {
            add(&mut builder, u, v);
        }
    }
    builder.build()
}

/// The degree-six Cayley graph on Z₄ × Z₄.
pub fn shrikhande_graph() -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(16, GraphKind::Undirected, 48);
    for x in 0..4 {
        for y in 0..4 {
            for (dx, dy) in [(1, 0), (0, 1), (1, 1)] {
                add(&mut builder, 4 * x + y, 4 * ((x + dx) % 4) + (y + dy) % 4);
            }
        }
    }
    builder.build()
}

pub fn hoffman_singleton_graph() -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(50, GraphKind::Undirected, 175);
    for i in 0..5 {
        for j in 0..5 {
            add(&mut builder, 5 * i + j, 5 * i + (j + 1) % 5);
            add(&mut builder, 25 + 5 * i + j, 25 + 5 * i + (j + 2) % 5);
            for k in 0..5 {
                add(&mut builder, 5 * i + j, 25 + 5 * k + (i * k + j) % 5);
            }
        }
    }
    builder.build()
}

pub fn truncated_cube_graph() -> Graph<()> {
    truncate(&cubical_graph())
}
pub fn truncated_tetrahedron_graph() -> Graph<()> {
    truncate(&tetrahedral_graph())
}

fn truncate(base: &Graph<()>) -> Graph<()> {
    let embedding = check_planarity(base).unwrap();
    let mut builder = GraphBuilder::new(2 * base.edge_count(), GraphKind::Undirected);
    for edge in 0..base.edge_count() {
        add(&mut builder, 2 * edge, 2 * edge + 1);
    }
    for u in base.nodes() {
        let around = embedding.clockwise_neighbors(u).unwrap();
        let dart = |id: usize| 2 * id + usize::from(base.edges()[id].from != u);
        for i in 0..around.len() {
            add(
                &mut builder,
                dart(around[i].edge),
                dart(around[(i + 1) % around.len()].edge),
            );
        }
    }
    builder.build()
}

/// Four identical sectors of the Hamiltonian drawing, each spanning three vertices.
pub fn chvatal_graph() -> Graph<()> {
    let mut b = GraphBuilder::new(12, GraphKind::Undirected);
    for v in 0..12 {
        add(&mut b, v, (v + 1) % 12);
        match v % 3 {
            0 => add(&mut b, v, (v + 4) % 12),
            2 => add(&mut b, v, (v + 3) % 12),
            _ => {}
        }
        if v < 6 && v % 3 != 2 {
            add(&mut b, v, v + 6);
        }
    }
    b.build()
}

pub fn krackhardt_kite_graph() -> Graph<()> {
    // Two wings share the central actor. The tail leaves through their tips.
    paths(
        10,
        &[
            &[2, 0, 1, 4, 6, 1, 3, 4],
            &[0, 3, 2, 5, 0],
            &[5, 3, 6, 5, 7, 6],
            &[7, 8, 9],
        ],
    )
}

pub fn sedgewick_maze_graph() -> Graph<()> {
    // The five-room circuit, a shortcut with a dead end, and a triangular room.
    paths(8, &[&[0, 2, 6, 4, 5, 0], &[0, 7, 4], &[7, 1], &[4, 3, 5]])
}

/// Three copies of Tutte's 15-vertex fragment surround one hub.
/// The first row of each fragment contains its two side ports (0 and 2).
pub fn tutte_graph() -> Graph<()> {
    let mut b = GraphBuilder::new(46, GraphKind::Undirected);
    let paths: &[&[usize]] = &[
        &[0, 1, 2, 7, 13, 14, 11, 8, 3, 0], // perimeter, ending at the left port
        &[3, 4, 5, 6, 7],
        &[8, 9, 10, 12, 11],
        &[13, 12],
        &[9, 4],
        &[5, 1],
        &[10, 6],
    ];
    for copy in 0..3 {
        let offset = 15 * copy;
        for path in paths {
            for edge in path.windows(2) {
                add(&mut b, offset + edge[0], offset + edge[1]);
            }
        }
        add(&mut b, 45, offset + 14);
        add(&mut b, offset + 2, 15 * ((copy + 1) % 3));
    }
    b.build()
}

/// The Margulis–Gabber–Galil expander on Zₙ × Zₙ. The graph is 8-regular
/// counting loops twice and retains all parallel edges.
pub fn margulis_gabber_galil_graph(n: usize) -> Result<Graph<()>, GraphError> {
    let vertices = n.checked_mul(n).ok_or(GraphError::GraphSizeOverflow)?;
    let edges = vertices
        .checked_mul(4)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(vertices, GraphKind::Undirected, edges);
    for x in 0..n {
        for y in 0..n {
            for (u, v) in [
                ((x + 2 * y) % n, y),
                ((x + 2 * y + 1) % n, y),
                (x, (y + 2 * x) % n),
                (x, (y + 2 * x + 1) % n),
            ] {
                add(&mut builder, x * n + y, u * n + v);
            }
        }
    }
    Ok(builder.build())
}

/// Dorogovtsev–Goltsev–Mendes deterministic growth model. At each
/// iteration a new vertex is attached to both endpoints of every old edge.
pub fn dorogovtsev_goltsev_mendes_graph(iterations: usize) -> Result<Graph<()>, GraphError> {
    let exponent = u32::try_from(iterations).map_err(|_| GraphError::GraphSizeOverflow)?;
    let edges = 3_usize
        .checked_pow(exponent)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let n = edges / 2 + 2;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, edges);
    let mut endpoints = vec![(0, 1)];
    for (u, v) in endpoints.iter().copied() {
        add(&mut builder, u, v);
    }
    let mut next = 2;
    for _ in 0..iterations {
        let old = endpoints.len();
        for id in 0..old {
            let (u, v) = endpoints[id];
            add(&mut builder, u, next);
            add(&mut builder, v, next);
            endpoints.extend([(u, next), (v, next)]);
            next += 1;
        }
    }
    Ok(builder.build())
}
