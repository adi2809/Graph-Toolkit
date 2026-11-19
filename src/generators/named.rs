//! Named graphs and explicit algebraic constructions.
//!
//! Several fixed edge tables use the vertex numbering in NetworkX.
//! See `src/matching/LICENSE.networkx` for its BSD notice.

use std::collections::HashSet;

use super::{add, cubical_graph, cycle_graph, null_graph, tetrahedral_graph};
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

pub fn bull_graph() -> Graph<()> {
    from_pairs(5, &[(0, 1), (0, 2), (1, 2), (1, 3), (2, 4)])
}

pub fn diamond_graph() -> Graph<()> {
    from_pairs(4, &[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)])
}

pub fn house_graph() -> Graph<()> {
    house(false)
}

pub fn house_x_graph() -> Graph<()> {
    house(true)
}

fn house(diagonals: bool) -> Graph<()> {
    let base = [
        (0, 1),
        (0, 2),
        (1, 3),
        (2, 3),
        (2, 4),
        (3, 4),
        (0, 3),
        (1, 2),
    ];
    from_pairs(5, &base[..if diagonals { 8 } else { 6 }])
}

pub fn frucht_graph() -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(12, GraphKind::Undirected, 18);
    for edge in cycle_graph(7, GraphKind::Undirected).edges() {
        add(&mut builder, edge.from, edge.to);
    }
    for (u, v) in [
        (0, 7),
        (1, 7),
        (2, 8),
        (3, 9),
        (4, 9),
        (5, 10),
        (6, 10),
        (7, 11),
        (8, 11),
        (8, 9),
        (10, 11),
    ] {
        add(&mut builder, u, v);
    }
    builder.build()
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

pub fn chvatal_graph() -> Graph<()> {
    from_pairs(
        12,
        &[
            (0, 1),
            (0, 4),
            (0, 6),
            (0, 9),
            (1, 2),
            (1, 5),
            (1, 7),
            (2, 3),
            (2, 6),
            (2, 8),
            (3, 4),
            (3, 7),
            (3, 9),
            (4, 5),
            (4, 8),
            (5, 10),
            (5, 11),
            (6, 10),
            (6, 11),
            (7, 8),
            (7, 11),
            (8, 10),
            (9, 10),
            (9, 11),
        ],
    )
}

pub fn krackhardt_kite_graph() -> Graph<()> {
    from_pairs(
        10,
        &[
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 5),
            (1, 3),
            (1, 4),
            (1, 6),
            (2, 3),
            (2, 5),
            (3, 4),
            (3, 5),
            (3, 6),
            (4, 6),
            (5, 6),
            (5, 7),
            (6, 7),
            (7, 8),
            (8, 9),
        ],
    )
}

pub fn sedgewick_maze_graph() -> Graph<()> {
    from_pairs(
        8,
        &[
            (0, 2),
            (0, 5),
            (0, 7),
            (1, 7),
            (2, 6),
            (3, 4),
            (3, 5),
            (4, 5),
            (4, 6),
            (4, 7),
        ],
    )
}

/// Tutte's cubic non-Hamiltonian planar graph.
pub fn tutte_graph() -> Graph<()> {
    from_pairs(
        46,
        &[
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 4),
            (1, 26),
            (2, 10),
            (2, 11),
            (3, 18),
            (3, 19),
            (4, 5),
            (4, 33),
            (5, 6),
            (5, 29),
            (6, 7),
            (6, 27),
            (7, 8),
            (7, 14),
            (8, 9),
            (8, 38),
            (9, 10),
            (9, 37),
            (10, 39),
            (11, 12),
            (11, 39),
            (12, 13),
            (12, 35),
            (13, 14),
            (13, 15),
            (14, 34),
            (15, 16),
            (15, 22),
            (16, 17),
            (16, 44),
            (17, 18),
            (17, 43),
            (18, 45),
            (19, 20),
            (19, 45),
            (20, 21),
            (20, 41),
            (21, 22),
            (21, 23),
            (22, 40),
            (23, 24),
            (23, 27),
            (24, 25),
            (24, 32),
            (25, 26),
            (25, 31),
            (26, 33),
            (27, 28),
            (28, 29),
            (28, 32),
            (29, 30),
            (30, 31),
            (30, 33),
            (31, 32),
            (34, 35),
            (34, 38),
            (35, 36),
            (36, 37),
            (36, 39),
            (37, 38),
            (40, 41),
            (40, 44),
            (41, 42),
            (42, 43),
            (42, 45),
            (43, 44),
        ],
    )
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

fn from_pairs(n: usize, pairs: &[(usize, usize)]) -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, pairs.len());
    for &(u, v) in pairs {
        add(&mut builder, u, v);
    }
    builder.build()
}
