use std::collections::HashMap;

use super::add;
use crate::random::Random;
use crate::{Graph, GraphBuilder, GraphError, GraphKind};

#[derive(Clone, Debug)]
pub struct GeometricGraph {
    pub graph: Graph<()>,
    pub positions: Vec<[f64; 2]>,
}

/// Uniform random points in the unit square, with edges of length at most
/// radius. See `geometric_edges` for the spatial indexing strategy.
pub fn random_geometric_graph(
    n: usize,
    radius: f64,
    seed: u64,
) -> Result<GeometricGraph, GraphError> {
    let mut rng = Random(seed);
    let positions: Vec<_> = (0..n).map(|_| [rng.uniform(), rng.uniform()]).collect();
    let edges = geometric_edges(&positions, radius)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, edges.len());
    for (from, to) in edges {
        add(&mut builder, from, to);
    }
    Ok(GeometricGraph {
        graph: builder.build(),
        positions,
    })
}

/// Euclidean radius query for 2D points. A spatial hash examines only
/// adjacent cells. Repeated coordinates remain distinct vertices.
pub fn geometric_edges(
    positions: &[[f64; 2]],
    radius: f64,
) -> Result<Vec<(usize, usize)>, GraphError> {
    if !radius.is_finite() || radius < 0.0 || positions.iter().flatten().any(|v| !v.is_finite()) {
        return Err(GraphError::InvalidModelParameters);
    }
    let mut edges = Vec::new();
    if radius == 0.0 {
        let mut duplicates: HashMap<_, Vec<_>> = HashMap::new();
        for (to, &[x, y]) in positions.iter().enumerate() {
            let key = (float_bits(x), float_bits(y));
            let row = duplicates.entry(key).or_default();
            edges.extend(row.iter().map(|&from| (from, to)));
            row.push(to);
        }
        return Ok(edges);
    }

    // Stay below the range in which floating-point bin coordinates lose
    // integer precision. A wider two-cell neighborhood covers boundary
    // rounding when the exact quotient lies on a cell boundary.
    let bins: Option<Vec<_>> = positions
        .iter()
        .map(|&[x, y]| {
            let x = (x / radius).floor();
            let y = (y / radius).floor();
            (x.abs() <= (1_u64 << 50) as f64 && y.abs() <= (1_u64 << 50) as f64)
                .then_some((x as i64, y as i64))
        })
        .collect();
    if let Some(bins) = bins {
        let mut cells: HashMap<_, Vec<usize>> = HashMap::new();
        for (to, &(x, y)) in bins.iter().enumerate() {
            for dy in -2..=2 {
                for dx in -2..=2 {
                    let Some(candidates) = cells.get(&(x + dx, y + dy)) else {
                        continue;
                    };
                    for &from in candidates {
                        if within(positions[from], positions[to], radius) {
                            edges.push((from, to));
                        }
                    }
                }
            }
            cells.entry((x, y)).or_default().push(to);
        }
    } else {
        // Offsets too large for exact bin indexing still support exact
        // floating-point distance comparisons.
        for (to, &point) in positions.iter().enumerate() {
            for (from, &candidate) in positions[..to].iter().enumerate() {
                if within(candidate, point, radius) {
                    edges.push((from, to));
                }
            }
        }
    }
    Ok(edges)
}

fn float_bits(value: f64) -> u64 {
    if value == 0.0 { 0 } else { value.to_bits() }
}

fn within(a: [f64; 2], b: [f64; 2], radius: f64) -> bool {
    (a[0] - b[0]).hypot(a[1] - b[1]) <= radius
}
