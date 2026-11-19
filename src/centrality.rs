//! Vertex importance measures.

use std::error::Error;
use std::fmt;

use crate::Graph;

mod distance;
pub use distance::{
    betweenness_centrality, closeness_centrality, degree_centrality, harmonic_centrality,
};

/// Options for weighted, personalized PageRank.
#[derive(Clone, Copy, Debug)]
pub struct PageRankOptions<'a> {
    /// Probability of following a graph edge. Must be in [0, 1).
    pub damping: f64,
    /// Absolute L1 tolerance between consecutive rank vectors.
    pub tolerance: f64,
    pub max_iterations: usize,
    /// Nonnegative, finite teleport weights in vertex order. Defaults to uniform.
    pub personalization: Option<&'a [f64]>,
}

impl Default for PageRankOptions<'_> {
    fn default() -> Self {
        Self {
            damping: 0.85,
            tolerance: 1e-9,
            max_iterations: 200,
            personalization: None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub enum PageRankError {
    InvalidDamping,
    InvalidTolerance,
    InvalidPersonalization,
    DidNotConverge { iterations: usize, residual: f64 },
}

impl fmt::Display for PageRankError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::InvalidDamping => f.write_str("PageRank damping must be in [0, 1)"),
            Self::InvalidTolerance => f.write_str("PageRank tolerance must be finite and positive"),
            Self::InvalidPersonalization => {
                f.write_str("provide a finite nonnegative teleport weight per vertex with positive total weight")
            }
            Self::DidNotConverge { iterations, residual } => {
                write!(f, "PageRank did not converge after {iterations} iterations (L1 change {residual})")
            }
        }
    }
}

impl Error for PageRankError {}

/// Computes weighted PageRank by sparse power iteration in O(V + E) per step.
/// Ranks are in vertex order. Vertices with zero outgoing weight redistribute
/// their mass according to the normalized personalization vector.
pub fn pagerank(
    graph: &Graph<u64>,
    options: PageRankOptions<'_>,
) -> Result<Vec<f64>, PageRankError> {
    if !(0.0..1.0).contains(&options.damping) {
        return Err(PageRankError::InvalidDamping);
    }
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        return Err(PageRankError::InvalidTolerance);
    }
    let n = graph.node_count();
    let teleport = teleport_weights(n, options.personalization)?;
    if n == 0 || options.damping == 0.0 {
        return Ok(teleport);
    }
    let scale: Vec<_> = (0..n)
        .map(|node| {
            let sum: u128 = graph
                .out(node)
                .iter()
                .map(|arc| u128::from(graph.edges()[arc.edge].weight))
                .sum();
            if sum == 0 {
                0.0
            } else {
                options.damping / sum as f64
            }
        })
        .collect();
    let mut ranks = teleport.clone();
    let mut next = vec![0.0; n];
    let mut residual = f64::INFINITY;
    for _ in 0..options.max_iterations {
        let dangling: f64 = ranks
            .iter()
            .zip(&scale)
            .filter_map(|(&rank, &factor)| (factor == 0.0).then_some(rank))
            .sum();
        let teleport_mass = (1.0 - options.damping) + options.damping * dangling;
        for (rank, &weight) in next.iter_mut().zip(&teleport) {
            *rank = teleport_mass * weight;
        }
        for (node, (&rank, &factor)) in ranks.iter().zip(&scale).enumerate() {
            if factor == 0.0 {
                continue;
            }
            let per_unit = factor * rank;
            for arc in graph.out(node) {
                next[arc.node] += per_unit * graph.edges()[arc.edge].weight as f64;
            }
        }
        residual = ranks.iter().zip(&next).map(|(a, b)| (a - b).abs()).sum();
        if residual <= options.tolerance {
            return Ok(next);
        }
        std::mem::swap(&mut ranks, &mut next);
    }
    Err(PageRankError::DidNotConverge {
        iterations: options.max_iterations,
        residual,
    })
}

fn teleport_weights(n: usize, input: Option<&[f64]>) -> Result<Vec<f64>, PageRankError> {
    let Some(input) = input else {
        return Ok(if n == 0 {
            Vec::new()
        } else {
            vec![1.0 / n as f64; n]
        });
    };
    if input.len() != n
        || input
            .iter()
            .any(|&weight| !weight.is_finite() || weight < 0.0)
    {
        return Err(PageRankError::InvalidPersonalization);
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    let peak = input.iter().copied().fold(0.0, f64::max);
    if peak == 0.0 {
        return Err(PageRankError::InvalidPersonalization);
    }
    // Scale before summing so finite weights close to f64::MAX can be normalized.
    let mut weights: Vec<_> = input.iter().map(|weight| weight / peak).collect();
    let sum: f64 = weights.iter().sum();
    for weight in &mut weights {
        *weight /= sum;
    }
    Ok(weights)
}
