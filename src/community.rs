//! Louvain, Leiden, label propagation and weighted modularity.

mod leiden;
pub use leiden::{LeidenOptions, LeidenResult, leiden_communities, leiden_communities_by};

use std::collections::HashMap;
use std::error::Error;
use std::fmt;

use crate::random::Random;
use crate::{EdgeId, Graph, GraphBuilder, NONE, NodeId};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LouvainOptions {
    /// Positive resolution parameter. Higher values favor smaller communities.
    pub resolution: f64,
    /// Minimum change in modularity required to continue to another level.
    pub threshold: f64,
    /// Seed for the vertex-order shuffle at each level.
    pub seed: u64,
    /// Optional upper bound on aggregation levels. A supplied limit must be positive.
    pub max_levels: Option<usize>,
}

impl Default for LouvainOptions {
    fn default() -> Self {
        Self {
            resolution: 1.0,
            threshold: 1e-7,
            seed: 0,
            max_levels: None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LouvainResult {
    /// Compact community IDs, in original vertex order.
    pub community_of: Vec<usize>,
    pub modularity: f64,
    /// Modularity after each completed level. A singleton result has one value.
    pub level_modularities: Vec<f64>,
    /// Number of local-moving sweeps, including convergence checks.
    pub sweeps: usize,
}

impl LouvainResult {
    pub fn community_count(&self) -> usize {
        self.community_of
            .iter()
            .copied()
            .max()
            .map_or(0, |last| last + 1)
    }

    pub fn communities(&self) -> Vec<Vec<NodeId>> {
        let mut groups = vec![Vec::new(); self.community_count()];
        for (node, &community) in self.community_of.iter().enumerate() {
            groups[community].push(node);
        }
        groups
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum CommunityError {
    RequiresUndirected,
    InvalidTemperature,
    InvalidResolution,
    InvalidThreshold,
    InvalidLevelLimit,
    InvalidWeight { edge: EdgeId },
    InvalidPartition,
}

impl fmt::Display for CommunityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RequiresUndirected => f.write_str("Leiden requires an undirected graph"),
            Self::InvalidTemperature => {
                f.write_str("refinement temperature must be finite and positive")
            }
            Self::InvalidResolution => f.write_str("resolution must be finite and positive"),
            Self::InvalidThreshold => f.write_str("threshold must be finite and nonnegative"),
            Self::InvalidLevelLimit => f.write_str("the level limit must be positive"),
            Self::InvalidWeight { edge } => {
                write!(f, "edge {edge} has a non-finite or negative weight")
            }
            Self::InvalidPartition => f.write_str("provide one community ID in 0..V per vertex"),
        }
    }
}

impl Error for CommunityError {}

/// Computes a seed-reproducible Louvain partition. Both graph kinds,
/// parallel edges and loops are supported. Weights must be finite and
/// nonnegative. A graph with zero total weight yields singleton
/// communities of modularity zero.
///
/// A local-moving sweep scans O(V + E) incidence data. The number of
/// sweeps and aggregation levels depends on the graph and resolution.
pub fn louvain_communities(
    graph: &Graph<f64>,
    options: LouvainOptions,
) -> Result<LouvainResult, CommunityError> {
    louvain_communities_by(graph, |weight| *weight, options)
}

/// Louvain with an edge-payload accessor. For an unweighted graph, pass
/// `|_| 1.0`. The accessor is evaluated once per logical edge.
pub fn louvain_communities_by<W>(
    graph: &Graph<W>,
    weight: impl FnMut(&W) -> f64,
    options: LouvainOptions,
) -> Result<LouvainResult, CommunityError> {
    check_resolution(options.resolution)?;
    if !options.threshold.is_finite() || options.threshold < 0.0 {
        return Err(CommunityError::InvalidThreshold);
    }
    if options.max_levels == Some(0) {
        return Err(CommunityError::InvalidLevelLimit);
    }
    let n = graph.node_count();
    let mut result = LouvainResult {
        community_of: (0..n).collect(),
        modularity: 0.0,
        level_modularities: Vec::new(),
        sweeps: 0,
    };
    let Some(mut level) = normalise(graph, weight)? else {
        result.level_modularities.push(0.0);
        return Ok(result);
    };
    let resolution = options.resolution;
    result.modularity = level.score(&result.community_of, resolution);
    let mut rng = Random(options.seed);
    loop {
        let (mut assignment, sweeps) = level.move_vertices(resolution, &mut rng);
        result.sweeps += sweeps;
        let count = compact(&mut assignment);
        if count == level.graph.node_count() {
            if result.level_modularities.is_empty() {
                result.level_modularities.push(result.modularity);
            }
            return Ok(result);
        }
        let modularity = level.score(&assignment, resolution);
        let gain = modularity - result.modularity;
        // Preserve the previous partition if numerical error reverses a gain.
        if gain < 0.0 {
            if result.level_modularities.is_empty() {
                result.level_modularities.push(result.modularity);
            }
            return Ok(result);
        }
        for community in &mut result.community_of {
            *community = assignment[*community];
        }
        result.modularity = modularity;
        result.level_modularities.push(modularity);
        if gain <= options.threshold
            || count == 1
            || options.max_levels == Some(result.level_modularities.len())
        {
            return Ok(result);
        }
        level = level.contract(&assignment, count);
    }
}

/// Weighted modularity of a supplied partition. IDs may have gaps and
/// must be less than the vertex count. Zero total weight gives zero.
pub fn modularity(
    graph: &Graph<f64>,
    community_of: &[usize],
    resolution: f64,
) -> Result<f64, CommunityError> {
    modularity_by(graph, community_of, |weight| *weight, resolution)
}

pub fn modularity_by<W>(
    graph: &Graph<W>,
    community_of: &[usize],
    weight: impl FnMut(&W) -> f64,
    resolution: f64,
) -> Result<f64, CommunityError> {
    check_resolution(resolution)?;
    let n = graph.node_count();
    if community_of.len() != n || community_of.iter().any(|&id| id >= n) {
        return Err(CommunityError::InvalidPartition);
    }
    Ok(normalise(graph, weight)?.map_or(0.0, |level| level.score(community_of, resolution)))
}

fn check_resolution(resolution: f64) -> Result<(), CommunityError> {
    if !resolution.is_finite() || resolution <= 0.0 {
        Err(CommunityError::InvalidResolution)
    } else {
        Ok(())
    }
}

/// Normalise total logical edge weight to 1 in two passes. Dividing by the
/// largest weight first keeps the intermediate sum within floating range.
fn normalise<W>(
    graph: &Graph<W>,
    mut weight: impl FnMut(&W) -> f64,
) -> Result<Option<Level>, CommunityError> {
    let mut weights = Vec::with_capacity(graph.edge_count());
    let mut maximum = 0.0_f64;
    for (id, edge) in graph.edges().iter().enumerate() {
        let value = weight(&edge.weight);
        if !value.is_finite() || value < 0.0 {
            return Err(CommunityError::InvalidWeight { edge: id });
        }
        maximum = maximum.max(value);
        weights.push(value);
    }
    if maximum == 0.0 {
        return Ok(None);
    }
    let mut sum = 0.0;
    for amount in &mut weights {
        *amount /= maximum;
        sum += *amount;
    }
    let mut builder =
        GraphBuilder::with_edge_capacity(graph.node_count(), graph.kind(), graph.edge_count());
    for (edge, amount) in graph.edges().iter().zip(weights) {
        let normalized = amount / sum;
        if normalized > 0.0 {
            builder.add_edge(edge.from, edge.to, normalized).unwrap();
        }
    }
    Ok(Some(Level::new(builder.build())))
}

fn compact(assignment: &mut [usize]) -> usize {
    let mut index = vec![NONE; assignment.len()];
    let mut count = 0;
    for group in assignment {
        if index[*group] == NONE {
            index[*group] = count;
            count += 1;
        }
        *group = index[*group];
    }
    count
}

struct Level {
    graph: Graph<f64>,
    // Directed outgoing/incoming fractions. For undirected graphs the
    // two entries equal k(v)/(2*m). Loops contribute full undirected mass.
    degrees: Vec<[f64; 2]>,
}

impl Level {
    fn new(graph: Graph<f64>) -> Self {
        let mut degrees = vec![[0.0; 2]; graph.node_count()];
        for edge in graph.edges() {
            let (u, v, w) = (edge.from, edge.to, edge.weight);
            if graph.is_directed() {
                degrees[u][0] += w;
                degrees[v][1] += w;
            } else if u == v {
                degrees[u][0] += w;
            } else {
                degrees[u][0] += 0.5 * w;
                degrees[v][0] += 0.5 * w;
            }
        }
        if !graph.is_directed() {
            for degree in &mut degrees {
                degree[1] = degree[0];
            }
        }
        Self { graph, degrees }
    }

    fn score(&self, assignment: &[usize], resolution: f64) -> f64 {
        let mut total = vec![[0.0; 2]; assignment.len()];
        for (degree, &group) in self.degrees.iter().zip(assignment) {
            total[group][0] += degree[0];
            total[group][1] += degree[1];
        }
        let inner: f64 = self
            .graph
            .edges()
            .iter()
            .filter(|edge| assignment[edge.from] == assignment[edge.to])
            .map(|edge| edge.weight)
            .sum();
        let expected: f64 = total.iter().map(|x| x[0] * x[1]).sum();
        inner - resolution * expected.min(1.0)
    }

    fn move_vertices(&self, resolution: f64, rng: &mut Random) -> (Vec<usize>, usize) {
        let n = self.graph.node_count();
        let mut order: Vec<_> = (0..n).collect();
        rng.shuffle(&mut order);
        let mut assignment: Vec<_> = (0..n).collect();
        let mut totals = self.degrees.clone();
        let mut size = vec![1; n];
        let mut links = vec![0.0; n];
        let mut touched = Vec::new();
        let mut sweeps = 0;
        loop {
            let mut moved = false;
            sweeps += 1;
            for &node in &order {
                let original = assignment[node];
                let outgoing = self.graph.out(node);
                let incoming = if self.graph.is_directed() {
                    self.graph.inc(node)
                } else {
                    &[]
                };
                for arc in outgoing
                    .iter()
                    .chain(incoming)
                    .filter(|arc| arc.node != node)
                {
                    let community = assignment[arc.node];
                    if links[community] == 0.0 {
                        touched.push(community);
                    }
                    links[community] += self.graph.edges()[arc.edge].weight;
                }
                if touched.is_empty() {
                    continue;
                }

                let [out, incoming] = self.degrees[node];
                size[original] -= 1;
                if size[original] == 0 {
                    totals[original] = [0.0; 2];
                } else {
                    totals[original][0] = (totals[original][0] - out).max(0.0);
                    totals[original][1] = (totals[original][1] - incoming).max(0.0);
                }
                let score = |total: [f64; 2], adjacent: f64| {
                    let expected = resolution * (out * total[1] + incoming * total[0]);
                    (adjacent - expected, adjacent.max(expected))
                };
                let mut best = original;
                let (mut best_score, mut magnitude) = score(totals[original], links[original]);
                for community in touched.drain(..) {
                    let amount = std::mem::replace(&mut links[community], 0.0);
                    let (value, candidate_magnitude) = score(totals[community], amount);
                    let epsilon = 32.0 * f64::EPSILON * magnitude.max(candidate_magnitude);
                    if value > best_score + epsilon {
                        best = community;
                        best_score = value;
                        magnitude = candidate_magnitude;
                    }
                }
                size[best] += 1;
                totals[best][0] += out;
                totals[best][1] += incoming;
                assignment[node] = best;
                moved |= best != original;
            }
            if !moved {
                return (assignment, sweeps);
            }

            // Rebase degree sums between sweeps so accumulated rounding from
            // many removals and insertions cannot steer later moves.
            totals.fill([0.0; 2]);
            for (degree, &group) in self.degrees.iter().zip(&assignment) {
                totals[group][0] += degree[0];
                totals[group][1] += degree[1];
            }
        }
    }

    fn contract(&self, assignment: &[usize], count: usize) -> Self {
        let mut index = HashMap::new();
        let mut edges: Vec<(usize, usize, f64)> = Vec::new();
        for edge in self.graph.edges() {
            let (mut from, mut to) = (assignment[edge.from], assignment[edge.to]);
            if !self.graph.is_directed() && from > to {
                std::mem::swap(&mut from, &mut to);
            }
            let id = *index.entry((from, to)).or_insert(edges.len());
            if let Some(row) = edges.get_mut(id) {
                row.2 += edge.weight;
            } else {
                edges.push((from, to, edge.weight));
            }
        }
        let mut builder = GraphBuilder::with_edge_capacity(count, self.graph.kind(), edges.len());
        for (u, v, w) in edges {
            builder.add_edge(u, v, w).unwrap();
        }
        Self::new(builder.build())
    }
}

/// Asynchronous unweighted label propagation on an undirected graph. Each
/// parallel edge casts one vote; loops are ignored. Keeping the current label
/// on ties makes each change strictly increase the number of agreeing edges.
pub fn label_propagation_communities<W>(
    graph: &Graph<W>,
    seed: u64,
) -> Result<Vec<usize>, crate::GraphError> {
    graph.require(crate::GraphKind::Undirected)?;
    let n = graph.node_count();
    let mut order: Vec<_> = graph.nodes().collect();
    let mut labels = order.clone();
    let mut votes = vec![0; n];
    let mut used = Vec::new();
    Random(seed).shuffle(&mut order);
    let mut active = std::collections::VecDeque::from(order);
    let mut queued = vec![true; n];
    while let Some(u) = active.pop_front() {
        queued[u] = false;
        for arc in graph.out(u).iter().filter(|a| a.node != u) {
            let label = labels[arc.node];
            if votes[label] == 0 {
                used.push(label);
            }
            votes[label] += 1;
        }
        let mut best = labels[u];
        for &label in &used {
            if votes[label] > votes[best] {
                best = label;
            }
        }
        for label in used.drain(..) {
            votes[label] = 0;
        }
        if best == labels[u] {
            continue;
        }
        labels[u] = best;
        // Only neighbors with a different label can be destabilized.
        for arc in graph.out(u) {
            if labels[arc.node] != best && !queued[arc.node] {
                queued[arc.node] = true;
                active.push_back(arc.node);
            }
        }
    }
    compact(&mut labels);
    Ok(labels)
}
