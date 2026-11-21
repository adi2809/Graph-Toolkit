//! Undirected weighted Leiden with normalized modularity.
use super::{CommunityError, Level, check_resolution, compact, normalise};
use crate::{Graph, NONE, random::Random};
use std::collections::{BTreeSet, VecDeque};

#[derive(Clone, Copy, Debug)]
pub struct LeidenOptions {
    /// Modularity resolution gamma.
    pub resolution: f64,
    /// Softmax temperature in normalized modularity units.
    pub temperature: f64,
    pub seed: u64,
}
impl Default for LeidenOptions {
    fn default() -> Self {
        Self {
            resolution: 1.,
            temperature: 0.01,
            seed: 0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct LeidenResult {
    pub community_of: Vec<usize>,
    pub modularity: f64,
    /// Number of refinements that reduced the aggregate graph.
    pub contractions: usize,
}

// Accumulate incident mass by group in linear time in the vertex degree.
struct Links {
    mass: Vec<f64>,
    groups: Vec<usize>,
}
impl Links {
    fn new(n: usize) -> Self {
        Self {
            mass: vec![0.; n],
            groups: vec![],
        }
    }
    fn collect(
        &mut self,
        level: &Level,
        node: usize,
        labels: &[usize],
        allow: impl Fn(usize) -> bool,
    ) {
        for g in self.groups.drain(..) {
            self.mass[g] = 0.;
        }
        for a in level.graph.out(node) {
            if a.node == node || !allow(a.node) {
                continue;
            }
            let g = labels[a.node];
            if self.mass[g] == 0. {
                self.groups.push(g);
            }
            self.mass[g] += level.graph.edges()[a.edge].weight;
        }
    }
}

fn visit_order(n: usize, rng: &mut Random) -> Vec<usize> {
    let mut order: Vec<_> = (0..n).collect();
    rng.shuffle(&mut order);
    order
}
fn expected(gamma: f64, a: f64, b: f64) -> f64 {
    (2. * a * b) * gamma
}

fn move_nodes(level: &Level, labels: &mut [usize], gamma: f64, rng: &mut Random) {
    let n = labels.len();
    let mut size = vec![0; n];
    let mut volume = vec![0.; n];
    for (u, &g) in labels.iter().enumerate() {
        size[g] += 1;
        volume[g] += level.degrees[u][0];
    }
    let mut empty: BTreeSet<_> = (0..n).filter(|&g| size[g] == 0).collect();
    let mut queue = VecDeque::from(visit_order(n, rng));
    let mut queued = vec![true; n];
    let mut links = Links::new(n);
    while let Some(u) = queue.pop_front() {
        queued[u] = false;
        let old = labels[u];
        let degree = level.degrees[u][0];
        links.collect(level, u, labels, |_| true);
        size[old] -= 1;
        if size[old] == 0 {
            volume[old] = 0.;
            empty.insert(old);
        } else {
            volume[old] = (volume[old] - degree).max(0.);
        }
        let mut best = old;
        let mut score = links.mass[old] - expected(gamma, degree, volume[old]);
        for &g in links.groups.iter().chain(empty.first()) {
            let gain = links.mass[g] - expected(gamma, degree, volume[g]);
            // Absolute threshold in normalized modularity units guards cancellation.
            if gain > score + 32. * f64::EPSILON {
                best = g;
                score = gain;
            }
        }
        empty.remove(&best);
        size[best] += 1;
        volume[best] += degree;
        labels[u] = best;
        if old != best {
            for a in level.graph.out(u) {
                if labels[a.node] != best && !queued[a.node] {
                    queue.push_back(a.node);
                    queued[a.node] = true;
                }
            }
        }
    }
}

fn refine(level: &Level, coarse: &[usize], options: LeidenOptions, rng: &mut Random) -> Vec<usize> {
    let n = coarse.len();
    let mut labels: Vec<_> = (0..n).collect();
    let mut size = vec![1; n];
    let mut volume: Vec<_> = level.degrees.iter().map(|d| d[0]).collect();
    let mut parent_volume = vec![0.; n];
    let mut boundary = vec![0.; n];
    for (u, &c) in coarse.iter().enumerate() {
        parent_volume[c] += volume[u];
        for a in level.graph.out(u) {
            if a.node != u && coarse[a.node] == c {
                boundary[u] += level.graph.edges()[a.edge].weight;
            }
        }
    }
    let mut links = Links::new(n);
    let mut choices = Vec::new();
    let mut any_candidate = None;
    let mut merged = false;
    for u in visit_order(n, rng) {
        let c = coarse[u];
        let eligible = |mass: f64, cut: f64| {
            cut >= expected(options.resolution, mass, (parent_volume[c] - mass).max(0.))
        };
        if size[u] != 1 || !eligible(volume[u], boundary[u]) {
            continue;
        }
        links.collect(level, u, &labels, |v| coarse[v] == c);
        choices.clear();
        choices.push((u, 0.)); // Staying is an admissible zero-gain choice.
        let mut peak = 0.0f64;
        for &g in &links.groups {
            let delta = links.mass[g] - expected(options.resolution, volume[u], volume[g]);
            if delta >= 0. && eligible(volume[g], boundary[g]) {
                choices.push((g, delta));
                peak = peak.max(delta);
                any_candidate = Some((u, g));
            }
        }
        let mut sum = 0.;
        for (_, gain) in &mut choices {
            *gain = ((*gain - peak) / options.temperature).exp();
            sum += *gain;
        }
        let mut draw = rng.uniform() * sum;
        let target = choices
            .iter()
            .find(|(_, w)| {
                draw -= w;
                draw < 0.
            })
            .map_or(u, |x| x.0);
        if target == u {
            continue;
        }
        merged = true;
        labels[u] = target;
        size[target] += 1;
        size[u] = 0;
        volume[target] += volume[u];
        boundary[target] = (boundary[target] + boundary[u] - 2. * links.mass[target]).max(0.);
    }
    // If every draw stayed in place, take one recorded admissible merge.
    // The labels in this case are still singletons.
    if !merged {
        if let Some((u, g)) = any_candidate {
            labels[u] = g;
        }
    }
    labels
}

fn connected_parts(level: &Level, labels: &[usize]) -> Vec<usize> {
    let mut parts = vec![NONE; labels.len()];
    let mut stack = Vec::new();
    for u in 0..labels.len() {
        if parts[u] != NONE {
            continue;
        }
        parts[u] = u;
        stack.push(u);
        while let Some(v) = stack.pop() {
            for a in level.graph.out(v) {
                if parts[a.node] == NONE && labels[a.node] == labels[u] {
                    parts[a.node] = u;
                    stack.push(a.node);
                }
            }
        }
    }
    parts
}

/// One multilevel Leiden pass for finite, nonnegative weights on an undirected
/// graph. Refined groups aggregate vertices while retaining their coarse labels.
/// Returns connected communities with respect to positive-weight edges.
pub fn leiden_communities_by<W>(
    graph: &Graph<W>,
    weight: impl FnMut(&W) -> f64,
    options: LeidenOptions,
) -> Result<LeidenResult, CommunityError> {
    if graph.is_directed() {
        return Err(CommunityError::RequiresUndirected);
    }
    check_resolution(options.resolution)?;
    if !options.temperature.is_finite() || options.temperature <= 0. {
        return Err(CommunityError::InvalidTemperature);
    }
    let mut originals: Vec<_> = (0..graph.node_count()).collect();
    let Some(mut level) = normalise(graph, weight)? else {
        return Ok(LeidenResult {
            community_of: originals,
            modularity: 0.,
            contractions: 0,
        });
    };
    let mut labels = originals.clone();
    let mut rng = Random(options.seed);
    let mut contractions = 0;
    loop {
        move_nodes(&level, &mut labels, options.resolution, &mut rng);
        if compact(&mut labels) < labels.len() {
            let mut refined = refine(&level, &labels, options, &mut rng);
            let count = compact(&mut refined);
            if count < labels.len() {
                let mut inherited = vec![0; count];
                for (u, &g) in refined.iter().enumerate() {
                    inherited[g] = labels[u];
                }
                // Parent IDs may exceed the size of the new graph.
                let mut translation = vec![NONE; labels.len()];
                let mut next = 0;
                for g in &mut inherited {
                    if translation[*g] == NONE {
                        translation[*g] = next;
                        next += 1;
                    }
                    *g = translation[*g];
                }
                for id in &mut originals {
                    *id = refined[*id];
                }
                level = level.contract(&refined, count);
                labels = inherited;
                contractions += 1;
                continue;
            }
            labels = connected_parts(&level, &labels);
        }
        let modularity = level.score(&labels, options.resolution);
        for id in &mut originals {
            *id = labels[*id];
        }
        compact(&mut originals);
        return Ok(LeidenResult {
            community_of: originals,
            modularity,
            contractions,
        });
    }
}

pub fn leiden_communities(
    graph: &Graph<f64>,
    options: LeidenOptions,
) -> Result<LeidenResult, CommunityError> {
    leiden_communities_by(graph, |w| *w, options)
}
