use std::collections::{BinaryHeap, HashSet};

use super::{add, advance, candidate_count, complete_graph, empty_graph};
use crate::{Graph, GraphBuilder, GraphError, GraphKind};

use crate::random::Random;

/// Independent edge sampling using geometric skips in expected O(V + E) time.
pub fn fast_gnp_random_graph(
    n: usize,
    probability: f64,
    kind: GraphKind,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    validate_prob(probability)?;
    if probability == 0.0 || n < 2 {
        return Ok(empty_graph(n, kind));
    }
    if probability == 1.0 {
        return complete_graph(n, kind);
    }
    let count = candidate_count(n, kind)?;
    let mut builder = GraphBuilder::new(n, kind);
    let mut rng = Random(seed);
    let mut row = PairRow::new(n, kind);
    sample_pairs(count, probability, &mut rng, |index| {
        let (a, b) = row.pair(index);
        add(&mut builder, a, b);
    });
    Ok(builder.build())
}

pub use fast_gnp_random_graph as gnp_random_graph;
pub use fast_gnp_random_graph as erdos_renyi_graph;

fn sample_pairs(count: usize, probability: f64, rng: &mut Random, mut emit: impl FnMut(usize)) {
    if probability == 0.0 {
        return;
    }
    if probability == 1.0 {
        for index in 0..count {
            emit(index);
        }
        return;
    }
    let log_absence = (-probability).ln_1p();
    let mut index = 0;
    while index < count {
        let skip = ((-rng.uniform()).ln_1p() / log_absence) as usize;
        if skip >= count - index {
            break;
        }
        index += skip;
        emit(index);
        index += 1;
    }
}

struct PairRow {
    n: usize,
    from: usize,
    start: usize,
    kind: GraphKind,
}

impl PairRow {
    fn new(n: usize, kind: GraphKind) -> Self {
        Self {
            n,
            from: 0,
            start: 0,
            kind,
        }
    }

    fn pair(&mut self, index: usize) -> (usize, usize) {
        if self.kind == GraphKind::Directed {
            let from = index / (self.n - 1);
            let rank = index % (self.n - 1);
            return (from, rank + usize::from(rank >= from));
        }
        let mut width = self.n - self.from - 1;
        while index - self.start >= width {
            self.start += width;
            self.from += 1;
            width -= 1;
        }
        (self.from, self.from + 1 + index - self.start)
    }
}

/// Uniform simple graph with exactly m edges. Samples the smaller of the
/// selected set and its complement using Floyd's sampling algorithm.
pub fn gnm_random_graph(
    n: usize,
    m: usize,
    kind: GraphKind,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    let count = candidate_count(n, kind)?;
    if m > count {
        return Err(GraphError::InvalidModelParameters);
    }
    let complement = m > count / 2;
    let take = m.min(count - m);
    let mut rng = Random(seed);
    let mut chosen = HashSet::with_capacity(take);
    for limit in count - take..count {
        let item = rng.below(limit + 1);
        if !chosen.insert(item) {
            chosen.insert(limit);
        }
    }
    let mut builder = GraphBuilder::with_edge_capacity(n, kind, m);
    let mut row = PairRow::new(n, kind);
    let mut insert = |index| {
        let (from, to) = row.pair(index);
        add(&mut builder, from, to);
    };
    if complement {
        (0..count)
            .filter(|i| !chosen.contains(i))
            .for_each(&mut insert);
    } else {
        let mut chosen: Vec<_> = chosen.into_iter().collect();
        chosen.sort_unstable();
        chosen.into_iter().for_each(&mut insert);
    }
    Ok(builder.build())
}

pub use gnm_random_graph as dense_gnm_random_graph;

/// Seeded preferential attachment. A Fenwick tree samples weighted targets
/// without replacement in O(m log V) time, where m is the final edge count.
pub fn barabasi_albert_graph(
    n: usize,
    attachments: usize,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    if attachments == 0 || attachments >= n {
        return Err(GraphError::InvalidModelParameters);
    }
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    let mut rng = Random(seed);
    let mut fenwick = Fenwick(vec![0; n + 1]);
    let mut degree = vec![0; n];
    degree[0] = attachments;
    fenwick.add(0, attachments, false);
    for (vertex, value) in degree.iter_mut().enumerate().take(attachments + 1).skip(1) {
        add(&mut builder, 0, vertex);
        *value = 1;
        fenwick.add(vertex, 1, false);
    }
    let mut mass = attachments
        .checked_mul(2)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut selected = Vec::with_capacity(attachments);
    for vertex in attachments + 1..n {
        for _ in 0..attachments {
            let target = fenwick.select(rng.below(mass));
            selected.push(target);
            mass -= degree[target];
            fenwick.add(target, degree[target], true);
            add(&mut builder, vertex, target);
        }
        for target in selected.drain(..) {
            degree[target] += 1;
            mass += degree[target];
            fenwick.add(target, degree[target], false);
        }
        degree[vertex] = attachments;
        mass += attachments;
        fenwick.add(vertex, attachments, false);
    }
    Ok(builder.build())
}

struct Fenwick(Vec<usize>);

impl Fenwick {
    fn add(&mut self, vertex: usize, value: usize, remove: bool) {
        let mut i = vertex + 1;
        while let Some(entry) = self.0.get_mut(i) {
            if remove {
                *entry -= value;
            } else {
                *entry += value;
            }
            i += i & i.wrapping_neg();
        }
    }

    fn select(&self, mut weight: usize) -> usize {
        let mut position = 0;
        let mut bit = 1 << (self.0.len() - 1).ilog2();
        while bit != 0 {
            let next = position | bit;
            if let Some(&mass) = self.0.get(next) {
                if mass <= weight {
                    weight -= mass;
                    position = next;
                }
            }
            bit >>= 1;
        }
        position
    }
}

/// Uniform labeled tree generated from a uniform Prüfer sequence in O(V).
pub fn random_labeled_tree(n: usize, seed: u64) -> Graph<()> {
    if n < 2 {
        return empty_graph(n, GraphKind::Undirected);
    }
    let mut random = Random(seed);
    let mut degree = vec![1; n];
    let sequence: Vec<_> = (0..n - 2)
        .map(|_| {
            let node = random.below(n);
            degree[node] += 1;
            node
        })
        .collect();
    decode_prufer(degree, sequence)
}

/// Builds a labeled tree with the specified degrees in O(V) time.
/// Returns None if the sequence cannot describe a tree.
pub fn degree_sequence_tree(degrees: &[usize]) -> Option<Graph<()>> {
    let n = degrees.len();
    if n == 1 && degrees[0] == 0 {
        return Some(empty_graph(1, GraphKind::Undirected));
    }
    if n < 2
        || degrees.contains(&0)
        || degrees
            .iter()
            .try_fold(0usize, |sum, &d| sum.checked_add(d - 1))?
            != n - 2
    {
        return None;
    }
    let sequence = degrees
        .iter()
        .enumerate()
        .flat_map(|(node, &degree)| std::iter::repeat_n(node, degree - 1));
    Some(decode_prufer(degrees.to_vec(), sequence))
}

fn decode_prufer(mut degree: Vec<usize>, sequence: impl IntoIterator<Item = usize>) -> Graph<()> {
    let n = degree.len();
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, n - 1);
    let mut scan = 0;
    while degree[scan] != 1 {
        scan += 1;
    }
    let mut leaf = scan;
    for node in sequence {
        add(&mut builder, leaf, node);
        degree[leaf] = 0;
        degree[node] -= 1;
        if degree[node] == 1 && node < scan {
            leaf = node;
        } else {
            scan += 1;
            while degree[scan] != 1 {
                scan += 1;
            }
            leaf = scan;
        }
    }
    let other = degree
        .iter()
        .enumerate()
        .find(|&(node, &d)| node != leaf && d == 1)
        .unwrap()
        .0;
    add(&mut builder, leaf, other);
    builder.build()
}

/// Uniform pairing of shuffled stubs. Loops and parallel edges are retained.
pub fn configuration_model(degrees: &[usize], seed: u64) -> Result<Graph<()>, GraphError> {
    let count = degrees.iter().try_fold(0_usize, |sum, &degree| {
        sum.checked_add(degree).ok_or(GraphError::GraphSizeOverflow)
    })?;
    if count % 2 != 0 {
        return Err(GraphError::InvalidModelParameters);
    }
    let mut stubs = Vec::with_capacity(count);
    for (node, &degree) in degrees.iter().enumerate() {
        stubs.extend(std::iter::repeat_n(node, degree));
    }
    let mut rng = Random(seed);
    rng.shuffle(&mut stubs);
    let mut builder =
        GraphBuilder::with_edge_capacity(degrees.len(), GraphKind::Undirected, count / 2);
    for endpoints in stubs.chunks_exact(2) {
        add(&mut builder, endpoints[0], endpoints[1]);
    }
    Ok(builder.build())
}

/// Realize an undirected simple degree sequence or return None.
pub fn havel_hakimi_graph(degrees: &[usize]) -> Option<Graph<()>> {
    let n = degrees.len();
    if degrees.iter().any(|&degree| degree >= n) {
        return None;
    }
    let mut remaining: BinaryHeap<_> = degrees
        .iter()
        .enumerate()
        .filter_map(|(node, &degree)| (degree > 0).then_some((degree, std::cmp::Reverse(node))))
        .collect();
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    let mut updates = Vec::new();
    while let Some((degree, std::cmp::Reverse(node))) = remaining.pop() {
        if degree > remaining.len() {
            return None;
        }
        for _ in 0..degree {
            let (value, neighbor) = remaining.pop().unwrap();
            add(&mut builder, node, neighbor.0);
            if value > 1 {
                updates.push((value - 1, neighbor));
            }
        }
        remaining.extend(updates.drain(..));
    }
    Some(builder.build())
}

/// Stochastic block model with contiguous block vertex IDs and geometric
/// skips within block pairs. Undirected probability matrices must be symmetric.
pub fn stochastic_block_model(
    sizes: &[usize],
    probabilities: &[Vec<f64>],
    kind: GraphKind,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    let blocks = sizes.len();
    if probabilities.len() != blocks || probabilities.iter().any(|row| row.len() != blocks) {
        return Err(GraphError::InvalidModelParameters);
    }
    let mut start = Vec::with_capacity(blocks + 1);
    start.push(0usize);
    for &size in sizes {
        start.push(
            start
                .last()
                .unwrap()
                .checked_add(size)
                .ok_or(GraphError::GraphSizeOverflow)?,
        );
    }
    let mut rng = Random(seed);
    let mut builder = GraphBuilder::new(start[blocks], kind);
    for (a, row) in probabilities.iter().enumerate() {
        for (b, &p) in row.iter().enumerate() {
            validate_prob(p)?;
            if kind == GraphKind::Undirected {
                if p != probabilities[b][a] {
                    return Err(GraphError::InvalidModelParameters);
                }
                if b < a {
                    continue;
                }
            }
            if p == 0.0 || sizes[a] == 0 || sizes[b] == 0 {
                continue;
            }
            if a == b {
                let count = candidate_count(sizes[a], kind)?;
                let mut row = PairRow::new(sizes[a], kind);
                sample_pairs(count, p, &mut rng, |index| {
                    let (u, v) = row.pair(index);
                    add(&mut builder, start[a] + u, start[b] + v);
                });
            } else {
                let count = sizes[a]
                    .checked_mul(sizes[b])
                    .ok_or(GraphError::GraphSizeOverflow)?;
                sample_pairs(count, p, &mut rng, |index| {
                    add(
                        &mut builder,
                        start[a] + index / sizes[b],
                        start[b] + index % sizes[b],
                    );
                });
            }
        }
    }
    Ok(builder.build())
}

/// Rewire each clockwise ring-lattice edge independently. Odd degree
/// parameters use the next lower even degree as in NetworkX.
pub fn watts_strogatz_graph(
    n: usize,
    degree: usize,
    probability: f64,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    small_world(n, degree, probability, seed, true)
}

/// Ring lattice with independently sampled shortcuts, keeping the lattice edges.
pub fn newman_watts_strogatz_graph(
    n: usize,
    degree: usize,
    probability: f64,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    small_world(n, degree, probability, seed, false)
}

/// Samples Watts–Strogatz graphs up to `attempts` times. Returns None
/// if no sampled graph was connected.
pub fn connected_watts_strogatz_graph(
    n: usize,
    degree: usize,
    probability: f64,
    attempts: usize,
    seed: u64,
) -> Result<Option<Graph<()>>, GraphError> {
    validate_prob(probability)?;
    if n == 0 || degree > n {
        return Err(GraphError::InvalidModelParameters);
    }
    if n > 1 && degree < 2 && degree < n {
        return Ok(None);
    }
    let mut rng = Random(seed);
    for attempt in 0..attempts {
        let next = if attempt == 0 { seed } else { rng.next() };
        let graph = watts_strogatz_graph(n, degree, probability, next)?;
        if crate::connectivity::is_connected(&graph)? {
            return Ok(Some(graph));
        }
    }
    Ok(None)
}

fn small_world(
    n: usize,
    degree: usize,
    probability: f64,
    seed: u64,
    rewire: bool,
) -> Result<Graph<()>, GraphError> {
    validate_prob(probability)?;
    if n == 0 || degree > n {
        return Err(GraphError::InvalidModelParameters);
    }
    if degree == n {
        return complete_graph(n, GraphKind::Undirected);
    }
    let mut neighbors = vec![HashSet::new(); n];
    let mut arcs = Vec::new();
    for offset in 1..=degree / 2 {
        for a in 0..n {
            let b = advance(a, offset, n);
            neighbors[a].insert(b);
            neighbors[b].insert(a);
            arcs.push((a, b));
        }
    }
    let mut rng = Random(seed);
    for id in 0..arcs.len() {
        let (a, b) = arcs[id];
        if rng.uniform() >= probability || neighbors[a].len() == n - 1 {
            continue;
        }
        let target = if neighbors[a].len() < n / 2 {
            loop {
                let target = rng.below(n);
                if target != a && !neighbors[a].contains(&target) {
                    break target;
                }
            }
        } else {
            let missing = n - 1 - neighbors[a].len();
            (0..n)
                .filter(|&v| v != a && !neighbors[a].contains(&v))
                .nth(rng.below(missing))
                .unwrap()
        };
        if rewire {
            neighbors[a].remove(&b);
            neighbors[b].remove(&a);
            arcs[id] = (a, target);
        } else {
            arcs.push((a, target));
        }
        neighbors[a].insert(target);
        neighbors[target].insert(a);
    }
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, arcs.len());
    for (from, to) in arcs {
        add(&mut builder, from, to);
    }
    Ok(builder.build())
}

fn validate_prob(probability: f64) -> Result<(), GraphError> {
    if !(0.0..=1.0).contains(&probability) {
        return Err(GraphError::InvalidModelParameters);
    }
    Ok(())
}

/// Chung–Lu expected-degree graph. Samples decreasing edge probabilities by
/// geometric skipping and rejection, preserving the original vertex IDs.
pub fn expected_degree_graph(
    expected: &[f64],
    self_loops: bool,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    if expected.iter().any(|&x| !x.is_finite() || x < 0.0) {
        return Err(GraphError::InvalidModelParameters);
    }
    let maximum = expected.iter().copied().fold(0.0, f64::max);
    if maximum == 0.0 {
        return Ok(empty_graph(expected.len(), GraphKind::Undirected));
    }
    let sum: f64 = expected.iter().map(|x| x / maximum).sum();
    let mut order: Vec<_> = (0..expected.len()).collect();
    order.sort_unstable_by(|&a, &b| expected[b].total_cmp(&expected[a]).then_with(|| a.cmp(&b)));
    let mut rng = Random(seed);
    let mut builder = GraphBuilder::new(expected.len(), GraphKind::Undirected);
    for (i, &source) in order.iter().enumerate() {
        let mut j = i + usize::from(!self_loops);
        let scale = (expected[source] / maximum) / sum;
        while let Some(&target) = order.get(j) {
            let upper = (scale * expected[target]).min(1.0);
            if upper == 0.0 {
                break;
            }
            if upper < 1.0 {
                let skip = ((-rng.uniform()).ln_1p() / (-upper).ln_1p()) as usize;
                if skip >= order.len() - j {
                    break;
                }
                j += skip;
            }
            let target = order[j];
            let actual = (scale * expected[target]).min(1.0);
            if actual == upper || rng.uniform() < actual / upper {
                add(&mut builder, source, target);
            }
            j += 1;
        }
    }
    Ok(builder.build())
}

/// Stochastic block model with equal within-block and between-block
/// probabilities across all respective blocks.
pub fn random_partition_graph(
    sizes: &[usize],
    p_in: f64,
    p_out: f64,
    kind: GraphKind,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    validate_prob(p_in)?;
    validate_prob(p_out)?;
    let matrix: Vec<_> = (0..sizes.len())
        .map(|i| {
            (0..sizes.len())
                .map(|j| if i == j { p_in } else { p_out })
                .collect()
        })
        .collect();
    stochastic_block_model(sizes, &matrix, kind, seed)
}

pub fn planted_partition_graph(
    communities: usize,
    size: usize,
    p_in: f64,
    p_out: f64,
    kind: GraphKind,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    random_partition_graph(&vec![size; communities], p_in, p_out, kind, seed)
}

/// A k-fold lift formed by an independent random matching between the
/// fibers of each logical edge. Copy i of vertex v has ID v*k+i.
/// Each logical edge, including loops and parallel edges, contributes k
/// lifted edges with the original payload and direction.
pub fn random_k_lift<W: Clone>(
    graph: &Graph<W>,
    k: usize,
    seed: u64,
) -> Result<Graph<W>, GraphError> {
    if k == 0 {
        return Err(GraphError::InvalidModelParameters);
    }
    let n = graph
        .node_count()
        .checked_mul(k)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let capacity = graph
        .edge_count()
        .checked_mul(k)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, graph.kind(), capacity);
    if capacity == 0 {
        return Ok(builder.build());
    }
    let mut rng = Random(seed);
    let mut targets: Vec<_> = (0..k).collect();
    for edge in graph.edges() {
        rng.shuffle(&mut targets);
        for (copy, &target) in targets.iter().enumerate() {
            builder
                .add_edge(
                    edge.from * k + copy,
                    edge.to * k + target,
                    edge.weight.clone(),
                )
                .unwrap();
        }
    }
    Ok(builder.build())
}

/// Uniform pairing of outgoing and incoming stubs. The sequences must
/// describe the same number of vertices and have the same degree sum.
pub fn directed_configuration_model(
    in_degrees: &[usize],
    out_degrees: &[usize],
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    if in_degrees.len() != out_degrees.len() {
        return Err(GraphError::InvalidModelParameters);
    }
    let sum = |degrees: &[usize]| {
        degrees.iter().try_fold(0_usize, |sum, &degree| {
            sum.checked_add(degree).ok_or(GraphError::GraphSizeOverflow)
        })
    };
    let count = sum(in_degrees)?;
    if count != sum(out_degrees)? {
        return Err(GraphError::InvalidModelParameters);
    }
    let mut incoming = Vec::with_capacity(count);
    for (node, &degree) in in_degrees.iter().enumerate() {
        incoming.extend(std::iter::repeat_n(node, degree));
    }
    Random(seed).shuffle(&mut incoming);
    let mut builder =
        GraphBuilder::with_edge_capacity(in_degrees.len(), GraphKind::Directed, count);
    let outgoing = out_degrees
        .iter()
        .enumerate()
        .flat_map(|(node, &degree)| std::iter::repeat_n(node, degree));
    for (from, to) in outgoing.zip(incoming) {
        add(&mut builder, from, to);
    }
    Ok(builder.build())
}
