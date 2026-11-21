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
    let mut degrees = crate::structures::FenwickTree::new(n);
    degrees
        .replace(0, attachments)
        .ok_or(GraphError::GraphSizeOverflow)?;
    for vertex in 1..=attachments {
        add(&mut builder, 0, vertex);
        degrees
            .replace(vertex, 1)
            .ok_or(GraphError::GraphSizeOverflow)?;
    }
    let mut selected = Vec::with_capacity(attachments);
    for vertex in attachments + 1..n {
        for _ in 0..attachments {
            let target = degrees.select(rng.below(degrees.total())).unwrap();
            let old = degrees.replace(target, 0).unwrap();
            selected.push((target, old));
            add(&mut builder, vertex, target);
        }
        for (target, old) in selected.drain(..) {
            degrees
                .replace(target, old + 1)
                .ok_or(GraphError::GraphSizeOverflow)?;
        }
        degrees
            .replace(vertex, attachments)
            .ok_or(GraphError::GraphSizeOverflow)?;
    }
    Ok(builder.build())
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
/// parameters use the next lower even degree.
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

/// Clustered configuration model. Each pair gives a vertex's number of
/// independent-edge stubs and triangle stubs. A triangle stub contributes
/// two degree incidences. Loops and parallel edges are retained.
pub fn random_clustered_graph(
    degrees: &[(usize, usize)],
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    let (single, triangle) = degrees
        .iter()
        .try_fold((0_usize, 0_usize), |(s, t), &(ds, dt)| {
            Some((s.checked_add(ds)?, t.checked_add(dt)?))
        })
        .ok_or(GraphError::GraphSizeOverflow)?;
    if single % 2 != 0 || triangle % 3 != 0 {
        return Err(GraphError::InvalidModelParameters);
    }
    let capacity = (single / 2)
        .checked_add(triangle)
        .ok_or(GraphError::GraphSizeOverflow)?;
    std::alloc::Layout::array::<crate::Edge<()>>(capacity)
        .map_err(|_| GraphError::GraphSizeOverflow)?;
    let mut singles = Vec::with_capacity(single);
    let mut triangles = Vec::with_capacity(triangle);
    for (u, &(ds, dt)) in degrees.iter().enumerate() {
        singles.extend(std::iter::repeat_n(u, ds));
        triangles.extend(std::iter::repeat_n(u, dt));
    }
    let mut rng = Random(seed);
    rng.shuffle(&mut singles);
    rng.shuffle(&mut triangles);
    let mut builder =
        GraphBuilder::with_edge_capacity(degrees.len(), GraphKind::Undirected, capacity);
    for pair in singles.chunks_exact(2) {
        add(&mut builder, pair[0], pair[1]);
    }
    for group in triangles.chunks_exact(3) {
        add(&mut builder, group[0], group[1]);
        add(&mut builder, group[1], group[2]);
        add(&mut builder, group[2], group[0]);
    }
    Ok(builder.build())
}

/// Directed growing network with linear preferential attachment. Every new
/// vertex adds one edge to an earlier vertex. Runs in O(V) time.
pub fn gn_graph(n: usize, seed: u64) -> Graph<()> {
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Directed, n.saturating_sub(1));
    if n >= 2 {
        add(&mut builder, 1, 0);
        let mut stubs = vec![0, 1];
        let mut rng = Random(seed);
        for node in 2..n {
            let target = stubs[rng.below(stubs.len())];
            add(&mut builder, node, target);
            stubs.extend([node, target]);
        }
    }
    builder.build()
}

/// Directed growing network with redirection. A sampled attachment can
/// redirect to its parent's vertex with the given probability.
pub fn gnr_graph(n: usize, probability: f64, seed: u64) -> Result<Graph<()>, GraphError> {
    validate_prob(probability)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Directed, n.saturating_sub(1));
    let mut parent = vec![0];
    let mut rng = Random(seed);
    for node in 1..n {
        let mut target = rng.below(node);
        if target != 0 && rng.uniform() < probability {
            target = parent[target];
        }
        parent.push(target);
        add(&mut builder, node, target);
    }
    Ok(builder.build())
}

/// Directed growing network with copying. Each new vertex connects to a
/// sampled older vertex and its out-neighbors. Runs in O(V + E) time.
pub fn gnc_graph(n: usize, seed: u64) -> Graph<()> {
    let mut builder = GraphBuilder::new(n, GraphKind::Directed);
    let mut parent = vec![0];
    let mut rng = Random(seed);
    for node in 1..n {
        let mut target = rng.below(node);
        parent.push(target);
        loop {
            add(&mut builder, node, target);
            if target == 0 {
                break;
            }
            target = parent[target];
        }
    }
    builder.build()
}

/// Independent edges between two contiguous parts, using geometric skips.
/// Left IDs are 0..left and right IDs start at left.
pub fn bipartite_gnp_random_graph(
    left: usize,
    right: usize,
    probability: f64,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    validate_prob(probability)?;
    let n = left
        .checked_add(right)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let pairs = left
        .checked_mul(right)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut b = GraphBuilder::new(n, GraphKind::Undirected);
    sample_pairs(pairs, probability, &mut Random(seed), |i| {
        add(&mut b, i / right, left + i % right)
    });
    Ok(b.build())
}

/// Union of d independent uniform perfect matchings on two parts of size n.
/// Parallel edges count toward degree. Sampling costs O(n + nd).
pub fn random_regular_bipartite_multigraph(
    n: usize,
    d: usize,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    let mut b = GraphBuilder::with_edge_capacity(
        n.checked_mul(2).ok_or(GraphError::GraphSizeOverflow)?,
        GraphKind::Undirected,
        n.checked_mul(d).ok_or(GraphError::GraphSizeOverflow)?,
    );
    if n == 0 {
        return Ok(b.build());
    }
    let mut perm: Vec<_> = (n..2 * n).collect();
    let mut rng = Random(seed);
    for _ in 0..d {
        rng.shuffle(&mut perm);
        for (left, &right) in perm.iter().enumerate() {
            add(&mut b, left, right);
        }
    }
    Ok(b.build())
}

/// Binomial random intersection graph. Each vertex independently takes each
/// attribute with probability p. Sharing any attribute creates one simple edge.
pub fn random_intersection_graph(
    n: usize,
    attributes: usize,
    p: f64,
    seed: u64,
) -> Result<Graph<()>, GraphError> {
    validate_prob(p)?;
    if attributes == 0 || p == 0.0 || n < 2 {
        return Ok(empty_graph(n, GraphKind::Undirected));
    }
    if p == 1.0 {
        return complete_graph(n, GraphKind::Undirected);
    }
    let mut b = GraphBuilder::new(n, GraphKind::Undirected);
    let mut seen = HashSet::new();
    let mut members = Vec::new();
    let mut rng = Random(seed);
    for _ in 0..attributes {
        members.clear();
        sample_pairs(n, p, &mut rng, |v| {
            for &u in &members {
                if seen.insert((u, v)) {
                    add(&mut b, u, v);
                }
            }
            members.push(v);
        });
    }
    Ok(b.build())
}

/// Uniform simple labelled d-regular graph, conditional on returning Some.
/// Rejection of loops and repeated pairs from a uniformly shuffled stub
/// pairing is efficient for small degrees. Dense graphs use complements.
/// None means the trial budget expired. Infeasible parameters return an error.
pub fn random_regular_graph(
    n: usize,
    d: usize,
    seed: u64,
    trials: usize,
) -> Result<Option<Graph<()>>, GraphError> {
    if d > n.saturating_sub(1) || (n % 2 == 1 && d % 2 == 1) {
        return Err(GraphError::InvalidModelParameters);
    }
    if d == 0 {
        return Ok(Some(empty_graph(n, GraphKind::Undirected)));
    }
    if d == n - 1 {
        return complete_graph(n, GraphKind::Undirected).map(Some);
    }
    if trials == 0 {
        return Ok(None);
    }
    let inverse = d > n - 1 - d;
    let sparse = d.min(n - 1 - d);
    let size = n.checked_mul(sparse).ok_or(GraphError::GraphSizeOverflow)?;
    let mut stubs: Vec<_> = (0..size).map(|i| i / sparse).collect();
    let mut pairs = HashSet::with_capacity(size / 2);
    let mut rng = Random(seed);
    for _ in 0..trials {
        rng.shuffle(&mut stubs);
        pairs.clear();
        if !stubs
            .chunks_exact(2)
            .all(|p| p[0] != p[1] && pairs.insert((p[0].min(p[1]), p[0].max(p[1]))))
        {
            continue;
        }
        let mut b = GraphBuilder::new(n, GraphKind::Undirected);
        if inverse {
            for u in 0..n {
                for v in u + 1..n {
                    if !pairs.contains(&(u, v)) {
                        add(&mut b, u, v);
                    }
                }
            }
        } else {
            for p in stubs.chunks_exact(2) {
                add(&mut b, p[0], p[1]);
            }
        }
        return Ok(Some(b.build()));
    }
    Ok(None)
}
