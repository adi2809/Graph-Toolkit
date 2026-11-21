use super::{add, advance, empty_graph};
use crate::{Graph, GraphBuilder, GraphError, GraphKind};

/// Kneser graph on the k-element subsets of 0..n. Vertices are numbered
/// in colexicographic order and edges join disjoint subsets. The generator
/// enumerates subsets of each complement instead of comparing all vertex pairs.
pub fn kneser_graph(n: usize, k: usize) -> Result<Graph<()>, GraphError> {
    if k > n {
        return Err(GraphError::InvalidModelParameters);
    }
    let vertices = choose(n, k)?;
    std::alloc::Layout::array::<usize>(
        vertices
            .checked_add(1)
            .ok_or(GraphError::GraphSizeOverflow)?,
    )
    .map_err(|_| GraphError::GraphSizeOverflow)?;
    if k == 0 || k > n - k {
        return Ok(empty_graph(vertices, GraphKind::Undirected));
    }
    let degree = choose(n - k, k)?;
    let edges = usize::try_from(vertices as u128 * degree as u128 / 2)
        .map_err(|_| GraphError::GraphSizeOverflow)?;
    std::alloc::Layout::array::<crate::Edge<()>>(edges)
        .map_err(|_| GraphError::GraphSizeOverflow)?;

    let width = k + 1;
    let cells = n
        .checked_add(1)
        .and_then(|rows| rows.checked_mul(width))
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut binomial = vec![0_usize; cells];
    for a in 0..=n {
        binomial[a * width] = 1;
        for b in 1..=k.min(a) {
            binomial[a * width + b] =
                binomial[(a - 1) * width + b - 1] + binomial[(a - 1) * width + b];
        }
    }
    let mut builder = GraphBuilder::with_edge_capacity(vertices, GraphKind::Undirected, edges);
    let mut left: Vec<_> = (0..k).collect();
    let mut right = left.clone();
    let mut complement = Vec::with_capacity(n - k);
    for u in 0..vertices {
        complement.clear();
        let mut next = 0;
        for &value in &left {
            complement.extend(next..value);
            next = value + 1;
        }
        complement.extend(next..n);
        for (i, value) in right.iter_mut().enumerate() {
            *value = i;
        }
        loop {
            let v = right
                .iter()
                .enumerate()
                .map(|(i, &x)| binomial[complement[x] * width + i + 1])
                .sum();
            if u < v {
                add(&mut builder, u, v);
            }
            if !next_combination(&mut right, complement.len()) {
                break;
            }
        }
        next_combination(&mut left, n);
    }
    Ok(builder.build())
}

fn choose(n: usize, k: usize) -> Result<usize, GraphError> {
    let mut answer = 1usize;
    for i in 0..k.min(n - k) {
        answer = usize::try_from(answer as u128 * (n - i) as u128 / (i + 1) as u128)
            .map_err(|_| GraphError::GraphSizeOverflow)?;
    }
    Ok(answer)
}

fn next_combination(subset: &mut [usize], limit: usize) -> bool {
    for i in 0..subset.len() {
        if subset[i] + 1 < subset.get(i + 1).copied().unwrap_or(limit) {
            subset[i] += 1;
            for (j, value) in subset[..i].iter_mut().enumerate() {
                *value = j;
            }
            return true;
        }
    }
    false
}

/// Paley graph on a prime field. Directed edges are differences that are
/// nonzero quadratic residues. The undirected version coalesces opposite arcs.
pub fn paley_graph(prime: usize, kind: GraphKind) -> Result<Graph<()>, GraphError> {
    check_prime(prime)?;
    let residues = prime / 2;
    let symmetric = prime == 2 || prime % 4 == 1;
    let arcs = prime as u128 * residues as u128;
    let edges = if kind == GraphKind::Undirected && symmetric {
        arcs / 2
    } else {
        arcs
    };
    let capacity = usize::try_from(edges).map_err(|_| GraphError::GraphSizeOverflow)?;
    std::alloc::Layout::array::<crate::Edge<()>>(capacity)
        .map_err(|_| GraphError::GraphSizeOverflow)?;
    let steps: Vec<_> = (1..=residues).map(|i| mul_mod(i, i, prime)).collect();
    let mut builder = GraphBuilder::with_edge_capacity(prime, kind, capacity);
    for u in 0..prime {
        for &step in &steps {
            let v = advance(u, step, prime);
            if kind == GraphKind::Directed || !symmetric || u < v {
                add(&mut builder, u, v);
            }
        }
    }
    Ok(builder.build())
}

/// The chordal-cycle multigraph on a prime field. Each vertex emits both
/// cycle incidences and an inverse chord. The inverse of zero is defined
/// as zero. Parallel edges and loops are retained.
pub fn chordal_cycle_graph(prime: usize) -> Result<Graph<()>, GraphError> {
    check_prime(prime)?;
    let capacity = prime.checked_mul(3).ok_or(GraphError::GraphSizeOverflow)?;
    std::alloc::Layout::array::<crate::Edge<()>>(capacity)
        .map_err(|_| GraphError::GraphSizeOverflow)?;
    let mut inverse = vec![0; prime];
    inverse[1] = 1;
    for i in 2..prime {
        inverse[i] = prime - mul_mod(prime / i, inverse[prime % i], prime);
    }
    let mut builder = GraphBuilder::with_edge_capacity(prime, GraphKind::Undirected, capacity);
    for (u, &v) in inverse.iter().enumerate() {
        add(&mut builder, u, advance(u, prime - 1, prime));
        add(&mut builder, u, advance(u, 1, prime));
        add(&mut builder, u, v);
    }
    Ok(builder.build())
}

fn mul_mod(a: usize, b: usize, p: usize) -> usize {
    (a as u128 * b as u128 % p as u128) as usize
}

fn check_prime(p: usize) -> Result<(), GraphError> {
    if p == 2 {
        return Ok(());
    }
    if p < 2 || p % 2 == 0 {
        return Err(GraphError::InvalidModelParameters);
    }
    // Sinclair's seven Miller–Rabin bases are deterministic below 2^64.
    let twos = (p - 1).trailing_zeros();
    let odd = (p - 1) >> twos;
    'base: for witness in [2, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022] {
        let mut factor = witness % p;
        if factor == 0 {
            continue;
        }
        let mut exponent = odd;
        let mut value = 1;
        while exponent != 0 {
            if exponent & 1 != 0 {
                value = mul_mod(value, factor, p);
            }
            factor = mul_mod(factor, factor, p);
            exponent >>= 1;
        }
        if value == 1 || value == p - 1 {
            continue;
        }
        for _ in 1..twos {
            value = mul_mod(value, value, p);
            if value == p - 1 {
                continue 'base;
            }
        }
        return Err(GraphError::InvalidModelParameters);
    }
    Ok(())
}

/// Harary graph H(k,n) with the minimum number of edges for k-connectivity.
/// For an odd k and odd n, one vertex has degree k+1.
pub fn hkn_harary_graph(k: usize, n: usize) -> Result<Graph<()>, GraphError> {
    if k == 0 || k >= n {
        return Err(GraphError::InvalidModelParameters);
    }
    if k == 1 {
        return Ok(super::path_graph(n, GraphKind::Undirected));
    }
    let count = usize::try_from((k as u128 * n as u128).div_ceil(2))
        .map_err(|_| GraphError::GraphSizeOverflow)?;
    std::alloc::Layout::array::<crate::Edge<()>>(count)
        .map_err(|_| GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Undirected, count);
    for step in 1..=k / 2 {
        for u in 0..n {
            add(&mut builder, u, advance(u, step, n));
        }
    }
    if k % 2 == 1 {
        for u in 0..n.div_ceil(2) {
            add(&mut builder, u, u + n / 2);
        }
    }
    Ok(builder.build())
}

/// Graph with vertex and logical edge IDs arranged as a directed tournament.
/// Each unordered pair has one uniformly sampled orientation.
pub fn random_tournament(n: usize, seed: u64) -> Result<Graph<()>, GraphError> {
    if n == 0 {
        return Ok(empty_graph(0, GraphKind::Directed));
    }
    let capacity = super::candidate_count(n, GraphKind::Undirected)?;
    std::alloc::Layout::array::<crate::Edge<()>>(capacity)
        .map_err(|_| GraphError::GraphSizeOverflow)?;
    let mut builder = GraphBuilder::with_edge_capacity(n, GraphKind::Directed, capacity);
    let mut rng = crate::random::Random(seed);
    let mut bits = 0;
    let mut left = 0;
    for u in 0..n {
        for v in u + 1..n {
            if left == 0 {
                bits = rng.next();
                left = 64;
            }
            if bits & 1 == 0 {
                add(&mut builder, u, v);
            } else {
                add(&mut builder, v, u);
            }
            bits >>= 1;
            left -= 1;
        }
    }
    Ok(builder.build())
}
