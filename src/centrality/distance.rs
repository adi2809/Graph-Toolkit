use crate::Graph;

/// Degree centrality, using total (in plus out) degree on directed graphs.
/// Parallel incidences and self-loops may produce a value above one.
pub fn degree_centrality<W>(graph: &Graph<W>) -> Vec<f64> {
    let n = graph.node_count();
    if n <= 1 {
        return vec![1.0; n];
    }
    graph
        .nodes()
        .map(|node| graph.degree(node).unwrap() as f64 / (n - 1) as f64)
        .collect()
}

/// Sum of reciprocal incoming hop distances to each vertex. Unreachable
/// vertices and the vertex itself contribute zero.
pub fn harmonic_centrality<W>(graph: &Graph<W>) -> Vec<f64> {
    distance_centrality(graph, |_, _, harmonic| harmonic)
}

/// Closeness from incoming hop distances. The result is adjusted for
/// the fraction of the graph that can reach each vertex.
pub fn closeness_centrality<W>(graph: &Graph<W>) -> Vec<f64> {
    let n = graph.node_count().saturating_sub(1);
    distance_centrality(graph, |count, sum, _| {
        if count == 0 {
            0.0
        } else {
            (count as f64 / sum as f64) * (count as f64 / n as f64)
        }
    })
}

fn distance_centrality<W>(graph: &Graph<W>, score: impl Fn(usize, u128, f64) -> f64) -> Vec<f64> {
    let n = graph.node_count();
    let mut stamp = vec![usize::MAX; n];
    let mut queue = Vec::new();
    graph
        .nodes()
        .map(|center| {
            queue.clear();
            queue.push(center);
            stamp[center] = center;
            let (mut first, mut distance) = (0, 0);
            let (mut sum, mut harmonic) = (0, 0.0);
            while first < queue.len() {
                let end = queue.len();
                distance += 1;
                for i in first..end {
                    for arc in graph.inc(queue[i]) {
                        if stamp[arc.node] != center {
                            stamp[arc.node] = center;
                            queue.push(arc.node);
                        }
                    }
                }
                let added = queue.len() - end;
                sum += added as u128 * distance as u128;
                harmonic += added as f64 / distance as f64;
                first = end;
            }
            score(queue.len() - 1, sum, harmonic)
        })
        .collect()
}

/// Brandes' exact unweighted betweenness algorithm in O(V(V + E)) time
/// and O(V + E) space. Paths that use different parallel edges are counted
/// separately. Endpoints never contribute to their own betweenness.
pub fn betweenness_centrality<W>(graph: &Graph<W>, normalized: bool) -> Vec<f64> {
    let n = graph.node_count();
    let mut result = vec![0.0; n];
    if n < 3 {
        return result;
    }
    let mut distance = vec![0; n];
    let mut stamp = vec![usize::MAX; n];
    let mut count = vec![PathCount::ZERO; n];
    let mut dependency = vec![0.0; n];
    let mut queue = Vec::new();
    for source in graph.nodes() {
        queue.clear();
        queue.push(source);
        stamp[source] = source;
        distance[source] = 0;
        count[source] = PathCount {
            fraction: 1.0,
            exponent: 0,
        };
        dependency[source] = 0.0;
        let mut next = 0;
        while let Some(&node) = queue.get(next) {
            next += 1;
            let hops = distance[node] + 1;
            let paths = count[node];
            for arc in graph.out(node) {
                let target = arc.node;
                if stamp[target] != source {
                    stamp[target] = source;
                    count[target] = PathCount::ZERO;
                    dependency[target] = 0.0;
                    distance[target] = hops;
                    queue.push(target);
                }
                if distance[target] == hops {
                    count[target].add(paths);
                }
            }
        }

        for &node in queue.iter().rev() {
            if node == source {
                continue;
            }
            let factor = 1.0 + dependency[node];
            result[node] += dependency[node];
            for arc in graph.inc(node) {
                let parent = arc.node;
                if stamp[parent] == source && distance[parent] + 1 == distance[node] {
                    dependency[parent] += count[parent].ratio(count[node]) * factor;
                }
            }
        }
    }
    let scale = if normalized {
        1.0 / ((n - 1) as f64 * (n - 2) as f64)
    } else if graph.is_directed() {
        1.0
    } else {
        0.5
    };
    for value in &mut result {
        *value *= scale;
    }
    result
}

// A normalized floating-point count with a wider exponent. Long chains of
// diamonds can have far more shortest paths than a finite f64 can represent.
#[derive(Clone, Copy)]
struct PathCount {
    fraction: f64,
    exponent: i64,
}

impl PathCount {
    const ZERO: Self = Self {
        fraction: 0.0,
        exponent: 0,
    };

    fn add(&mut self, mut other: Self) {
        if self.fraction == 0.0 {
            *self = other;
            return;
        }
        if self.exponent < other.exponent {
            std::mem::swap(self, &mut other);
        }
        self.fraction += other.fraction * exp2(other.exponent - self.exponent);
        if self.fraction >= 2.0 {
            self.fraction *= 0.5;
            self.exponent += 1;
        }
    }

    fn ratio(self, other: Self) -> f64 {
        (self.fraction / other.fraction) * exp2(self.exponent - other.exponent)
    }
}

fn exp2(exponent: i64) -> f64 {
    if exponent == 0 {
        1.0
    } else {
        (exponent as f64).exp2()
    }
}
