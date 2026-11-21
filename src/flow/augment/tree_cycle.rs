//! Tree-cycle descent for the ℓ1 interior-point approach.
//! Fractional descent guides a feasible integral circulation. The core-cycle
//! solver continues from its residual capacities and certifies exact recovery.
use super::{MaxFlow, Network, Residual, residual_neighbors, validate_terminals};
use crate::{
    Graph, GraphError, NONE,
    flow::{CycleRatioWeight, EmbeddedCycleHierarchy},
};

/// Spanning trees whose fundamental cycles are scored at each descent step.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TreeCycleBasis {
    ResidualBfs,
    /// Score both the residual BFS tree and a multiscale low-stretch forest.
    ResidualAndLowStretch,
    /// Maintain a signed core-cycle hierarchy for quantized barrier metrics.
    /// Each returned direction has a global ratio certificate on that graph.
    /// If a certificate exceeds the integer precision, barrier descent ends
    /// and exact optimization continues from the current feasible flow.
    DynamicHierarchy,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TreeCycleFlow {
    pub flow: MaxFlow,
    /// Number of accepted fractional potential-decreasing updates.
    pub descent_steps: usize,
    /// Feasible integral throughput already present when dyadic optimization starts.
    pub value_before_completion: u128,
}

#[derive(Clone, Copy)]
struct Step {
    edge: usize,
    reverse: bool,
}
impl Step {
    fn flipped(self) -> Self {
        Self {
            edge: self.edge,
            reverse: !self.reverse,
        }
    }
    fn sign(self) -> f64 {
        if self.reverse { -1.0 } else { 1.0 }
    }
}

#[derive(Clone, Copy)]
struct Barrier {
    flow: f64,
    capacity: f64,
    length: f64,
    gradient: f64,
}

impl Barrier {
    fn new(capacity: f64) -> Self {
        Self {
            flow: 0.0,
            capacity,
            length: 0.0,
            gradient: 0.0,
        }
    }
    fn slacks(&self, margin: f64) -> (f64, f64) {
        (margin + self.flow, (self.capacity - self.flow) + margin)
    }
    fn refresh(&mut self, margin: f64, alpha: f64) {
        let (lower, upper) = self.slacks(margin);
        let a = lower.powf(-1.0 - alpha);
        let b = upper.powf(-1.0 - alpha);
        self.length = a + b;
        self.gradient = alpha * (b - a);
    }
}

// The tree first covers all integral residual-reachable vertices, then spans
// their entire weak component. In particular the closing arc has a usable
// integral fundamental cycle whenever the sink remains reachable.
struct Basis {
    parent: Vec<usize>,
    incoming: Vec<Step>,
    depth: Vec<usize>,
    ancestors: Vec<Vec<usize>>,
    length: Vec<f64>,
    gradient: Vec<f64>,
}

impl Basis {
    fn empty(n: usize, s: usize) -> Self {
        let mut tree = Self {
            parent: vec![NONE; n],
            incoming: vec![
                Step {
                    edge: NONE,
                    reverse: false
                };
                n
            ],
            depth: vec![0; n],
            ancestors: vec![],
            length: vec![0.; n],
            gradient: vec![0.; n],
        };
        tree.parent[s] = s;
        tree
    }

    fn attach(&mut self, parent: usize, child: usize, step: Step, edges: &[Barrier]) {
        self.parent[child] = parent;
        self.incoming[child] = step;
        self.depth[child] = self.depth[parent] + 1;
        self.length[child] = self.length[parent] + edges[step.edge].length;
        self.gradient[child] = self.gradient[parent] + step.sign() * edges[step.edge].gradient;
    }

    fn build(state: &Network<'_>, edges: &[Barrier], s: usize, t: usize) -> Option<Self> {
        let mut tree = Self::empty(state.graph.node_count(), s);
        let mut queue = vec![s];
        let mut height = 0;
        for residual_only in [true, false] {
            let mut head = 0;
            while let Some(&v) = queue.get(head) {
                head += 1;
                for a in residual_neighbors(state.graph, v) {
                    if tree.parent[a.to] != NONE
                        || state.graph.edges()[a.edge].weight == 0
                        || (residual_only && state.capacity(a) == 0)
                    {
                        continue;
                    }
                    let step = Step {
                        edge: a.edge,
                        reverse: a.reversed != 0,
                    };
                    tree.attach(v, a.to, step, edges);
                    height = height.max(tree.depth[a.to]);
                    queue.push(a.to);
                }
            }
            if tree.parent[t] == NONE {
                return None;
            }
        }
        tree.index_ancestors(height);
        Some(tree)
    }

    fn index_ancestors(&mut self, height: usize) {
        self.ancestors.push(self.parent.clone());
        while (1usize << self.ancestors.len()) <= height {
            let prev = self.ancestors.last().unwrap();
            self.ancestors.push(
                prev.iter()
                    .map(|&p| if p == NONE { NONE } else { prev[p] })
                    .collect(),
            );
        }
    }

    // Quantization affects only the tree topology. Both bases are scored
    // with the same unrounded lengths and gradients.
    fn low_stretch(g: &Graph<u64>, edges: &[Barrier], s: usize) -> Self {
        let smallest = edges[..g.edge_count()]
            .iter()
            .filter(|e| e.capacity > 0.)
            .map(|e| e.length)
            .fold(f64::INFINITY, f64::min);
        let mut input = Vec::with_capacity(g.edge_count());
        for (id, e) in g.edges().iter().enumerate() {
            if e.weight != 0 {
                let length = (edges[id].length / smallest)
                    .ceil()
                    .clamp(1., u64::MAX as f64) as u64;
                input.push((e.from, e.to, length, id));
            }
        }
        let ids = crate::connectivity::low_stretch::build(g.node_count(), input.into_iter());
        let mut b = Graph::builder(g.node_count(), crate::GraphKind::Undirected);
        for id in ids {
            let e = &g.edges()[id];
            b.add_edge(e.from, e.to, id).unwrap();
        }
        let forest = b.build();
        let mut tree = Self::empty(g.node_count(), s);
        let mut stack = vec![s];
        let mut height = 0;
        while let Some(u) = stack.pop() {
            for a in forest.out(u) {
                if tree.parent[a.node] != NONE {
                    continue;
                }
                let edge = forest.edges()[a.edge].weight;
                tree.attach(
                    u,
                    a.node,
                    Step {
                        edge,
                        reverse: g.edges()[edge].from != u,
                    },
                    edges,
                );
                height = height.max(tree.depth[a.node]);
                stack.push(a.node);
            }
        }
        tree.index_ancestors(height);
        tree
    }

    fn lca(&self, mut a: usize, mut b: usize) -> usize {
        if self.depth[a] < self.depth[b] {
            std::mem::swap(&mut a, &mut b);
        }
        let gap = self.depth[a] - self.depth[b];
        for (i, jump) in self.ancestors.iter().enumerate() {
            if (gap >> i) & 1 != 0 {
                a = jump[a];
            }
        }
        if a == b {
            return a;
        }
        for jump in self.ancestors.iter().rev() {
            if jump[a] != jump[b] {
                a = jump[a];
                b = jump[b];
            }
        }
        self.parent[a]
    }

    fn best_cycle(
        &self,
        g: &Graph<u64>,
        edges: &[Barrier],
        s: usize,
        t: usize,
    ) -> (f64, Vec<Step>) {
        let mut best = 0.0;
        let mut chord = None;
        let endpoints = g.edges().iter().map(|e| (e.from, e.to)).chain([(t, s)]);
        for (id, (u, v)) in endpoints.enumerate() {
            if u == v
                || edges[id].capacity == 0.0
                || self.parent[u] == NONE
                || self.incoming[u].edge == id
                || self.incoming[v].edge == id
            {
                continue;
            }
            let base = self.lca(u, v);
            let length = edges[id].length
                + (self.length[u] - self.length[base])
                + (self.length[v] - self.length[base]);
            let gradient = edges[id].gradient + (self.gradient[u] - self.gradient[v]);
            let quality = gradient.abs() / length;
            if quality.is_finite() && quality > best {
                best = quality;
                chord = Some((id, u, v, base, gradient > 0.0));
            }
        }
        let Some((edge, mut u, mut v, base, reverse)) = chord else {
            return (0., vec![]);
        };
        let mut cycle = vec![Step {
            edge,
            reverse: false,
        }];
        while v != base {
            cycle.push(self.incoming[v].flipped());
            v = self.parent[v];
        }
        while u != base {
            cycle.push(self.incoming[u]);
            u = self.parent[u];
        }
        if reverse {
            for step in &mut cycle {
                *step = step.flipped();
            }
        }
        (best, cycle)
    }
}

struct SignedBarrier {
    hierarchy: EmbeddedCycleHierarchy,
    // Leave zero-capacity edges out of the relaxed direction problem.
    original: Vec<usize>,
    local: Vec<usize>,
    changed: Vec<usize>,
}

fn quantize(edge: Barrier, scale: f64, alpha: f64) -> CycleRatioWeight {
    CycleRatioWeight {
        length: (edge.length * scale).round().max(1.) as u64,
        cost: (edge.gradient * (scale / alpha)).round() as i64,
    }
}

impl SignedBarrier {
    fn build(
        graph: &Graph<u64>,
        s: usize,
        t: usize,
        edges: &[Barrier],
        scale: f64,
        alpha: f64,
    ) -> Option<Self> {
        let mut b = Graph::builder(graph.node_count(), crate::GraphKind::Undirected);
        let mut original = Vec::new();
        let mut local = vec![NONE; edges.len()];
        for (id, e) in edges.iter().enumerate().filter(|(_, e)| e.capacity > 0.) {
            let (from, to) = graph.edge(id).map_or((t, s), |arc| (arc.from, arc.to));
            local[id] = original.len();
            original.push(id);
            b.add_edge(from, to, quantize(*e, scale, alpha)).ok()?;
        }
        let hierarchy = EmbeddedCycleHierarchy::from_owned(b.build()).ok()?;
        Some(Self {
            hierarchy,
            original,
            local,
            changed: Vec::new(),
        })
    }
    fn direction(&self) -> Option<Vec<Step>> {
        let factor = 1 + self.original.len().max(1).ilog2();
        let cycle = self.hierarchy.certified_cycle(factor).ok()??;
        Some(
            cycle
                .steps
                .into_iter()
                .map(|a| Step {
                    edge: self.original[a.edge],
                    reverse: a.reversed,
                })
                .collect(),
        )
    }
}

struct Interior {
    signed: Option<SignedBarrier>,
    edges: Vec<Barrier>,
    margin: f64,
    alpha: f64,
    target: f64,
    coefficient: f64,
}
impl Interior {
    fn new(g: &Graph<u64>, bound: u128) -> Self {
        let mut edges: Vec<_> = g
            .edges()
            .iter()
            .map(|e| Barrier::new(e.weight as f64))
            .collect();
        edges.push(Barrier::new(bound as f64));
        let m = edges.len() as f64;
        let largest = edges.iter().map(|e| e.capacity).fold(2., f64::max);
        Self {
            edges,
            signed: None,
            margin: 0.25 / m,
            alpha: 1. / (1000. * (m * largest).ln()),
            target: bound as f64 + 1.,
            coefficient: 20. * m,
        }
    }
    fn descent(
        &mut self,
        state: &Network<'_>,
        s: usize,
        t: usize,
        method: TreeCycleBasis,
    ) -> Option<Vec<Step>> {
        let closing = self.edges.len() - 1;
        let gap = self.target - self.edges[closing].flow;
        if gap <= 0.0 {
            return None;
        }
        let scale = (1u64 << 28) as f64 * self.margin;
        if let Some(signed) = &mut self.signed {
            let mut changes = Vec::with_capacity(signed.changed.len());
            for id in signed.changed.drain(..) {
                let edge = &mut self.edges[id];
                edge.refresh(self.margin, self.alpha);
                if id == closing {
                    edge.gradient -= self.coefficient / gap;
                }
                changes.push((signed.local[id], quantize(*edge, scale, self.alpha)));
            }
            signed.hierarchy.update(&changes).ok()?;
        } else {
            for edge in &mut self.edges {
                edge.refresh(self.margin, self.alpha);
            }
            self.edges[closing].gradient -= self.coefficient / gap;
        }
        let cycle = if method == TreeCycleBasis::DynamicHierarchy {
            if self.signed.is_none() {
                self.signed = Some(SignedBarrier::build(
                    state.graph,
                    s,
                    t,
                    &self.edges,
                    scale,
                    self.alpha,
                )?);
            }
            self.signed.as_ref()?.direction()?
        } else {
            let basis = Basis::build(state, &self.edges, s, t)?;
            let (residual_score, mut cycle) = basis.best_cycle(state.graph, &self.edges, s, t);
            drop(basis);
            if method == TreeCycleBasis::ResidualAndLowStretch {
                let weighted = Basis::low_stretch(state.graph, &self.edges, s);
                let (weighted_score, weighted_cycle) =
                    weighted.best_cycle(state.graph, &self.edges, s, t);
                if weighted_score > residual_score {
                    cycle = weighted_cycle;
                }
            }
            cycle
        };
        if cycle.is_empty() {
            return None;
        }
        let mut amount = f64::INFINITY;
        for step in &cycle {
            let (lower, upper) = self.edges[step.edge].slacks(self.margin);
            amount = amount.min(if step.reverse { lower } else { upper });
            if step.edge == closing && !step.reverse {
                amount = amount.min(gap);
            }
        }
        // Backtrack on actual potential differences. log1p and expm1 keep a
        // small local change visible even on a graph with a large total barrier.
        for _ in 0..64 {
            amount *= 0.5;
            let mut change = 0.;
            let mut interior = true;
            for step in &cycle {
                let e = self.edges[step.edge];
                let (lower, upper) = e.slacks(self.margin);
                let delta = step.sign() * amount;
                let trial = Barrier {
                    flow: e.flow + delta,
                    ..e
                };
                let (lo, hi) = trial.slacks(self.margin);
                interior &= lo > 0. && hi > 0. && trial.flow != e.flow;
                change +=
                    lower.powf(-self.alpha) * (-self.alpha * (delta / lower).ln_1p()).exp_m1();
                change +=
                    upper.powf(-self.alpha) * (-self.alpha * (-delta / upper).ln_1p()).exp_m1();
                if step.edge == closing {
                    change += self.coefficient * (-delta / gap).ln_1p();
                }
            }
            if interior && change.is_finite() && change < 0. {
                for step in &cycle {
                    self.edges[step.edge].flow += step.sign() * amount;
                }
                if let Some(signed) = &mut self.signed {
                    signed.changed.extend(cycle.iter().map(|step| step.edge));
                }
                return Some(cycle);
            }
        }
        None
    }
}

/// Exact hybrid for studying fundamental-cycle interior-point directions.
/// At most `max_steps` descent steps precede exact residual completion.
/// Each static step rebuilds the selected bases and scores their fundamental cycles.
/// The core-cycle solver continues from the resulting feasible residual state.
pub fn max_flow_tree_cycle(
    graph: &Graph<u64>,
    source: usize,
    sink: usize,
    max_steps: usize,
    basis: TreeCycleBasis,
) -> Result<TreeCycleFlow, GraphError> {
    validate_terminals(graph, source, sink)?;
    let mut state = Network::new(graph);
    let out = graph
        .out(source)
        .iter()
        .filter(|a| a.node != source)
        .map(|a| graph.edges()[a.edge].weight as u128)
        .sum::<u128>();
    let inc = graph
        .inc(sink)
        .iter()
        .filter(|a| a.node != sink)
        .map(|a| graph.edges()[a.edge].weight as u128)
        .sum();
    let bound = out.min(inc);
    let mut steps = 0;
    if max_steps != 0 && bound != 0 {
        let mut interior = Interior::new(graph, bound);
        let closing = graph.edge_count();
        while steps < max_steps && state.value < bound {
            let Some(cycle) = interior.descent(&state, source, sink, basis) else {
                break;
            };
            steps += 1;
            // The integral iterate advances along the same signed circulation
            // as far as original residual capacities permit.
            let mut amount = u64::MAX;
            for a in &cycle {
                if a.edge != closing {
                    amount = amount.min(state.capacity[a.edge][usize::from(a.reverse)]);
                } else if a.reverse {
                    amount = amount.min(state.value.min(u128::from(u64::MAX)) as u64);
                }
            }
            for a in cycle {
                if a.edge != closing {
                    let edge = &graph.edges()[a.edge];
                    state.push(
                        Residual {
                            edge: a.edge,
                            to: if a.reverse { edge.from } else { edge.to },
                            reversed: usize::from(a.reverse),
                        },
                        amount,
                    );
                } else if a.reverse {
                    state.value -= u128::from(amount);
                } else {
                    state.value += u128::from(amount);
                }
            }
        }
    }
    let value_before_completion = state.value;
    Ok(TreeCycleFlow {
        flow: crate::flow::core_cycles::finish_from_residual(
            graph,
            source,
            sink,
            state.capacity.into_iter(),
            state.value,
        )?
        .flow,
        descent_steps: steps,
        value_before_completion,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::GraphKind;

    fn potential(interior: &Interior) -> f64 {
        let closing = interior.edges.last().unwrap().flow;
        interior.coefficient * (interior.target - closing).ln()
            + interior
                .edges
                .iter()
                .map(|e| {
                    let (lo, hi) = e.slacks(interior.margin);
                    assert!(lo > 0. && hi > 0.);
                    lo.powf(-interior.alpha) + hi.powf(-interior.alpha)
                })
                .sum::<f64>()
    }

    #[test]
    fn accepted_directions_are_circulations_and_reduce_the_potential() {
        let mut b = Graph::builder(7, GraphKind::Directed);
        for (u, v, c) in [
            (0, 1, 8),
            (0, 2, 5),
            (1, 3, 3),
            (3, 6, 9),
            (2, 3, 9),
            (1, 4, 4),
            (4, 2, 2),
            (5, 4, 3),
            (5, 1, 4),
            (2, 6, 7),
        ] {
            b.add_edge(u, v, c).unwrap();
        }
        let g = b.build();
        let state = Network::new(&g);
        for basis in [
            TreeCycleBasis::ResidualAndLowStretch,
            TreeCycleBasis::DynamicHierarchy,
        ] {
            let mut descent = Interior::new(&g, 13);
            let mut before = potential(&descent);
            let mut accepted = 0;
            for _ in 0..8 {
                let Some(cycle) = descent.descent(&state, 0, 6, basis) else {
                    break;
                };
                let mut balance = vec![0; g.node_count()];
                for step in cycle {
                    let (u, v) = g.edge(step.edge).map_or((6, 0), |e| (e.from, e.to));
                    let sign = if step.reverse { -1 } else { 1 };
                    balance[u] += sign;
                    balance[v] -= sign;
                }
                assert!(balance.iter().all(|&x| x == 0));
                let after = potential(&descent);
                assert!(after < before);
                before = after;
                accepted += 1;
            }
            assert!(accepted > 1);
        }
    }
}
