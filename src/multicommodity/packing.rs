//! Fractional path packing with multiplicative capacity prices.
use super::{
    Commodity, ConcurrentFlow, ConcurrentFlowError, FractionalEdgeFlow, validate_commodities,
    verify,
};
use crate::{Graph, NONE, queue::IndexedHeap};
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq)]
pub struct ConcurrentFlowBounds {
    /// The common multiplier and the feasible allocation attaining it.
    pub flow: ConcurrentFlow,
    /// Shortest-path dual upper bound, using floating-point arithmetic.
    pub upper_bound: f64,
    pub rounds: usize,
    pub reached_tolerance: bool,
}

struct Router {
    distance: Vec<f64>,
    parent: Vec<usize>,
    heap: IndexedHeap<u64>,
}
impl Router {
    fn path(
        &mut self,
        g: &Graph<u64>,
        commodity: Commodity,
        price: &[f64],
        path: &mut Vec<usize>,
    ) -> Option<f64> {
        let Commodity { source, sink, .. } = commodity;
        self.distance.fill(f64::INFINITY);
        self.distance[source] = 0.;
        self.heap.clear();
        self.heap.decrease(source, 0);
        path.clear();
        // Nonnegative finite IEEE-754 lengths are ordered by their bits.
        while let Some((bits, u)) = self.heap.pop() {
            let distance = f64::from_bits(bits);
            if u == sink {
                let mut v = u;
                while v != source {
                    let edge = self.parent[v];
                    path.push(edge);
                    v = g.edges()[edge].from;
                }
                return Some(distance);
            }
            for a in g.out(u) {
                if g.edges()[a.edge].weight == 0 {
                    continue;
                }
                let next = distance + price[a.edge];
                if next < self.distance[a.node] {
                    self.distance[a.node] = next;
                    self.parent[a.node] = a.edge;
                    self.heap.decrease(a.node, next.to_bits());
                }
            }
        }
        None
    }
}

/// Multiplicative-weights concurrent flow. Each round routes the demand
/// vector on current shortest paths, scaled so the busiest edge receives
/// one unit of congestion. The cumulative allocation is rescaled to fit.
/// Stops when the feasible multiplier reaches (1-relative_gap) times the
/// best dual upper bound, or returns both bounds after max_rounds.
/// This method requires only the standard library.
pub fn approximate_concurrent_flow(
    graph: &Graph<u64>,
    commodities: &[Commodity],
    relative_gap: f64,
    max_rounds: usize,
) -> Result<ConcurrentFlowBounds, ConcurrentFlowError> {
    validate_commodities(graph, commodities)?;
    if !(relative_gap > 0. && relative_gap < 1.) || max_rounds == 0 {
        return Err(ConcurrentFlowError::InvalidApproximationParameters);
    }
    let m = graph.edge_count();
    let n = graph.node_count();
    let mut mass: Vec<_> = graph
        .edges()
        .iter()
        .map(|e| if e.weight > 0 { 1. } else { 0. })
        .collect();
    let mut price = vec![0.; m];
    let mut routed = vec![0.; m];
    let mut congestion = vec![0.; m];
    let mut allocation = vec![BTreeMap::<usize, f64>::new(); commodities.len()];
    let mut paths = vec![Vec::new(); commodities.len()];
    let mut router = Router {
        distance: vec![0.; n],
        parent: vec![NONE; n],
        heap: IndexedHeap::new(n),
    };
    let mut total = 0.;
    let mut upper = f64::INFINITY;
    for round in 1..=max_rounds {
        routed.fill(0.);
        let mut dual_cost = 0.;
        for (i, e) in graph.edges().iter().enumerate() {
            price[i] = if e.weight == 0 {
                0.
            } else {
                mass[i] / e.weight as f64
            };
            dual_cost += price[i] * e.weight as f64;
        }
        let mut dual_demand = 0.;
        for (path, &commodity) in paths.iter_mut().zip(commodities) {
            let Some(distance) = router.path(graph, commodity, &price, path) else {
                return Ok(ConcurrentFlowBounds {
                    flow: ConcurrentFlow {
                        fraction: 0.,
                        commodity_flows: vec![vec![]; commodities.len()],
                    },
                    upper_bound: 0.,
                    rounds: round,
                    reached_tolerance: true,
                });
            };
            dual_demand += distance * commodity.demand as f64;
            for &edge in path.iter() {
                routed[edge] += commodity.demand as f64;
            }
        }
        // For any lengths, every feasible flow pays at least its multiplier
        // times dual_demand, and pays at most the capacity-weighted dual_cost.
        if dual_demand > 0. {
            upper = upper.min(dual_cost / dual_demand);
        }
        let peak = routed
            .iter()
            .zip(graph.edges())
            .filter(|(_, e)| e.weight != 0)
            .map(|(x, e)| x / e.weight as f64)
            .fold(0., f64::max);
        let step = 1. / peak;
        total += step;
        let mut maximum = 0.0f64;
        let mut largest_mass = 0.0f64;
        for (i, e) in graph.edges().iter().enumerate() {
            if e.weight == 0 {
                continue;
            }
            let increment = step * routed[i] / e.weight as f64;
            mass[i] *= 1. + relative_gap * increment;
            largest_mass = largest_mass.max(mass[i]);
            congestion[i] += increment;
            maximum = maximum.max(congestion[i]);
        }
        for ((path, c), used) in paths.iter().zip(commodities).zip(&mut allocation) {
            let delta = step * c.demand as f64;
            for &e in path {
                *used.entry(e).or_default() += delta;
            }
        }
        let fraction = total / maximum;
        if !fraction.is_finite() || !upper.is_finite() || fraction > upper * (1. + 1e-7) {
            return Err(ConcurrentFlowError::NumericalAccuracy);
        }
        let reached = fraction >= (1. - relative_gap) * upper;
        if reached || round == max_rounds {
            let flow = ConcurrentFlow {
                fraction,
                commodity_flows: allocation
                    .into_iter()
                    .map(|map| {
                        map.into_iter()
                            .map(|(edge, x)| FractionalEdgeFlow {
                                edge,
                                amount: x / maximum,
                            })
                            .collect()
                    })
                    .collect(),
            };
            verify(graph, commodities, &flow)?;
            return Ok(ConcurrentFlowBounds {
                flow,
                upper_bound: upper.max(fraction),
                rounds: round,
                reached_tolerance: reached,
            });
        }
        for w in &mut mass {
            *w /= largest_mass;
        }
    }
    unreachable!("round limit is positive")
}
