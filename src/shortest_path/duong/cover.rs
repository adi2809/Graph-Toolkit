use super::StepGraph;
use crate::queue::RadixHeap;

pub(super) struct Child {
    pub nodes: Vec<usize>,
    pub diameter: u128,
}

struct Ball {
    vertices: Vec<usize>,
    core_len: usize,
    core_mass: usize,
}

enum Outcome {
    Light(Ball),
    Heavy,
}

struct Search {
    heap: RadixHeap,
    distance: Vec<u128>,
    touched: Vec<usize>,
    vertices: Vec<usize>,
    mass: usize,
    core_mass: usize,
    core_len: usize,
    shell: usize,
    limit: u128,
    current: Option<(u128, usize, usize)>,
    reverse: bool,
}

impl Search {
    fn new(n: usize, reverse: bool) -> Self {
        Self {
            heap: RadixHeap::new(n),
            distance: vec![u128::MAX; n],
            touched: Vec::new(),
            vertices: Vec::new(),
            mass: 0,
            core_mass: 0,
            core_len: 0,
            shell: 0,
            limit: 0,
            current: None,
            reverse,
        }
    }

    fn start(&mut self, center: usize) {
        self.heap.reset();
        for node in self.touched.drain(..) {
            self.distance[node] = u128::MAX;
        }
        self.vertices.clear();
        self.mass = 0;
        self.core_mass = 0;
        self.core_len = 0;
        self.shell = 0;
        self.limit = 0;
        self.current = None;
        self.distance[center] = 0;
        self.touched.push(center);
        self.heap.decrease(center, u128::MAX, 0);
    }

    // Alternating elementary steps charges both competing searches to the
    // ball that is actually accepted. A single step examines at most one arc.
    fn step(&mut self, graph: &StepGraph, active: &[bool], mass: &[usize], limit: u128) -> bool {
        if let Some((cost, node, cursor)) = self.current {
            let adjacency = if self.reverse {
                graph.inc(node)
            } else {
                graph.out(node)
            };
            if let Some(arc) = adjacency.get(cursor) {
                self.current = Some((cost, node, cursor + 1));
                if active[arc.node] {
                    let candidate = cost + graph.edges()[arc.edge].weight.cost.max(0) as u128;
                    if candidate < self.distance[arc.node] {
                        if self.distance[arc.node] == u128::MAX {
                            self.touched.push(arc.node);
                        }
                        self.heap
                            .decrease(arc.node, self.distance[arc.node], candidate);
                        self.distance[arc.node] = candidate;
                    }
                }
                return true;
            }
            self.current = None;
        }
        if self
            .heap
            .peek_key(&self.distance)
            .is_some_and(|cost| cost <= limit)
        {
            let (cost, node) = self.heap.pop(&self.distance).unwrap();
            self.vertices.push(node);
            self.mass += mass[node];
            self.current = Some((cost, node, 0));
            return true;
        }
        false
    }

    fn grow(
        &mut self,
        graph: &StepGraph,
        active: &[bool],
        mass: &[usize],
        scale: &LocalScale,
    ) -> Option<Outcome> {
        if self.step(graph, active, mass, self.limit) {
            return None;
        }
        if self.shell > 0
            && (self.mass - self.core_mass) as u128 * (24 * scale.log_mass) as u128
                <= self.core_mass as u128 * (scale.total / self.core_mass).ilog2() as u128
        {
            return Some(Outcome::Light(Ball {
                vertices: std::mem::take(&mut self.vertices),
                core_len: self.core_len,
                core_mass: self.core_mass,
            }));
        }
        if self.mass >= scale.total.div_ceil(8) {
            return Some(Outcome::Heavy);
        }
        self.core_len = self.vertices.len();
        self.core_mass = self.mass;
        self.shell += 1;
        debug_assert!(self.shell <= scale.shells);
        self.limit = radius(scale.diameter, self.shell, scale.shells * 32);
        None
    }

    fn heavy_ball(
        &mut self,
        graph: &StepGraph,
        mass: &[usize],
        center: usize,
        scale: &LocalScale,
    ) -> (Vec<bool>, Vec<bool>) {
        self.start(center);
        let active = vec![true; mass.len()];
        let shells = 384 * scale.log_mass;
        let denominator = shells * 32;
        let upper = radius(scale.diameter, 2 * shells, denominator);
        while self.step(graph, &active, mass, upper) {}
        let mut cursor = self
            .vertices
            .partition_point(|&u| self.distance[u] <= scale.diameter / 32);
        let mut best = (usize::MAX, 0, 0);
        let mut shell = shells + 1;
        while shell <= 2 * shells {
            let begin = cursor;
            let bound = radius(scale.diameter, shell, denominator);
            let mut amount = 0;
            while cursor < self.vertices.len() && self.distance[self.vertices[cursor]] <= bound {
                amount += mass[self.vertices[cursor]];
                cursor += 1;
            }
            best = best.min((amount, begin, cursor));
            if amount == 0 {
                break;
            }
            shell += 1;
        }
        let mut core = vec![false; mass.len()];
        let mut pad = vec![false; mass.len()];
        for (i, &u) in self.vertices[..best.2].iter().enumerate() {
            core[u] = i < best.1;
            pad[u] = true;
        }
        (core, pad)
    }
}

fn radius(diameter: u128, numerator: usize, denominator: usize) -> u128 {
    let denominator = denominator as u128;
    let numerator = numerator as u128;
    (diameter / denominator) * numerator + (diameter % denominator) * numerator / denominator
}

struct LocalScale {
    diameter: u128,
    shells: usize,
    log_mass: usize,
    total: usize,
}

/// Duong's Lemma 4.9. A light split removes a quarter of the inherited
/// mass. A heavy split certifies a smaller weak diameter for its center.
pub(super) fn decompose(
    graph: &StepGraph,
    mass: &[usize],
    diameter: u128,
    log_mass: usize,
) -> (Vec<Child>, usize) {
    let n = graph.node_count();
    let scale = LocalScale {
        diameter,
        shells: 384 * log_mass * log_mass.next_power_of_two().ilog2().max(1) as usize,
        log_mass,
        total: mass.iter().sum(),
    };
    let mut active = [vec![true; n], vec![true; n]];
    let mut removed = [0, 0];
    let mut search = [Search::new(n, false), Search::new(n, true)];
    let mut children = [Vec::new(), Vec::new()];
    let mut center = 0;
    loop {
        while !active[0][center] || !active[1][center] {
            center += 1;
        }
        for state in &mut search {
            state.start(center);
        }
        let mut heavy = [false; 2];
        let light = loop {
            let mut accepted = None;
            for direction in 0..2 {
                if heavy[direction] {
                    continue;
                }
                match search[direction].grow(graph, &active[direction], mass, &scale) {
                    Some(Outcome::Heavy) => heavy[direction] = true,
                    Some(Outcome::Light(ball)) => {
                        accepted = Some((direction, ball));
                        break;
                    }
                    None => {}
                }
            }
            if accepted.is_some() || heavy == [true, true] {
                break accepted;
            }
        };

        if let Some((direction, ball)) = light {
            for &u in &ball.vertices[..ball.core_len] {
                active[direction][u] = false;
            }
            children[direction].push(Child {
                nodes: ball.vertices,
                diameter,
            });
            removed[direction] += ball.core_mass;
            if removed[direction] >= scale.total.div_ceil(4) {
                let rest = active[direction]
                    .iter()
                    .enumerate()
                    .filter_map(|(u, &is_active)| is_active.then_some(u))
                    .collect();
                children[direction].push(Child {
                    nodes: rest,
                    diameter,
                });
                let mut result = std::mem::take(&mut children[direction]);
                if direction == 1 {
                    result.reverse();
                }
                return (result, 32 * scale.shells);
            }
        } else {
            let (out_core, out_pad) = search[0].heavy_ball(graph, mass, center, &scale);
            let (in_core, in_pad) = search[1].heavy_ball(graph, mass, center, &scale);
            let mut parts = [Vec::new(), Vec::new(), Vec::new()];
            for u in 0..n {
                if out_pad[u] && !in_core[u] {
                    parts[0].push(u);
                }
                if out_pad[u] && in_pad[u] {
                    parts[1].push(u);
                }
                if !out_core[u] {
                    parts[2].push(u);
                }
            }
            let result = parts
                .into_iter()
                .zip([diameter, diameter / 2, diameter])
                .filter(|(nodes, _)| !nodes.is_empty())
                .map(|(nodes, diameter)| Child { nodes, diameter })
                .collect();
            return (result, 32 * 384 * scale.log_mass);
        }
    }
}
