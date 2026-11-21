use super::{MaxFlow, Network, Residual, residual_at, residual_neighbors, validate_terminals};
use crate::{Graph, GraphError, NONE, structures::LinkCutForest};

fn flush(
    state: &mut Network<'_>,
    forest: &mut LinkCutForest,
    links: &mut [Option<Residual>],
    child: usize,
) {
    let arc = links[child].take().unwrap();
    let remaining = forest.cut(child).unwrap();
    state.push(arc, state.capacity(arc) - remaining);
}

/// Sleator–Tarjan dynamic-tree blocking flow. A path augmentation updates
/// its entire root path lazily. Each admissible arc is linked at most once
/// per BFS phase. O(V(E+V) log(V+2)) time and O(V+E) space.
pub fn max_flow_dynamic_trees(
    graph: &Graph<u64>,
    source: usize,
    sink: usize,
) -> Result<MaxFlow, GraphError> {
    validate_terminals(graph, source, sink)?;
    let mut state = Network::new(graph);
    let n = graph.node_count();
    let mut level = vec![NONE; n];
    let mut queue = Vec::new();
    let mut forest = LinkCutForest::new(n);
    let mut links: Vec<Option<Residual>> = vec![None; n];
    let mut current = vec![0; n];
    loop {
        state.levels(source, &mut level, &mut queue);
        if level[sink] == NONE {
            return Ok(state.finish(source));
        }
        current.fill(0);
        loop {
            let tip = forest.root(source);
            if tip == sink {
                let (capacity, child) = forest.path_minimum(source).unwrap();
                let pushed = forest.subtract_path(source, capacity);
                debug_assert!(pushed);
                state.value += u128::from(capacity);
                flush(&mut state, &mut forest, &mut links, child);
                continue;
            }
            let mut next = None;
            while let Some(arc) = residual_at(graph, tip, current[tip]) {
                current[tip] += 1;
                if state.capacity(arc) != 0 && level[arc.to] == level[tip] + 1 {
                    next = Some(arc);
                    break;
                }
            }
            if let Some(arc) = next {
                let linked = forest.link(tip, arc.to, state.capacity(arc));
                debug_assert!(linked);
                links[tip] = Some(arc);
                continue;
            }
            if tip == source {
                break;
            }
            // A dead end can have incoming tree arcs from several branches.
            // Flush all of them before discarding it from the layered network.
            level[tip] = NONE;
            for arc in residual_neighbors(graph, tip) {
                if links[arc.to].is_some_and(|a| a.to == tip) {
                    flush(&mut state, &mut forest, &mut links, arc.to);
                }
            }
        }
        for u in 0..n {
            if links[u].is_some() {
                flush(&mut state, &mut forest, &mut links, u);
            }
        }
    }
}
