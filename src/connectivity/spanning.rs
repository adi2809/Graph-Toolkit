//! Minimum-total-weight and minimax spanning structures.

use crate::shortest_path::bottleneck_search;
use crate::{EdgeId, Graph, GraphError, GraphKind, NONE, NodeId};

use super::{DisjointSet, kruskal};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SpanningTree<W> {
    pub edges: Vec<EdgeId>,
    pub total_weight: i128,
    pub bottleneck: Option<W>,
}

/// Returns an undirected minimum-total-weight spanning tree using Kruskal.
/// A disconnected nonempty graph yields `None`. An empty tree has weight zero.
pub fn minimum_spanning_tree<W: Copy + Ord + Into<i128>>(
    graph: &Graph<W>,
) -> Result<Option<SpanningTree<W>>, GraphError> {
    let (edges, components) = kruskal(graph)?;
    if components > 1 {
        return Ok(None);
    }
    let bottleneck = edges.last().map(|&id| graph.edges()[id].weight);
    let weights = edges.iter().map(|&id| graph.edges()[id].weight.into());
    let mut positive = weights.clone().filter(|&weight| weight > 0);
    let mut negative = weights.filter(|&weight| weight < 0);
    let mut total_weight = 0_i128;
    loop {
        // Alternate signs when possible so cancellation never causes a
        // false intermediate overflow, even for i128 edge weights.
        let term = if total_weight >= 0 {
            negative.next().or_else(|| positive.next())
        } else {
            positive.next().or_else(|| negative.next())
        };
        let Some(term) = term else { break };
        total_weight = total_weight
            .checked_add(term)
            .ok_or(GraphError::WeightOverflow)?;
    }
    Ok(Some(SpanningTree {
        edges,
        total_weight,
        bottleneck,
    }))
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BottleneckTree<W> {
    pub edges: Vec<EdgeId>,
    /// Largest selected edge weight. A tree without edges has no bottleneck.
    pub bottleneck: Option<W>,
}

#[derive(Clone, Copy)]
struct Candidate {
    from: NodeId,
    to: NodeId,
    original: EdgeId,
}

/// Finds an undirected spanning tree whose heaviest edge is as light as possible.
/// Median selection and contraction take O(V + E α(V)) time with union-find.
/// A disconnected nonempty graph yields `None`.
pub fn minimum_bottleneck_spanning_tree<W: Copy + Ord>(
    graph: &Graph<W>,
) -> Result<Option<BottleneckTree<W>>, GraphError> {
    graph.require(GraphKind::Undirected)?;
    let mut remaining = graph.node_count();
    let mut candidates: Vec<_> = graph
        .edges()
        .iter()
        .enumerate()
        .filter_map(|(original, edge)| {
            (edge.from != edge.to).then_some(Candidate {
                from: edge.from,
                to: edge.to,
                original,
            })
        })
        .collect();
    let mut selected = Vec::new();

    while remaining > 1 {
        if candidates.len() < remaining - 1 {
            return Ok(None);
        }
        let middle = candidates.len() / 2;
        let (_, middle, _) = candidates
            .select_nth_unstable_by_key(middle, |edge| graph.edges()[edge.original].weight);
        let pivot = graph.edges()[middle.original].weight;
        let mut sets = DisjointSet::new(remaining);
        let mut forest = Vec::new();
        for edge in &candidates {
            if graph.edges()[edge.original].weight < pivot
                && sets.union(edge.from, edge.to).is_some()
            {
                forest.push(edge.original);
            }
        }
        if forest.len() == remaining - 1 {
            // The optimum lies strictly below the pivot.
            candidates.retain(|edge| graph.edges()[edge.original].weight < pivot);
            continue;
        }

        for edge in &candidates {
            if graph.edges()[edge.original].weight == pivot
                && sets.union(edge.from, edge.to).is_some()
            {
                forest.push(edge.original);
            }
        }
        let groups = remaining - forest.len();
        selected.extend(forest);
        if groups == 1 {
            break;
        }

        // All edges at or below the pivot can be contracted. Only heavier
        // edges between the new components can affect the answer.
        let mut names = vec![NONE; remaining];
        let mut next_name = 0;
        let compact: Vec<_> = (0..remaining)
            .map(|node| {
                let root = sets.find(node);
                if names[root] == NONE {
                    names[root] = next_name;
                    next_name += 1;
                }
                names[root]
            })
            .collect();
        candidates.retain_mut(|edge| {
            edge.from = compact[edge.from];
            edge.to = compact[edge.to];
            edge.from != edge.to && graph.edges()[edge.original].weight > pivot
        });
        remaining = groups;
    }
    let bottleneck = selected.iter().map(|&id| graph.edges()[id].weight).max();
    Ok(Some(BottleneckTree {
        edges: selected,
        bottleneck,
    }))
}

/// Minimum-bottleneck spanning arborescence rooted at `root`, oriented outward.
/// Returns `None` if the root cannot reach every vertex.
pub fn minimum_bottleneck_arborescence<W: Copy + Ord>(
    graph: &Graph<W>,
    root: NodeId,
) -> Result<Option<BottleneckTree<W>>, GraphError> {
    graph.require(GraphKind::Directed)?;
    graph.validate_node(root)?;
    let (costs, predecessors) = bottleneck_search(graph, root, None);
    if costs
        .iter()
        .enumerate()
        .any(|(node, cost)| node != root && cost.is_none())
    {
        return Ok(None);
    }
    let bottleneck = costs.into_iter().flatten().max();
    let edges = predecessors
        .into_iter()
        .filter_map(|(_, edge)| (edge != NONE).then_some(edge))
        .collect();
    Ok(Some(BottleneckTree { edges, bottleneck }))
}

/// Wilson's loop-erased random walks, with one root per connected component.
/// Samples a uniform spanning tree independently in each component. Parallel
/// edges are distinct choices. Payloads are ignored. Returns original edge IDs.
pub fn random_spanning_forest<W>(graph: &Graph<W>, seed: u64) -> Result<Vec<EdgeId>, GraphError> {
    let components = super::connected_components(graph)?;
    let mut tree = vec![false; graph.node_count()];
    let mut rooted = vec![false; components.count];
    for (v, &c) in components.component_of.iter().enumerate() {
        if !rooted[c] {
            tree[v] = true;
            rooted[c] = true;
        }
    }
    let mut rng = crate::random::Random(seed);
    let mut successor = vec![None; graph.node_count()];
    let mut result = Vec::with_capacity(graph.node_count() - components.count);
    for root in graph.nodes() {
        let mut v = root;
        while !tree[v] {
            let neighbors = graph.out(v);
            let arc = neighbors[rng.below(neighbors.len())];
            successor[v] = Some(arc);
            v = arc.node;
        }
        // Last exits erase every cycle in the walk, including self-loops.
        v = root;
        while !tree[v] {
            tree[v] = true;
            let arc = successor[v].unwrap();
            result.push(arc.edge);
            v = arc.node;
        }
    }
    Ok(result)
}
