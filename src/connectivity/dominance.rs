use crate::{Graph, GraphError, GraphKind, NONE, NodeId};

/// Immediate dominators of vertices reachable from start. Unreachable
/// vertices have None, and the start vertex is its own dominator.
///
/// Implements the simple Lengauer–Tarjan algorithm in O((V + E) log V)
/// time with iterative DFS and path compression.
pub fn immediate_dominators<W>(
    graph: &Graph<W>,
    start: NodeId,
) -> Result<Vec<Option<NodeId>>, GraphError> {
    graph.require(GraphKind::Directed)?;
    graph.validate_node(start)?;
    let n = graph.node_count();
    let mut rank = vec![NONE; n];
    rank[start] = 0;
    let mut order = vec![start];
    let mut parent = vec![NONE];
    let mut stack = vec![(start, 0)];
    while let Some((node, position)) = stack.last_mut() {
        if let Some(arc) = graph.out(*node).get(*position) {
            *position += 1;
            if rank[arc.node] == NONE {
                rank[arc.node] = order.len();
                order.push(arc.node);
                parent.push(rank[*node]);
                stack.push((arc.node, 0));
            }
        } else {
            stack.pop();
        }
    }

    let size = order.len();
    let mut sets = LinkEval {
        semi: (0..size).collect(),
        label: (0..size).collect(),
        ancestor: vec![NONE; size],
        scratch: Vec::new(),
    };
    // Each vertex needs exactly one bucket entry. Intrusive links avoid
    // allocating small vectors per semidominator.
    let mut buckets = vec![NONE; size];
    let mut next = vec![NONE; size];
    let mut dominator = vec![0; size];
    for node in (1..size).rev() {
        for arc in graph.inc(order[node]) {
            let predecessor = rank[arc.node];
            if predecessor != NONE {
                let candidate = sets.eval(predecessor);
                sets.semi[node] = sets.semi[node].min(sets.semi[candidate]);
            }
        }
        let semi = sets.semi[node];
        next[node] = buckets[semi];
        buckets[semi] = node;

        let link = parent[node];
        sets.ancestor[node] = link;
        let mut member = std::mem::replace(&mut buckets[link], NONE);
        while member != NONE {
            let candidate = sets.eval(member);
            dominator[member] = if sets.semi[candidate] < sets.semi[member] {
                candidate
            } else {
                link
            };
            member = next[member];
        }
    }

    let mut result = vec![None; n];
    result[start] = Some(start);
    for node in 1..size {
        if dominator[node] != sets.semi[node] {
            dominator[node] = dominator[dominator[node]];
        }
        result[order[node]] = Some(order[dominator[node]]);
    }
    Ok(result)
}

struct LinkEval {
    semi: Vec<usize>,
    label: Vec<usize>,
    ancestor: Vec<usize>,
    scratch: Vec<usize>,
}

impl LinkEval {
    fn eval(&mut self, node: usize) -> usize {
        let mut current = node;
        while self.ancestor[current] != NONE && self.ancestor[self.ancestor[current]] != NONE {
            self.scratch.push(current);
            current = self.ancestor[current];
        }
        while let Some(current) = self.scratch.pop() {
            let ancestor = self.ancestor[current];
            if self.semi[self.label[ancestor]] < self.semi[self.label[current]] {
                self.label[current] = self.label[ancestor];
            }
            self.ancestor[current] = self.ancestor[ancestor];
        }
        self.label[node]
    }
}

/// Each reachable graph vertex's dominance frontier, in increasing vertex
/// order. A frontier can include its owning vertex for loop headers.
pub fn dominance_frontiers<W>(
    graph: &Graph<W>,
    start: NodeId,
) -> Result<Vec<Vec<NodeId>>, GraphError> {
    let dominator = immediate_dominators(graph, start)?;
    let mut frontier = vec![Vec::new(); graph.node_count()];
    let mut seen = vec![NONE; graph.node_count()];
    for node in graph.nodes() {
        let Some(immediate) = dominator[node] else {
            continue;
        };
        // For the root, dominance is never strict with respect to itself.
        // The mathematical frontier includes root for reachable incoming edges.
        let stop = if node == start { NONE } else { immediate };
        for arc in graph.inc(node) {
            let mut previous = arc.node;
            while previous != stop && dominator[previous].is_some() && seen[previous] != node {
                seen[previous] = node;
                frontier[previous].push(node);
                if previous == start {
                    break;
                }
                previous = dominator[previous].unwrap();
            }
        }
    }
    Ok(frontier)
}
