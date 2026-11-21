use super::{StrongComponents, WeakComponents, condensation};
use crate::{Graph, GraphError, NONE, NodeId, graph::validate_node, random::Random};

fn same_component(labels: &[usize], source: usize, target: usize) -> Result<bool, GraphError> {
    validate_node(source, labels.len())?;
    validate_node(target, labels.len())?;
    Ok(labels[source] == labels[target])
}
impl WeakComponents {
    /// Constant-time connectivity query with edge direction ignored.
    pub fn connected(&self, source: NodeId, target: NodeId) -> Result<bool, GraphError> {
        same_component(&self.component_of, source, target)
    }
}
impl StrongComponents {
    /// Constant-time mutual reachability query.
    pub fn strongly_connected(&self, source: NodeId, target: NodeId) -> Result<bool, GraphError> {
        same_component(&self.component_of, source, target)
    }
}

#[derive(Clone, Copy)]
struct Label {
    min_reachable: usize,
    first_in_subtree: usize,
    finish: usize,
}

/// Exact static reachability with GRAIL interval filters on the SCC DAG.
/// With a fixed number of dimensions construction and space are linear.
/// DFS-tree ancestry certifies positive queries. Ambiguous queries use pruned
/// search, taking O(C+A) worst-case for a fixed dimension count on a DAG of C
/// components and A arcs. Zero dimensions enables just the topological filter.
pub struct ReachabilityIndex {
    component: Vec<usize>,
    dag: Graph<()>,
    labels: Vec<Vec<Label>>,
}
impl ReachabilityIndex {
    pub fn new<W>(graph: &Graph<W>, dimensions: usize, seed: u64) -> Self {
        let condensed = condensation(graph);
        let dag = condensed.graph;
        let mut rng = Random(seed);
        let mut order: Vec<_> = dag.nodes().collect();
        let mut children: Vec<Vec<_>> = dag
            .nodes()
            .map(|u| dag.out(u).iter().map(|arc| arc.node).collect())
            .collect();
        let mut labels = Vec::new();
        for _ in 0..dimensions {
            rng.shuffle(&mut order);
            for list in &mut children {
                rng.shuffle(list);
            }
            let mut label = vec![
                Label {
                    min_reachable: NONE,
                    first_in_subtree: NONE,
                    finish: NONE
                };
                order.len()
            ];
            let mut clock = 0;
            let mut stack = Vec::new();
            for &root in &order {
                if label[root].finish != NONE {
                    continue;
                }
                stack.push((root, 0));
                label[root].first_in_subtree = clock;
                while let Some((u, i)) = stack.last_mut() {
                    if let Some(&v) = children[*u].get(*i) {
                        *i += 1;
                        if label[v].finish == NONE {
                            label[v].first_in_subtree = clock;
                            stack.push((v, 0));
                        }
                    } else {
                        let (u, _) = stack.pop().unwrap();
                        label[u].finish = clock;
                        label[u].min_reachable = children[u]
                            .iter()
                            .fold(clock, |low, &v| low.min(label[v].min_reachable));
                        clock += 1;
                    }
                }
            }
            labels.push(label);
        }
        Self {
            component: condensed.components.component_of,
            dag,
            labels,
        }
    }

    fn classify(&self, u: usize, v: usize) -> Option<bool> {
        if u == v {
            return Some(true);
        }
        if u > v {
            return Some(false);
        }
        for dim in &self.labels {
            let (a, b) = (dim[u], dim[v]);
            if b.finish > a.finish || b.min_reachable < a.min_reachable {
                return Some(false);
            }
            if b.finish >= a.first_in_subtree {
                return Some(true);
            }
        }
        None
    }

    pub fn reachable(&self, source: NodeId, target: NodeId) -> Result<bool, GraphError> {
        validate_node(source, self.component.len())?;
        validate_node(target, self.component.len())?;
        let (u, v) = (self.component[source], self.component[target]);
        if let Some(answer) = self.classify(u, v) {
            return Ok(answer);
        }
        let mut stack = vec![u];
        let mut visited = std::collections::HashSet::from([u]);
        while let Some(next) = stack.pop() {
            for arc in self.dag.out(next) {
                match self.classify(arc.node, v) {
                    Some(true) => return Ok(true),
                    None if visited.insert(arc.node) => stack.push(arc.node),
                    _ => {}
                }
            }
        }
        Ok(false)
    }
}

/// O(1) reachability queries using a chain decomposition of the SCC DAG.
/// Practical for narrow DAGs. The table stores one component ID per chain
/// per component, and construction is O(V+E+k(C+A)) for k chains.
pub struct ChainReachabilityIndex {
    component: Vec<usize>,
    chain: Vec<usize>,
    width: usize,
    earliest: Vec<usize>,
}

impl ChainReachabilityIndex {
    /// Returns None when the dense label table exceeds table_budget_bytes.
    /// The condensation and path cover need additional linear working memory.
    pub fn with_budget<W>(graph: &Graph<W>, table_budget_bytes: usize) -> Option<Self> {
        let condensed = condensation(graph);
        let dag = condensed.graph;
        let count = dag.node_count();
        let mut chain = vec![NONE; count];
        let mut width = 0;
        for u in dag.nodes() {
            if chain[u] == NONE {
                chain[u] = width;
                width += 1;
            }
            // Extending to the earliest unclaimed successor keeps the chain
            // available to nearby vertices in the topological order.
            if let Some(next) = dag
                .out(u)
                .iter()
                .map(|a| a.node)
                .filter(|&v| chain[v] == NONE)
                .min()
            {
                chain[next] = chain[u];
            }
        }
        let cells = count.checked_mul(width)?;
        if cells > table_budget_bytes / size_of::<usize>() {
            return None;
        }
        let mut earliest = vec![NONE; cells];
        for u in dag.nodes().rev() {
            let (prefix, successors) = earliest.split_at_mut((u + 1) * width);
            let row = &mut prefix[u * width..];
            for arc in dag.out(u) {
                let v = arc.node;
                if row[chain[v]] <= v {
                    continue;
                }
                let offset = (v - u - 1) * width;
                for (entry, &candidate) in row.iter_mut().zip(&successors[offset..offset + width]) {
                    *entry = (*entry).min(candidate);
                }
            }
            row[chain[u]] = u;
        }
        Some(Self {
            component: condensed.components.component_of,
            chain,
            width,
            earliest,
        })
    }
    pub fn chain_count(&self) -> usize {
        self.width
    }
    pub fn reachable(&self, source: usize, target: usize) -> Result<bool, GraphError> {
        validate_node(source, self.component.len())?;
        validate_node(target, self.component.len())?;
        let (a, b) = (self.component[source], self.component[target]);
        Ok(self.earliest[a * self.width + self.chain[b]] <= b)
    }
}
