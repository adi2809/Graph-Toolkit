# Graph Library

A Rust library for graph traversal, routing, matching, and network optimization. Build an immutable graph from a mutable edge list, then reuse it across queries. The graph uses compressed sparse rows (CSR) with contiguous adjacency and stable edge IDs. Directed graphs also index incoming arcs. Parallel edges and self-loops are supported.

```rust
use graph_library::{Graph, GraphKind};
use graph_library::shortest_path::bidirectional_dijkstra;
use graph_library::traversal::bfs;

# fn main() -> Result<(), graph_library::GraphError> {
let mut builder = Graph::builder(4, GraphKind::Directed);
builder.add_edge(0, 1, 4_u64)?;
let first = builder.add_edge(0, 2, 1)?;
let second = builder.add_edge(2, 1, 1)?;
builder.add_edge(1, 3, 3)?;
let graph = builder.build();

assert_eq!(bfs(&graph, 0)?, [0, 1, 2, 3]);
let path = bidirectional_dijkstra(&graph, 0, 3)?.unwrap();
assert_eq!(path.cost, 5);
assert_eq!(path.nodes, [0, 2, 1, 3]);
assert_eq!(&path.edges[..2], &[first, second]);
# Ok(())
# }
```

| Module | Algorithms and data structures |
| --- | --- |
| `traversal` | Insertion-ordered BFS, direction-optimizing BFS, and iterative DFS. Multi-source BFS returns a shortest-path forest with source ownership. |
| `grid` | Implicit grids with bit-packed obstacles, four- or eight-way movement, multi-source distance fields, A*, and four-way jump-point search. |
| `shortest_path` | Indexed-radix and Fibonacci-heap Dijkstra, bidirectional Dijkstra, A*, ALT, 0–1 BFS, signed Bellman–Ford, Duong’s local padded solver, Johnson, blocked Floyd–Warshall, indexed DAG routing, minimax paths, and k-shortest simple paths. |
| `matching` | Edmonds and Gabow cardinality matching, exact weighted blossom, bottleneck matching, sparse bipartite assignment, edge covers, and validation. Kuhn, Hopcroft–Karp, and seeded loop-erased random walks return minimum vertex covers. |
| `flow` | Ford–Fulkerson, capacity scaling, push–relabel, Dinic, dynamic-tree blocking flow, cofacial planar flow, dynamic signed cycle hierarchies, terminal core cycles with exact dyadic recovery, Edmonds–Karp, Boykov–Kolmogorov, cost and capacity scaling for transshipment, network simplex, minimum-ratio cycles, cycle canceling, successive shortest paths, Stoer–Wagner, and Gomory–Hu. |
| `multicommodity` | Maximum concurrent fractional flow for directed graphs. Path packing returns a feasible allocation and a dual upper bound. The optional sparse LP optimizes the common demand multiplier. |
| `connectivity` | Interval and chain reachability indexes, union–find, offline dynamic connectivity, components, chordal recognition with a maximum clique and minimum coloring, edge cuts, dominators, Eulerian traversals, cycle bases, core numbers, DAG ordering and reduction, Kruskal and randomized filtered spanning forests, low-stretch spanning forests, uniform random spanning forests, and bottleneck arborescences. |
| `centrality` | Weighted personalized PageRank, degree, harmonic, closeness, and exact or pivot-sampled unweighted Brandes betweenness. |
| `clustering` | Degree-oriented triangle enumeration, clustering coefficients, transitivity, sampled mean local clustering, and k-truss. |
| `functions` | CSR graph inspection, typed attribute maps, and stable-ID borrowed subgraphs. |
| `planarity` | Linear-time left–right testing, rotation systems, face boundaries, and Kuratowski subdivisions. |
| `isomorphism` | VF2++ whole-graph and vertex-induced subgraph isomorphism, optional node labels, and lazy enumeration of mappings. |
| `community` | Active-queue label propagation, weighted undirected Leiden, and weighted Louvain and modularity for both graph kinds. |
| `generators` | Classical and named graphs, triangular and hexagonal lattices, expanders, graph lifts, algebraic graphs, tournaments, seeded random network models, bipartite regular multigraphs, and binomial intersection graphs. |
| `structures` | Union–find, rollback union–find, link–cut forests with path minima and path updates, indexed Fibonacci heaps and radix queues, bounded-arity trees, segment trees, static range minima, and Fenwick frequency trees. |

Nonnegative routing uses `u64` edge weights and `u128` distances. Signed routing uses `i64` weights and `i128` distances. Single-commodity capacities are `u64`, with `u128` flow and cut aggregates. `randomized_maximum_matching` randomizes the Hopcroft–Karp warm start. `random_walk_maximum_matching` finds alternating paths by loop-erased random walks and uses a bounded retry policy on arbitrary inputs. Fractional multicommodity flow uses floating point and can report numerical accuracy failures for ill-conditioned inputs.

The library requires Rust 1.86 or newer. The default `multicommodity` feature uses the pure-Rust [microlp solver](https://docs.rs/microlp/0.8.0/microlp/). Build with `--no-default-features` for a standard-library-only crate.

```sh
cargo test --locked
cargo test --locked --no-default-features
```

## Flow solver selection

For an exact directed maximum flow, `flow::max_flow` runs highest-label push–relabel. `max_flow_st_planar` uses a dual shortest-path computation when the terminals can share a face. `max_flow_core_cycles` maintains routes through contracted forests, then recovers an integer flow and a minimum cut. `max_flow_tree_cycle` can first run signed barrier descent with `TreeCycleBasis::DynamicHierarchy`, then continue from that flow through the same exact solver. Both use only the standard library. The general almost-linear running-time bound remains open in this implementation. [Sparse and dense benchmarks](benches/flow.md)

For repeated minimum-ratio-cycle queries on an undirected graph with signed integer costs and positive lengths, build an `EmbeddedCycleHierarchy` and apply batches with `update`. The inexpensive `improving_cycle` returns an embedded witness. Use `certified_cycle(factor)` for a global guarantee, or `certified_cycle(1)` for the exact minimum improving ratio.

For supplies with signed arc costs, `minimum_cost_flow_capacity_scaling` works through successively finer capacity thresholds. `minimum_cost_flow_minimum_mean_cycle` repeatedly cancels a residual cycle of minimum mean cost. Both accept positive supplies at senders, negative supplies at receivers, and zero supplies for minimum-cost circulation. An infeasible balance vector returns `Ok(None)`. Use `minimum_cost_flow_with_lower_bounds` when some edges must carry a minimum amount.

```rust
use graph_library::{Graph, GraphKind};
use graph_library::flow::{ArcCost, minimum_cost_flow_capacity_scaling};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let mut network = Graph::builder(2, GraphKind::Directed);
network.add_edge(0, 1, ArcCost { capacity: 5, unit_cost: -2 })?;
let optimum = minimum_cost_flow_capacity_scaling(&network.build(), &[3, -3])?.unwrap();
assert_eq!(optimum.edge_flows, [3]);
assert_eq!(optimum.total_cost, -6);
# Ok(())
# }
```

## Choosing a solver or query index

For integral max flow, `max_flow_capacity_scaling` handles a large spread of capacities through thresholded residual searches. Plain `max_flow_ford_fulkerson` can take one augmentation per unit of flow. `max_flow_network_simplex` converts max flow to a zero-supply circulation and returns the same per-edge flows and minimum-cut certificate as the other max-flow APIs. Use `network_simplex` directly for supplies and costs.

For many connectivity queries on one snapshot, build weak or strong components once and query their labels in O(1). For one-way reachability, `ReachabilityIndex` uses linear space at fixed dimension and searches only when its intervals leave a query unresolved. `ChainReachabilityIndex` offers O(1) queries when the condensation has a small path cover. Its constructor accepts a byte budget for the label table. For undirected queries under a changing edge-weight limit, `ThresholdConnectivityIndex` answers in O(1) after Kruskal preprocessing.

```rust
use graph_library::{Graph, GraphKind};
use graph_library::connectivity::{ChainReachabilityIndex, ReachabilityIndex};

# fn main() -> Result<(), graph_library::GraphError> {
let mut b = Graph::builder(5, GraphKind::Directed);
b.add_path(&[0, 1, 2, 3], ())?;
b.add_edge(2, 1, ())?;
let g = b.build();
let intervals = ReachabilityIndex::new(&g, 3, 42);
assert!(intervals.reachable(0, 3)?);
assert!(!intervals.reachable(3, 0)?);
let chains = ChainReachabilityIndex::with_budget(&g, 4096).unwrap();
assert!(chains.reachable(2, 1)?);
assert!(!chains.reachable(0, 4)?);
# Ok(())
# }
```

The [algorithm survey](docs/algorithm-survey.md) covers implementation provenance, primary sources, method selection, and the workload measurements for these indexes.

## Signed paths and bottleneck objectives

`bellman_ford` returns paths from one source and reports negative cycles reachable from that source. `johnson` and `floyd_warshall` return all-pairs distances and paths, and report a negative cycle anywhere in the graph. Cycle errors include the ordered vertices, edge IDs, and exact negative cost. `find_negative_cycle` checks for a witness without requesting an all-pairs matrix. `dag_shortest_paths` and `dag_longest_path` use a linear topological pass. Use `DagShortestPathIndex` to reuse the topological order across signed point-to-point and multi-source queries.

Bellman–Ford processes strongly connected components in topological order and relaxes only changed labels within a component. Johnson reweights the graph and reuses the indexed radix queue across sources. Floyd–Warshall uses cache tiles and breaks equal-distance ties by edge count so zero-cost cycles do not affect path reconstruction. Both all-pairs functions materialize O(V²) output.

```rust
use graph_library::{Graph, GraphKind};
use graph_library::shortest_path::{bellman_ford, johnson};

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let mut builder = Graph::builder(3, GraphKind::Directed);
builder.add_edge(0, 1, -4_i64)?;
builder.add_edge(1, 2, 1)?;
let graph = builder.build();
assert_eq!(bellman_ford(&graph, 0)?.distance(2), Some(-3));
assert_eq!(johnson(&graph)?.path(0, 2).unwrap().nodes, [0, 1, 2]);
# Ok(())
# }
```

`minimum_spanning_tree` minimizes total weight on an undirected graph. `minimum_bottleneck_path`, `minimum_bottleneck_spanning_tree`, and `minimum_bottleneck_arborescence` minimize the largest selected edge weight. For general matching, `minimum_bottleneck_matching` maximizes cardinality first. The exact-size and perfect variants require the stated cardinality.

The deterministic `duong_shortest_paths` function includes a bounded Bellman–Ford warm start, integer scaling, Duong’s local padded decomposition, and ordered merges. It returns the same exact path and witness types as `bellman_ford`. See the [Duong benchmark](benches/duong.md) for the measured tradeoffs and the [implementation notes](docs/negative-shortest-paths.md) for the distinction between the crate’s solver and the paper’s word-RAM analysis.

See the [signed-path benchmark](benches/signed_shortest_paths.md) and the [negative-edge shortest-path research survey](docs/negative-shortest-paths.md).

## General-graph matching

`maximum_general_matching` runs Edmonds' algorithm. `maximum_general_matching_gabow` uses Gabow's batched shortest augmenting paths and Gabow–Tarjan incremental tree set union. It uses O(V + E√V) time and O(V + E) space, including the vertex scan for graphs with isolates. Both functions return original edge IDs and ignore self-loops and weights.

The [matching benchmark](benches/README.md) compares the two functions on deterministic sparse, dense, and adversarial graph families.

```sh
cargo bench --locked --no-default-features --bench general_matching
```

### Weighted matching and assignment

`maximum_weight_matching(graph, false)` optimizes total weight. Pass `true` to maximize cardinality first. `minimum_cost_matching` minimizes the cost of a maximum-cardinality matching. `minimum_cost_perfect_matching` requires every vertex to be matched and returns `None` when this is infeasible. The specialized bipartite functions use a supplied partition and can require the smaller side to be fully matched.

Weights are signed `i64`, and totals use `i128`. The engines preserve original edge IDs, including when parallel edges have different weights. General weighted matching coalesces parallel edges, solves each component independently, and uses Gabow for uniform weights. The primal–dual solver records cyclic edges in contracted blossoms and opens them when their dual reaches zero or the solution is extracted. Its event queue stores predicted tight edges.

```rust
use graph_library::{Graph, GraphKind};
use graph_library::matching::{maximum_weight_matching, minimum_cost_perfect_matching};

# fn main() -> Result<(), graph_library::GraphError> {
let mut builder = Graph::builder(4, GraphKind::Undirected);
builder.add_edge(0, 1, -3_i64)?;
builder.add_edge(1, 2, 8)?;
builder.add_edge(2, 3, -2)?;
let graph = builder.build();

assert_eq!(maximum_weight_matching(&graph, false)?.total_weight, 8);
assert_eq!(maximum_weight_matching(&graph, true)?.matching.len(), 2);
assert_eq!(minimum_cost_perfect_matching(&graph)?.unwrap().total_weight, -5);
# Ok(())
# }
```

The [weighted benchmark](benches/weighted_matching.md) compares all three objectives with NetworkX. The [NetworkX coverage inventory](NETWORKX.md) lists supported algorithm families and gaps.

## Planarity, isomorphism, and communities

`check_planarity` returns the clockwise order of original edge IDs and the face boundary walks of the underlying simple undirected graph. `is_planar` skips embedding construction. `vf2pp_isomorphism` returns a left-to-right vertex bijection and accounts for edge direction, loops, and parallel edges. Optional labels and all-mapping iteration are available.

For a nonplanar graph, `kuratowski_witness` returns internally disjoint paths forming a subdivision of K₅ or K₃,₃, expressed with original vertex and edge IDs. `check_planarity_with_witness` returns either an embedding or such a subdivision. Extraction isolates a biconnected block, suppresses degree-two chains, and removes dispensable edges in batches. The resulting subdivision is inclusion-minimal.

`vf2pp_subgraph_isomorphism(host, pattern)` returns a vector from pattern vertices to host vertices. It preserves every edge and nonedge inside the image, including direction, loops, and multiplicity. Host edges leaving the image are unrestricted. Labelled and lazy all-mapping variants are available.

Louvain accepts finite, nonnegative weights. Its result contains original-vertex community IDs and the modularity at each aggregation level. The seed controls vertex order, so different seeds and implementations may find different partitions.

```rust
use graph_library::community::{LouvainOptions, louvain_communities_by};
use graph_library::generators::{caveman_graph, hexagonal_lattice_graph};
use graph_library::isomorphism::vf2pp_is_isomorphic;
use graph_library::planarity::check_planarity;

# fn main() -> Result<(), Box<dyn std::error::Error>> {
let lattice = hexagonal_lattice_graph(3, 4, false)?;
assert!(check_planarity(&lattice).is_some());
assert!(vf2pp_is_isomorphic(&lattice, &lattice));

let graph = caveman_graph(3, 5)?;
let result = louvain_communities_by(&graph, |_| 1.0, LouvainOptions::default())?;
assert_eq!(result.community_count(), 3);
# Ok(())
# }
```

The [network structure benchmark](benches/network_structure.md) covers planarity, whole-graph VF2++, and Louvain. The [witness and induced-match benchmark](benches/obstructions.md) compares Kuratowski extraction and induced matching on sparse, dense, and adversarial inputs.

## Algorithm references

- [Ahuja, Mehlhorn, Orlin, and Tarjan on radix heaps](https://www.princeton.edu/~alaink/Orf467F11/AhujaRadixHeap.pdf)
- [Duong’s deterministic negative-weight shortest-path algorithm](https://arxiv.org/html/2609.05590v2)
- [Johnson’s sparse shortest-path algorithm](https://doi.org/10.1145/321992.321993)
- [Camerini’s minimum bottleneck spanning tree algorithm](https://doi.org/10.1016/0020-0190(78)90030-3)
- [Goldberg and Harrelson’s landmark routing algorithm](https://www.microsoft.com/en-us/research/publication/computing-the-shortest-path-a-search-meets-graph-theory/)
- [Hopcroft and Karp’s bipartite matching algorithm](https://doi.org/10.1137/0202019)
- [Edmonds’ blossom algorithm](https://doi.org/10.4153/CJM-1965-045-4)
- [Gabow’s phased cardinality matching algorithm](https://arxiv.org/abs/1703.03998)
- [Mehlhorn and Nobahari’s canonical-path formulation and chain benchmark](https://arxiv.org/abs/2603.22909)
- [Gabow and Tarjan’s linear-time incremental tree set union](https://doi.org/10.1016/0022-0000(85)90014-5)
- [Batagelj and Zaveršnik’s core decomposition](https://arxiv.org/abs/cs/0310049)
- [Brandes’s left–right planarity test](https://acm.math.spbu.ru/~sk1/courses/1819s_au3/conspect/papers/brandes2010-planarity.pdf)
- [Jüttner and Madarasi’s VF2++ algorithm](https://egres.elte.hu/www/tr-18-03.html)
- [Blondel, Guillaume, Lambiotte, and Lefebvre’s Louvain algorithm](https://arxiv.org/abs/0803.0476)
- [NetworkX algorithm reference](https://networkx.org/documentation/stable/reference/algorithms/index.html)
- [Goldberg and Tarjan’s push–relabel algorithm](https://doi.org/10.1145/48014.61051)
- [Goldberg and Tarjan’s cost-scaling minimum-cost circulation algorithm](https://doi.org/10.1287/moor.15.3.430)
- [Stoer and Wagner’s global minimum cut algorithm](https://www.cs.dartmouth.edu/~ac/Teach/CS105-Winter05/Handouts/stoerwagner-mincut.pdf)
- [Gusfield’s all-pairs cut tree algorithm](https://doi.org/10.1137/0219009)
