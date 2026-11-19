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
| `traversal` | Insertion-ordered BFS and iterative DFS. Multi-source BFS returns a shortest-path forest with source ownership. |
| `grid` | Implicit grids with bit-packed obstacles, four- or eight-way movement, multi-source distance fields, and A*. |
| `shortest_path` | Indexed-radix Dijkstra, bidirectional Dijkstra, A*, ALT, 0–1 BFS, signed Bellman–Ford, Duong’s local padded solver, Johnson, blocked Floyd–Warshall, indexed DAG routing, minimax paths, and k-shortest simple paths. |
| `matching` | Edmonds and Gabow cardinality matching, exact weighted blossom, bottleneck matching, sparse bipartite assignment, edge covers, and validation. Deterministic and seeded Hopcroft–Karp also return minimum vertex covers. |
| `flow` | Push–relabel, Dinic, Edmonds–Karp, Boykov–Kolmogorov, cost scaling, network simplex, cycle canceling, successive shortest paths, Stoer–Wagner, and Gomory–Hu. |
| `multicommodity` | Maximum concurrent fractional flow for directed graphs. A sparse LP optimizes the common demand multiplier and returns a capacity- and conservation-checked allocation. |
| `connectivity` | Union–find, offline dynamic connectivity, components, edge cuts, dominators, Eulerian traversals, cycle bases, core numbers, DAG ordering and reduction, spanning trees and forests, and bottleneck arborescences. |
| `centrality` | Weighted personalized PageRank, degree, harmonic, closeness, and unweighted Brandes betweenness. |
| `clustering` | Degree-oriented triangle enumeration, clustering coefficients, transitivity, and k-truss. |
| `functions` | CSR graph inspection, typed attribute maps, and stable-ID borrowed subgraphs. |
| `planarity` | Linear-time left–right testing, rotation systems, and face boundaries. |
| `isomorphism` | VF2++ directed and undirected multigraph isomorphism, optional node labels, and lazy enumeration of mappings. |
| `community` | Seeded weighted Louvain community detection and modularity for directed and undirected graphs. |
| `generators` | Classical and named graphs, triangular and hexagonal lattices, expanders, graph lifts, and seeded random network models. |

Nonnegative routing uses `u64` edge weights and `u128` distances. Signed routing uses `i64` weights and `i128` distances. Single-commodity capacities are `u64`, with `u128` flow and cut aggregates. The randomized bipartite solver uses its seed to choose the greedy warm start and then runs exact augmenting phases. Fractional multicommodity flow uses floating point and can report numerical accuracy failures for ill-conditioned inputs.

The library requires Rust 1.86 or newer. The default `multicommodity` feature uses the pure-Rust [microlp solver](https://docs.rs/microlp/0.8.0/microlp/). Build with `--no-default-features` for a standard-library-only crate.

```sh
cargo test --locked
cargo test --locked --no-default-features
```

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

Weights are signed `i64`, and totals use `i128`. The engines preserve original edge IDs, including when parallel edges have different weights. General weighted matching uses a primal–dual blossom algorithm with O(V³ + E) worst-case time after coalescing parallel edges, solves independent components separately, and uses the Gabow matcher for uniform-weight components.

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

See the [network structure benchmark](benches/network_structure.md) for sparse, dense, and adversarial inputs and the NetworkX comparisons.

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
