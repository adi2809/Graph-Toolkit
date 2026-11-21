# NetworkX and Boost.Graph coverage

The Rust APIs use immutable CSR graphs with dense vertex IDs and typed edge payloads. User-defined edge and vertex attributes can live in separate typed maps. The algorithm, function, and generator references cover substantially more ground than this crate. [Algorithms](https://networkx.org/documentation/stable/reference/algorithms/index.html) · [Functions](https://networkx.org/documentation/stable/reference/functions.html) · [Generators](https://networkx.org/documentation/stable/reference/generators.html)

| Area | Implemented | Remaining examples |
| --- | --- | --- |
| Matching and covering | General weighted and unweighted matching, cardinality-first and bottleneck objectives, bipartite assignment and vertex covers, edge covers, greedy maximal matching, and validation. | Bipartite projections, extendability, and floating-point matching costs. |
| Flow and connectivity | Ford–Fulkerson, capacity scaling, push–relabel, Dinic, Edmonds–Karp, Boykov–Kolmogorov, cost and capacity scaling for minimum-cost flow, network simplex, minimum-mean-cycle cancellation, cofacial planar flow, dynamic-tree blocking flow, successive shortest path, residual cut certificates, global and all-pairs cuts, concurrent flow, local and global edge connectivity, terminal vertex cuts, disjoint paths, bidirectional and indexed reachability, and offline dynamic connectivity. | Global node connectivity, all-pairs disjoint-path queries, and more cut statistics. |
| Traversal and shortest paths | BFS, DFS, multi-source and four-way jump-point grid routing, 0–1 BFS, Dijkstra, bidirectional Dijkstra, A*, ALT, Bellman–Ford, Duong’s local decomposition solver, Johnson, Floyd–Warshall, DAG routing, minimax and k-shortest simple paths. | All-simple-path iterators, traversal classifications, time-dependent routing, and additional shortest-path variants. |
| Graph structure | Strong, weak, and biconnected components, condensation, bridges, articulation points, core numbers, topological orders and generations, DAG closure and reduction, dominators and dominance frontiers, Eulerian traversals, fundamental cycle bases, spanning trees and forests, bottleneck trees and arborescences, left–right planarity, planar rotation systems, Kuratowski subdivisions, and labelled VF2++ whole-graph and induced subgraph isomorphism. | Directed minimum-total-cost arborescences, minimum cycle bases, chordality, planar drawing, cliques, non-induced subgraph matching, graph polynomials, and additional DAG analyses. |
| Network analysis | PageRank, degree, harmonic, closeness, unweighted Brandes betweenness, triangles, local clustering, transitivity, k-truss, weighted undirected Leiden, chordal recognition and optimal chordal coloring, and weighted directed and undirected Louvain. | Additional community detection methods, spectral and current-flow centralities, link prediction, assortativity, similarity, sparsifiers, and approximation families. |
| Graph functions | Degree and density, nodes and edges, neighbor and non-edge iterators, self-loops, typed attributes, walk checks and weights, graph direction conversion, empty copies, borrowed subgraphs, and builder path, star, and cycle insertion. | Python-compatible mutation and attribute dictionaries. Broader interchange, file formats, drawing, and conversion APIs. |
| Generators | Classic paths, cycles, complete and multipartite graphs, wheels, stars, r-ary and binomial trees, barbells, lollipops, tadpoles, ladders, grids and hypercubes, circulants, generalized Petersen and several named graphs, clique communities, line and ego graphs, triangular and hexagonal lattices, LCF and named graphs, Margulis–Gabber–Galil expanders, random k-lifts, directed and clustered configuration, degree-sequence trees, small-world models, Kneser and Paley graphs, chordal-cycle and Harary graphs, random tournaments, growing directed networks, bipartite random graphs, and binomial intersection models. | Atlas and dataset generators, directed scale-free models, inverse line graphs, additional expanders, and further geometric models. |

The [algorithm and provenance survey](docs/algorithm-survey.md) describes the independent weighted matcher, named graph constructions, additional solvers, and large-graph indexes. NetworkX comparisons use an optional external Python installation.

## Graph function conventions

- `GraphBuilder::build` freezes structure. `Graph::is_frozen` reports that guarantee. Use `AttributeMap<T>` for additional typed data.
- `Graph::to_directed` and `Graph::to_undirected` create owned snapshots and preserve the individual payloads of parallel and antiparallel arcs.
- `induced_subgraph`, `restricted_view`, and `edge_subgraph` borrow the original graph and keep its IDs. Materialization returns the ID maps required to apply algorithms to the compact copy.
- Signed `path_weight` takes the cheapest parallel edge for each hop. Eulerian traversal operates on the subgraph of edge-bearing vertices. Unweighted betweenness counts parallel edges as distinct shortest-path choices.

## Structure algorithm conventions

- Planarity ignores direction, removes loops, and keeps the first edge of each unordered pair in the rotation system. A disconnected embedding represents outer boundaries separately per component. Kuratowski witnesses return original IDs for the paths of a K₅ or K₃,₃ subdivision.
- VF2++ preserves direction, loops, multiplicity, and optional vertex labels. Edge payloads are outside the equivalence relation. Induced matches use a pattern-to-host injection and require exact adjacency within its image.
- Louvain uses a seeded vertex order and nonnegative edge weights. The differential tests compare the modularity of the returned partition independently of NetworkX’s own heuristic result.

## Linked Boost.Graph overviews

| Overview | Rust entry points |
| --- | --- |
| [Topological sort](https://www.boost.org/doc/libs/latest/libs/graph/doc/html/graph/algorithms/topological_sort/topological_sort_overview.html) | `reverse_topological_sort` emits the iterative DFS finish order used by Boost. `topological_sort` emits source-first Kahn order. Both reject cycles. |
| [Connectivity](https://www.boost.org/doc/libs/latest/libs/graph/doc/html/graph/algorithms/connectivity/connectivity_overview.html) | `st_connected`, `connected_components`, `edge_connectivity`, `minimum_st_edge_cut`, and `DisjointSet`. |
| [Network flow](https://www.boost.org/doc/libs/latest/libs/graph/doc/html/graph/algorithms/network_flow/network_flow_overview.html) | `max_flow`, `max_flow_edmonds_karp`, `max_flow_boykov_kolmogorov`, `global_min_cut`, `minimum_cost_flow_cycle_canceling`, `min_cost_max_flow_successive_shortest_path`, and `flow_cost`. The earlier general matching solvers cover the two matching entries. The Rust flow API also offers Dinic, Goldberg–Tarjan cost scaling, and primal `network_simplex`. |
