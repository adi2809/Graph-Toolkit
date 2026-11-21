# Maximum flow on sparse, dense and planar inputs

## Standalone core cycles and signed barrier descent

All four solvers returned the same certified flow value on every instance. Push–relabel had the lowest median in each row. The standalone core-cycle solver beat Dinic on the smaller sparse instance. The larger sparse instance had overlapping timing ranges. On both dense instances Dinic was faster. The signed barrier variant accepted eight steps on each input and took longer end to end. [Median, minimum and maximum timings](flow_core_cycles.csv)

| Family | Vertices | Arcs | Dinic ms | Push–relabel ms | Core cycles ms | Signed barrier + core cycles ms | Maximum flow |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Sparse | 16,384 | 98,431 | 22.037 | 6.792 | 16.374 | 129.108 | 35,163 |
| Sparse | 65,536 | 393,727 | 168.152 | 10.943 | 165.462 | 759.551 | 132,511 |
| Dense | 512 | 130,236 | 2.886 | 1.351 | 17.509 | 114.241 | 118,588 |
| Dense | 1,024 | 523,531 | 19.673 | 9.415 | 92.361 | 499.239 | 245,988 |

The standalone solver ran until its rounding and cut certificates agreed. The signed run used the existing benchmark budget of eight barrier steps, then passed its feasible edge flows into the same standalone solver. The next table gives the standalone solver’s work and the signed variant’s integral throughput at the handoff.

| Family | Vertices | Terminal forest builds that routed flow | Augmentations | Fractional edges before rounding | Integral flow after signed barrier steps |
| --- | ---: | ---: | ---: | ---: | ---: |
| Sparse | 16,384 | 10 | 169 | 465 | 599 |
| Sparse | 65,536 | 20 | 689 | 1,852 | 599 |
| Dense | 512 | 2 | 561 | 826 | 518 |
| Dense | 1,024 | 3 | 1,157 | 1,643 | 718 |

## Static tree bases, dynamic-tree blocking flow and planar flow

These timings come from a separate run of the broader flow benchmark. The BFS and multiscale forest variants now continue from their integral warm starts through core-cycle optimization and dyadic recovery.

| Family | Vertices | Directed arcs | Dinic ms | Push–relabel ms | Dynamic-tree blocking flow ms | BFS cycles + core ms | Both static bases + core ms | Hassin planar ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Sparse | 16,384 | 98,431 | 22.896 | 6.978 | 44.919 | 70.540 | 162.180 | — |
| Sparse | 65,536 | 393,727 | 192.122 | 13.005 | 466.884 | 463.492 | 1067.571 | — |
| Dense | 512 | 130,236 | 2.718 | 1.354 | 3.779 | 62.421 | 125.915 | — |
| Dense | 1,024 | 523,531 | 18.809 | 9.481 | 29.597 | 268.204 | 549.750 | — |
| Bidirected grid | 16,384 | 65,024 | 0.998 | 0.636 | 1.860 | 136.163 | 174.583 | 13.482 |

Each hybrid timing includes its exact residual completion. The table below gives feasible integral throughput immediately before that completion. Both configurations accepted eight fractional descent steps in each run.

| Family | Vertices | BFS cycles | Both bases | Certified maximum |
| --- | ---: | ---: | ---: | ---: |
| Sparse | 16,384 | 3,364 | 599 | 35,163 |
| Sparse | 65,536 | 2,739 | 599 | 132,511 |
| Dense | 512 | 2,466 | 677 | 118,588 |
| Dense | 1,024 | 2,834 | 718 | 245,988 |
| Bidirected grid | 16,384 | 40 | 2 | 912 |

## Reproduction

Run either comparison on the shared seeded graph generators.

```sh
cargo bench --locked --no-default-features --bench flow_core_cycles
cargo bench --locked --no-default-features --bench flow
```

The runs used one solver thread, release optimization, Rust 1.100.0 and arm64 macOS 27.0.1. Each solver receives one warm-up and five measured runs in rotating order. Timings include building the solver's data structures, optimizing and extracting its result. The graph is built once outside the timer. Both harnesses print median, minimum and maximum milliseconds. The tree-cycle runs report the number of accepted barrier steps and their integral throughput. [Core-cycle comparison](flow_core_cycles.rs) · [Planar and static-basis comparison](flow.rs)

- Sparse instances have a directed backbone, five pseudorandom arcs per vertex on average, and arcs from the source and into the sink at regular intervals. Parallel edges and loops are allowed.
- Dense instances sample each ordered pair of distinct vertices at probability one half.
- The grid has 128 rows and columns with independent capacities in both directions. Source and sink occupy opposite corners. The planar solver was timed on this cofacial instance.
- Capacities range from 1 through 1,000. The generator restarts from the same fixed seed at each size. Repetitions measure timing variability on one instance per row.

After every run the harness checks all original edge capacities, conservation at every vertex, source and sink balance, exact residual reachability and equality between the cut capacity and the flow value. All methods agreed. See the [benchmark harness](flow.rs).

## Exactness and running time

`max_flow_core_cycles` augments through a hierarchy of contracted terminal forests at decreasing capacity thresholds. Pairing parallel core arcs retains the shortest currently eligible witness in each direction. All represented terminal routes have minimum hop count at that threshold. Every residual length bucket is between one and O(log(2+mU)) for maximum input capacity U, which bounds the route's approximation factor for the weighted terminal objective. Long routes are represented by one crossing and two forest paths. Heavy-path decomposition and a lazy segment tree update those paths and report edges when their eligibility or length bucket changes. Capacity scaling bounds the number of augmentations by O(m log(2+mU)). A forest rebuild or a metric repair can still take linear work.

For exact recovery, let D be the smallest power of two greater than the edge count. Shrinking each positive forward residual capacity by 1/D preserves a feasible integral warm start and lowers every cut capacity by less than one. Whenever the ceiling of the current feasible throughput meets a known original cut, the solver rounds the dyadic flow one bit at a time. At each bit it orients the odd-edge trails, taking the source-to-sink trail first if the terminals have odd degree. This rounds throughput upwards while preserving conservation and respecting the original capacities. A residual traversal supplies the final cut. Recovery takes O((V+E)(1+log D)) time and O(V+E) space. [Flow rounding literature](https://arxiv.org/abs/1507.08139) · [Recovery implementation](../src/flow/augment/recovery.rs)

`EmbeddedCycleHierarchy` handles arbitrary signed gradients and positive lengths on a fixed undirected topology. Each level contracts forest stars, combines parallel edges and embeds selected edges through two retained spokes. At construction the two-hop path costs at most four times the direct length and each spoke supports at most eight replacements per level. Signed closed walks lift to original edges. A dependency index propagates changed metrics, and a volume threshold amortizes topology rebuilds across distinct length changes. A forest potential gives a lower bound on every cycle ratio by cancellation of potential differences. When that bound is loose, `certified_cycle` makes Bellman–Ford threshold decisions to establish the requested approximation on the full graph. Each such decision can take O(VE). [Hierarchy implementation](../src/flow/core_cycles.rs)

`max_flow_tree_cycle` adds a return edge, computes barrier derivatives and backtracks on the actual change in potential. Its dynamic basis sends only changed edge metrics to the signed hierarchy. It certifies directions on the quantized instance with factor 1+floor(log₂ m′), where m′ counts the participating edges including the return edge. The BFS and multiscale static bases remain available. All three modes pass their feasible integral residual state to the terminal solver. The signed barrier stage may lower its potential without increasing that integral throughput. [Chen et al.](https://arxiv.org/abs/2203.00671) · [Alon, Karp, Peleg and West](https://dwest.web.illinois.edu/pubs/kserver.pdf)

The published deterministic m^(1+o(1)) bound for polynomially bounded integral data relies on a low-stretch forest collection, dynamic sparsification with bounded recourse and a shift-and-rebuild schedule across levels. This implementation can propagate a single edit through the whole hierarchy, uses full-graph decisions for some ratio certificates, and completes optimization via capacity scaling. The paper’s subpolynomial amortized update bound and end-to-end interior-point iteration bound remain unimplemented. [Van den Brand et al., Theorems 4.6–4.7](https://arxiv.org/pdf/2309.16629) · [Kyng et al.](https://arxiv.org/abs/2609.23852)

`max_flow_dynamic_trees` computes blocking flows with a splay-based link–cut forest. Tree edges carry residual capacities. A path-minimum query identifies a bottleneck, a lazy path update subtracts it, and cutting a tree edge writes the accumulated flow back to the residual graph. The solver takes O(V(E+V) log(V+2)) time including the per-phase vertex arrays. [Sleator and Tarjan](https://www.cs.cmu.edu/~sleator/papers/dynamic-trees.pdf)

`max_flow_st_planar` checks whether an extra return path through a new vertex admits a planar embedding. If so, shortest-path potentials on the dual yield a maximum flow on the original edges. Deterministic grouping of parallel edges, embedding, radix shortest paths on fixed-width integer keys, and certificate extraction take O((V+E) log(V+2)) time and O(V+E) space together. Noncofacial terminals or a nonplanar underlying graph give `Ok(None)`. [Hassin](https://www.tau.ac.il/~hassin/planarst.pdf)
