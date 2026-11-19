# Signed shortest-path benchmark

Run all three algorithms against the same graphs with Python and NetworkX installed.

```sh
cargo bench --locked --no-default-features --bench signed_shortest_paths -- --networkx
cargo test --locked --test signed_shortest_paths -- --ignored --nocapture
```

The harness validates distances after timing. Rust reports medians from seven runs, and NetworkX reports medians from three. Graph construction, serialization, and comparison are outside the timed intervals. Bellman–Ford computes one source and both other algorithms compute all pairs. NetworkX uses its predecessor-and-distance implementations for Bellman–Ford and Floyd–Warshall. Its Johnson API materializes complete vertex paths, whose costs the harness sums outside the timed interval. Rust returns distances with one predecessor per reachable destination.

The sparse input contains 4V distinct directed edges, including a directed ring. The dense input samples roughly two thirds of all directed non-loop edges. Both generate signed weights by adding potential differences to nonnegative costs, so cycle costs remain nonnegative. The worst-case family is a single strongly connected component with a long negative chain feeding a zero-cost directed clique. Every improvement along the chain forces another scan of the clique, exercising Θ(VE) work in active-frontier Bellman–Ford.

## Measured medians

Rust 1.100.0 and Python 3.9.6 with NetworkX 3.2.1 on arm64 macOS 27.0.1. All times are in milliseconds. Each cell lists Rust / NetworkX.

| Family | Vertices | Edges | Bellman–Ford | Johnson | Floyd–Warshall |
| --- | ---: | ---: | ---: | ---: | ---: |
| Sparse | 96 | 384 | 0.014 / 0.340 | 0.617 / 28.282 | 1.961 / 93.069 |
| Sparse | 192 | 768 | 0.036 / 0.835 | 2.725 / 113.529 | 14.016 / 708.014 |
| Sparse | 384 | 1,536 | 0.073 / 1.613 | 9.853 / 581.444 | 98.643 / 5,694.881 |
| Dense | 64 | 2,698 | 0.040 / 2.070 | 0.415 / 80.504 | 0.775 / 31.033 |
| Dense | 128 | 10,838 | 0.159 / 8.119 | 3.514 / 602.296 | 7.307 / 228.348 |
| Dense | 224 | 33,195 | 0.411 / 16.855 | 15.777 / 3,173.982 | 38.976 / 1,160.845 |
| Worst-case chain and clique | 64 | 1,056 | 0.062 / 0.975 | 0.204 / 30.255 | 0.414 / 24.633 |
| Worst-case chain and clique | 128 | 4,160 | 0.469 / 3.613 | 1.411 / 225.258 | 3.199 / 186.331 |
| Worst-case chain and clique | 192 | 9,312 | 1.604 / 7.984 | 4.316 / 744.002 | 10.626 / 608.197 |

All 27 algorithm-and-input combinations agreed with NetworkX. The separate oracle suite verifies distance matrices and negative-cycle status on 271 graphs, including parallel edges, loops, disconnected graphs, zero-cost cycles, and signed integer extremes.
