# Deterministic signed shortest paths benchmark

`duong_shortest_paths` and the crate’s SCC-aware `bellman_ford` run from vertex zero on the same input graphs. These are medians of seven optimized Rust runs on arm64 macOS using rustc 1.100.0. One warm-up and rotating algorithm order precede the measurements. The timer covers each public solver call, including its setup and output allocation. Graph generation and the check against Bellman–Ford are outside the timing window.

| Graph family | Vertices | Edges | Bellman–Ford (ms) | Duong (ms) |
| --- | ---: | ---: | ---: | ---: |
| Sparse | 1,024 | 4,096 | 0.087 | 0.087 |
| Sparse | 8,192 | 32,768 | 1.411 | 1.419 |
| Sparse | 32,768 | 131,072 | 8.509 | 8.550 |
| Dense | 128 | 10,838 | 0.146 | 0.146 |
| Dense | 256 | 43,390 | 0.498 | 0.500 |
| Dense | 512 | 174,556 | 1.651 | 1.665 |
| Worst-case chain and clique | 256 | 16,512 | 3.625 | 5.259 |
| Worst-case chain and clique | 512 | 65,792 | 28.929 | 22.581 |
| Worst-case chain and clique | 1,024 | 262,656 | 230.544 | 104.542 |

The random sparse and dense inputs assign feasible integer potentials to create signed edges with no negative cycle. The adversarial input is strongly connected. Each improvement along its negative chain triggers another scan of the zero-weight clique in SCC-aware Bellman–Ford, giving Θ(VE) work on this family. The Duong timings include the bounded Bellman–Ford warm start. The largest adversarial run was about 2.21× faster. The smallest adversarial run still favored Bellman–Ford.

The optional integration test compares exact distances and reachable-cycle detection with NetworkX 3.2.1 on 273 graphs, including two 81-vertex strongly connected cases long enough to exceed the warm-start limit. It also checks the returned paths and cycle witnesses against the original edge IDs and signed costs. A separate unit test compares 220 graphs directly against Bellman–Ford with the warm start disabled.

```sh
cargo bench --locked --no-default-features --bench signed_shortest_paths -- --duong
cargo test --release --locked --test signed_shortest_paths -- --ignored --nocapture
```

See the [implementation notes](../docs/negative-shortest-paths.md) for the local padded decomposition, the integer queue, and the arithmetic fallback.
