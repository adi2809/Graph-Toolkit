# Weighted matching benchmark

Run both implementations on the same deterministic inputs with Python and NetworkX installed.

```sh
cargo bench --locked --no-default-features --bench weighted_matching -- --networkx
cargo test --locked --test weighted_matching -- --ignored --nocapture
```

The benchmark builds the graphs before timing. Each Rust median uses seven runs, and each NetworkX median uses three. Both engines receive a warm-up run for each objective. Medians include each matching call and its result allocation. Serialization, graph construction, and result validation are outside the timed intervals. All 27 objective values agreed. Cardinalities also agreed for the two objectives that optimize cardinality first.

The sparse family is a connected graph with 3V distinct edges. The dense family has roughly two thirds of the possible edges. Both use weights from −1,000 to 1,000. The worst-case family is a complete graph with two ladders of disjoint, distinctly weighted edges and unit-weight background edges. One ladder is strongly positive and the other strongly negative. Each ladder forces a linear number of separate augmentations for its favored objective, while the remaining free vertices still have quadratic-size adjacency to scan. This stresses repeated augmentation on a dense graph.

## Independent solver

The replacement solver uses a primal–dual quotient graph and an indexed event heap. The table gives seven-run Rust medians from November 20, 2026 on the same fixtures. The earlier implementation was faster on these fixtures. For example, dense 256-vertex cardinality-first matching rose from 18.743 ms to 44.061 ms. The exhaustive objective tests and optional Python comparison cover the replacement independently.

| Family | Vertices | Maximum weight (ms) | Cardinality first (ms) | Minimum cost (ms) |
| --- | ---: | ---: | ---: | ---: |
| sparse | 128 | 0.599 | 0.947 | 0.981 |
| sparse | 256 | 2.109 | 3.456 | 4.058 |
| sparse | 512 | 9.391 | 12.963 | 12.962 |
| dense | 64 | 0.750 | 1.435 | 1.338 |
| dense | 128 | 4.671 | 8.677 | 7.989 |
| dense | 256 | 23.497 | 44.061 | 39.767 |
| worst-case ladders | 64 | 0.888 | 0.888 | 0.899 |
| worst-case ladders | 128 | 6.164 | 6.198 | 6.308 |
| worst-case ladders | 256 | 46.633 | 46.732 | 48.301 |

## Historical measurements before replacement

Rust 1.100.0 and Python 3.9.6 with NetworkX 3.2.1 on arm64 macOS 27.0.1. All times are in milliseconds. Each cell lists Rust / NetworkX.

| Family | Vertices | Edges | Maximum weight | Cardinality first | Minimum cost |
| --- | ---: | ---: | ---: | ---: | ---: |
| Sparse | 128 | 384 | 0.371 / 56.520 | 0.610 / 62.427 | 0.589 / 62.626 |
| Sparse | 256 | 768 | 1.098 / 197.300 | 2.070 / 224.703 | 2.506 / 259.704 |
| Sparse | 512 | 1,536 | 4.059 / 763.298 | 8.023 / 918.233 | 7.260 / 878.934 |
| Dense | 64 | 1,355 | 0.363 / 85.518 | 0.632 / 85.688 | 0.586 / 82.771 |
| Dense | 128 | 5,491 | 2.128 / 519.313 | 3.512 / 521.013 | 3.432 / 505.036 |
| Dense | 256 | 21,683 | 10.488 / 2,997.697 | 18.743 / 2,994.181 | 17.270 / 2,764.998 |
| Worst-case ladder | 64 | 2,016 | 0.436 / 84.212 | 0.442 / 84.278 | 0.436 / 84.904 |
| Worst-case ladder | 128 | 8,128 | 3.100 / 638.063 | 3.112 / 638.709 | 3.114 / 640.705 |
| Worst-case ladder | 256 | 32,640 | 22.466 / 5,835.777 | 22.582 / 5,753.955 | 21.945 / 5,420.135 |

The comparison measures the Rust and Python implementations on these fixtures. Integer optima can have multiple edge sets, so validation compares the objective and the cardinality where specified, then checks the reported original edges for validity.
