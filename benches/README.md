# General-graph matching benchmark

Run the optimized benchmark with Cargo.

```sh
cargo bench --locked --no-default-features --bench general_matching
```

The harness builds each graph once, warms both matchers, then records seven timed runs per matcher with alternating order. Times include the matcher and its result allocation. The harness validates cardinality, vertex uniqueness, and original edge IDs outside the timed interval. It checks dense and adversarial results against the perfect-matching upper bound.

Sparse graphs sample 3V distinct undirected edges. Dense graphs include each possible edge with probability one third. The adversarial family attaches many short chains and increasingly long chains to a clique, with vertex IDs chosen so the greedy pass leaves two free vertices per chain. This follows the short-and-long-chain construction in the published Gabow benchmark. [Mehlhorn and Nobahari, Section 5](https://arxiv.org/html/2603.22909v2#S5)

## Measured medians

Rust 1.100.0 on arm64 macOS 27.0.1. Ratios divide Edmonds time by Gabow time.

| Family | Vertices | Edges | Matching size | Edmonds (ms) | Gabow (ms) | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Sparse | 6,000 | 18,000 | 2,995 | 1.153 | 3.312 | 0.35× |
| Sparse | 24,000 | 72,000 | 11,964 | 7.510 | 18.386 | 0.41× |
| Sparse | 60,000 | 180,000 | 29,921 | 40.085 | 89.401 | 0.45× |
| Dense | 400 | 26,402 | 200 | 0.028 | 0.020 | 1.37× |
| Dense | 800 | 106,347 | 400 | 0.106 | 0.071 | 1.49× |
| Dense | 1,600 | 426,702 | 800 | 0.298 | 0.257 | 1.16× |
| Short and long chains | 2,370 | 6,524 | 1,185 | 4.603 | 6.179 | 0.74× |
| Short and long chains | 7,240 | 20,370 | 3,620 | 39.810 | 34.736 | 1.15× |
| Short and long chains | 19,086 | 53,921 | 9,543 | 262.492 | 143.252 | 1.83× |
| Short and long chains | 40,002 | 113,789 | 20,001 | 1,133.209 | 449.980 | 2.52× |

Edmonds led on all three sparse random inputs. Gabow led on the sampled dense inputs by 0.041 ms or less. The largest gap appeared on the adversarial input, where Gabow completed in about 0.45 s and Edmonds took 1.13 s.

The [Duong benchmark](duong.md) compares deterministic signed shortest paths with the SCC-aware Bellman–Ford baseline.

The [network structure benchmark](network_structure.md) compares planarity, VF2++, and Louvain with NetworkX.
