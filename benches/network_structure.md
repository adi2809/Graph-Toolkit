# Planarity, VF2++, and Louvain benchmarks

Measured on arm64 macOS with rustc 1.100.0 and NetworkX 3.2.1. Rust figures are medians of seven optimized runs. NetworkX figures are medians of three runs. Both runtimes receive a warm-up call. Timings include solver preprocessing and result construction. Graph generation, input conversion to CSR or NetworkX, and result verification are outside the timers.

## Planarity with embeddings

Both engines build a rotation system on planar inputs. The complete graphs exercise rejection by the planar edge bound. The Rust projection can stop as soon as it has counted enough distinct pairs. Consecutive insertion into triangular faces gives the stacked triangulations a deep DFS and exactly `3V - 6` edges, forcing a full traversal and embedding.

| Family | Vertices | Edges | Rust (ms) | NetworkX (ms) | Planar |
| --- | ---: | ---: | ---: | ---: | --- |
| Sparse grid | 4,096 | 8,064 | 0.819 | 222.668 | Yes |
| Sparse grid | 16,384 | 32,512 | 3.269 | 1,101.420 | Yes |
| Dense complete | 128 | 8,128 | 0.002 | 6.320 | No |
| Dense complete | 512 | 130,816 | 0.008 | 101.013 | No |
| Worst-case stacked triangulation | 4,096 | 12,282 | 1.042 | 316.034 | Yes |
| Worst-case stacked triangulation | 16,384 | 49,146 | 4.220 | 1,505.349 | Yes |

## VF2++ whole-graph isomorphism

Sparse and dense positive cases use seeded random graphs and independently relabel their vertices and reverse the insertion order of their edges. The adversarial negative pair is the 4×4 rook graph and the Shrikhande graph, which have the same strongly regular parameters and defeat degree and one-dimensional color refinement alone. The larger pair takes their Cartesian products with a triangle. These hard negative cases exercise the exact search.

| Family | Vertices per graph | Edges per graph | Rust (ms) | NetworkX (ms) | Isomorphic |
| --- | ---: | ---: | ---: | ---: | --- |
| Sparse undirected | 1,024 | 4,096 | 1.212 | 80.878 | Yes |
| Sparse directed | 8,192 | 32,768 | 16.122 | 6,668.418 | Yes |
| Dense undirected | 128 | 2,709 | 0.168 | 9.818 | Yes |
| Dense directed | 512 | 43,605 | 2.615 | 218.799 | Yes |
| Adversarial strongly regular pair | 16 | 48 | 0.039 | 4.394 | No |
| Adversarial strongly regular product | 48 | 192 | 0.165 | 16.632 | No |

## Louvain communities

The seeded stochastic block graphs use higher within-block connection probabilities. The rings have a single link between consecutive cliques and exercise repeated aggregation across many near-tied local moves. Both solvers use resolution 1, threshold `1e-7`, and seed 613. Each engine has its own vertex shuffle and tie breaking, so the seed reproduces runs within each engine.

| Family | Vertices | Edges | Rust (ms) | NetworkX (ms) | Rust modularity | NetworkX modularity | Rust levels |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Sparse blocks | 2,048 | 6,239 | 4.254 | 182.932 | 0.740616 | 0.745516 | 3 |
| Sparse directed blocks | 8,192 | 120,737 | 26.770 | 1,513.411 | 0.833558 | 0.833558 | 2 |
| Dense blocks | 256 | 13,846 | 1.266 | 69.793 | 0.217350 | 0.217350 | 1 |
| Dense directed blocks | 512 | 110,949 | 10.659 | 517.177 | 0.218712 | 0.218712 | 1 |
| Adversarial ring of cliques | 960 | 2,112 | 0.329 | 30.147 | 0.955418 | 0.955705 | 3 |
| Adversarial ring of cliques | 5,120 | 11,264 | 1.969 | 168.597 | 0.980699 | 0.980784 | 4 |

NetworkX independently scored every Rust partition and agreed with the reported score to within `1e-10` across the benchmarks. The largest difference between the two independently chosen partitions was 0.00490 on the 2,048-vertex sparse block input, where NetworkX found the higher score.

## Correctness checks and reproduction

The optional oracle checks all 1,253 NetworkX Graph Atlas graphs and another 180 directed and undirected multigraphs for planarity. Rust verifies both dart coverage and Euler’s relation for each connected planar component. The differential suite also runs 210 VF2++ comparisons, 55 independent Louvain modularity checks, and 125 named-graph and lattice comparisons. Native tests cover deep DFS and VF2++ stacks, labelled and unlabelled automorphism counts, weighted and parallel edges, and seeded generator invariants.

```sh
cargo bench --locked --no-default-features --bench network_structure -- --networkx
cargo test --release --locked --test network_structure -- --ignored --nocapture
```

The benchmark runs Rust alone when the `--networkx` flag is omitted. Set `NETWORKX_PYTHON` to select a Python executable with NetworkX installed.
