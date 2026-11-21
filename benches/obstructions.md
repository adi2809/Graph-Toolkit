# Kuratowski witnesses and induced subgraph isomorphism

Timings on arm64 macOS with rustc 1.100.0 and NetworkX 3.2.1. The Rust columns are medians of seven optimized runs and the NetworkX columns are medians of three runs. Both engines receive a warm-up. The timer covers the public algorithm call through result construction. Input generation, conversion, and verification are outside the timer.

## Kuratowski extraction

NetworkX runs `check_planarity(counterexample=True)`. Its extractor probes individual edge deletions. Rust caps dense projections at the first prefix that violates the planar edge bound, checks biconnected blocks, and suppresses degree-two chains while preserving original paths. It tests batches of compressed edges for deletion and splits batches whose deletion makes the graph planar.

The sparse graphs connect multiple planar Tutte blocks to a final K₃,₃. Complete graphs exercise the dense prefix. In the adversarial diamond graphs, each K₃,₃ edge becomes a chain of diamonds, and each diamond has a chord between its two internal vertices. The chords prevent initial suppression, and a witness must cross a growing number of diamonds. Extraction has a quadratic worst case. The witnesses can differ in length because they are inclusion-minimal subdivisions.

| Family | Vertices | Edges | Rust (ms) | NetworkX (ms) | Witness edges, Rust / NetworkX |
| --- | ---: | ---: | ---: | ---: | ---: |
| Sparse connected blocks | 190 | 289 | 0.049 | 653.931 | 9 / 9 |
| Sparse connected blocks | 374 | 569 | 0.068 | 2,537.251 | 9 / 9 |
| Dense complete | 32 | 496 | 0.150 | 164.250 | 10 / 10 |
| Dense complete | 64 | 2,016 | 0.290 | 1,942.978 | 10 / 10 |
| Worst-case diamond subdivision | 105 | 180 | 3.211 | 572.030 | 75 / 72 |
| Worst-case diamond subdivision | 213 | 360 | 12.592 | 2,323.500 | 147 / 144 |

## VF2++ induced subgraph isomorphism

The NetworkX 3.2.1 baseline is the exact VF2 `GraphMatcher.subgraph_is_isomorphic` API. The benchmark harness selects the induced VF2++ API automatically on NetworkX versions that provide it. Rust returns a mapping on successful cases. The sparse and dense inputs contain planted induced patterns. The adversarial inputs use a rook graph and a Shrikhande pattern, plus their Cartesian products with a triangle and extra isolated host vertices. They match the pattern and host core on degree and edge count.

| Family | Host vertices | Host edges | Pattern vertices | Pattern edges | Rust (ms) | NetworkX (ms) | Found |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Sparse planted motif | 584 | 1,115 | 8 | 10 | 0.192 | 20.460 | Yes |
| Sparse planted motif | 1,608 | 3,131 | 8 | 10 | 0.559 | 57.434 | Yes |
| Dense planted pattern | 96 | 3,040 | 7 | 14 | 0.025 | 0.829 | Yes |
| Dense planted pattern | 192 | 12,224 | 8 | 19 | 0.272 | 1.787 | Yes |
| Adversarial strongly regular | 28 | 48 | 16 | 48 | 0.072 | 66.626 | No |
| Adversarial strongly regular product | 64 | 192 | 48 | 192 | 0.266 | 331.280 | No |

## Verification

The differential suite checks 1,363 planarity decisions, including the complete NetworkX Graph Atlas. For all 302 nonplanar graphs it independently suppresses the degree-two vertices of both witnesses and checks that the kernels are K₅ or K₃,₃. Rust also verifies original IDs, path coverage, and disjoint interiors, including a 36,864-edge subdivision. The induced matcher is compared with NetworkX on 300 graph pairs, including labels, directed edges, loops, and parallel edges. A separate native test compares complete mapping sets against exhaustive injections.

```sh
cargo bench --locked --no-default-features --bench obstructions -- --networkx
cargo test --release --locked --test obstructions -- --ignored --nocapture
```

Omit the benchmark flag to run Rust alone. Set `NETWORKX_PYTHON` to select a Python executable with NetworkX installed.
