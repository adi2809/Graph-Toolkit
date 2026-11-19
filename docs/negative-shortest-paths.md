# Negative-edge shortest-path research

Surveyed on November 19, 2026. In the bounds below, V is the vertex count, E is the arc count, and W bounds the magnitude of negative integer edge weights. The results concern single-source shortest paths in directed graphs.

| Work | Weight model and guarantee | Technique |
| --- | --- | --- |
| Bernstein, Nanongkai, and Wulff-Nilsen, 2022 | Randomized O(E log⁸ V log W) for integer weights. | Directed low-diameter decomposition and scaling established the first combinatorial near-linear bound for general directed graphs. [Paper](https://arxiv.org/abs/2203.03456) |
| Bringmann, Cassis, and Fischer, 2023 | Randomized O(E log² V log(VW) log log V) for integer weights, including negative cycle finding. | Faster decomposition and a separate cycle-finding procedure remove most of the original logarithmic factors. [Paper](https://arxiv.org/abs/2304.05279) |
| Li, Mowry, and Rao, SODA 2026 | O((E + V log log V) log(VW) log V log log V) expected time for integer weights. | Faster randomized directed low-diameter decomposition. [Paper](https://arxiv.org/abs/2510.22721) |
| Haeupler, Jiang, and Saranurak, and independently Li, 2025–2026 | Deterministic near-linear time for integer weights. | The first line uses directed path covers. Li uses deterministic padded decompositions. [Path covers](https://arxiv.org/abs/2511.08551) · [Padded decompositions](https://arxiv.org/abs/2511.07859) |
| Hair, G. Z. Li, J. Li, and Zhang, July 2026 | Las Vegas randomized E^(1+o(1)) time with high probability for real-valued edge weights. | Shortcutting and sparsification reduce the number of negative hops while controlling the size of the auxiliary graphs. [Paper, Theorem 1.1](https://arxiv.org/pdf/2607.19346) |
| Duong, September 2026 | Deterministic O((E + V log log V) log² V log(VW)) for integer weights, with a negative-cycle witness. | Size-sensitive padded decompositions and BF-Dijkstra merging give a logarithmic improvement over Li’s previous deterministic bound. [Paper](https://arxiv.org/html/2609.05590v2) |

## Implementation evidence

Cassis, Karrenbauer, Nusser, and Rinaldi published a C++20 implementation of the 2023 algorithm and benchmarked it against Goldberg–Radzik. On their hardest instance family they report a gain of nearly two orders of magnitude. Their study also shows why SCC decomposition needs to be accounted for in a comparison. Some published Bellman–Ford worst cases are DAGs and become linear-time instances once that structure is exploited. [Algorithm engineering paper, Sections 5.1–5.5](https://arxiv.org/html/2502.11999v1#S5)

## Rust implementation

`bellman_ford` uses SCC decomposition and an active frontier. Its running time is linear on DAGs and O(V + VE) in the worst case. `duong_shortest_paths` first tries at most 48 improving rounds per SCC so easy inputs finish without scaling. It then works on the source-reachable subgraph. Each scaling step halves the largest negative reduced cost, solves overlapping padded pieces, and combines the child potentials through prefix and suffix transition chains. The nonpositive-arc condensation supplies a linear-time seed when its internal arcs all have zero cost.

The padded split uses Duong’s Lemma 4.9, selected by the paper when the ceiling of the log of graph mass is less than 256. Entering the size-sensitive peeling regime would require mass above `2^255`. The solver uses 128-bit radix queues in place of the theoretical integer queue. It uses a conservative arithmetic guard with an exact Bellman–Ford fallback on exceptionally large combinations of graph size and weight magnitude. The crate makes no claim to the paper’s unrestricted word-RAM bound.

A failed merge can report a longer walk of negative shifted cost. The solver closes the walk with a clipped return path when a weak-diameter bound is available, then extracts a negative simple cycle using original edge IDs. Exact distances follow from a feasible rounded potential and a final Dijkstra traversal. Verification covers both the public dispatch and the scaling path with the warm start disabled. [Duong’s algorithm and proof](https://arxiv.org/html/2609.05590v2)

See the [Duong benchmark](../benches/duong.md) for sparse, dense, and worst-case Bellman–Ford comparisons. Johnson’s all-pairs algorithm continues to use the SCC-aware Bellman–Ford potentials. Its explicit output matrix requires Ω(V²) time and space.
