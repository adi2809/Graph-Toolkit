"""NetworkX oracle for graph structure and modularity benchmarks."""

from itertools import islice
import statistics
import sys
import time

import networkx as nx


def graph(tokens):
    n, m, directed = [int(next(tokens)) for _ in range(3)]
    edges = []
    pairs = set()
    for _ in range(m):
        u, v = [int(next(tokens)) for _ in range(2)]
        w = float(next(tokens))
        edges.append((u, v, w))
        pairs.add((u, v) if directed or u <= v else (v, u))
    multi = len(pairs) < m
    kind = (nx.MultiDiGraph if multi else nx.DiGraph) if directed else (
        nx.MultiGraph if multi else nx.Graph
    )
    result = kind()
    result.add_nodes_from(range(n))
    result.add_weighted_edges_from(edges)
    return result


def measure(fn, repetitions):
    if repetitions > 1:
        fn()
    durations = []
    for _ in range(repetitions):
        start = time.perf_counter_ns()
        result = fn()
        durations.append(time.perf_counter_ns() - start)
    return int(statistics.median(durations)), result


def planar(g):
    return nx.check_planarity(g)[0]


def iso(a, b, labels):
    # VF2++ 3.2.1 returns False when either graph has no vertices.
    if not (len(a) and len(b)):
        return nx.is_isomorphic(a, b)
    return nx.vf2pp_is_isomorphic(a, b, node_label="label" if labels else None)


def induced(host, pattern, labels):
    newer = getattr(nx, "vf2pp_subgraph_is_isomorphic", None)
    if newer is not None:
        return newer(host, pattern, node_label="label" if labels else None)
    if host.is_directed():
        factory = nx.isomorphism.MultiDiGraphMatcher if host.is_multigraph() else nx.isomorphism.DiGraphMatcher
    else:
        factory = nx.isomorphism.MultiGraphMatcher if host.is_multigraph() else nx.isomorphism.GraphMatcher
    match = nx.isomorphism.categorical_node_match("label", None) if labels else None
    return factory(host, pattern, node_match=match).subgraph_is_isomorphic()


def edge_subdivision(witness):
    """Suppress degree-two vertices and check the exact kernel independently."""
    if witness.is_multigraph() or nx.number_of_selfloops(witness) or not witness:
        return False
    degrees = dict(witness.degree())
    branches = [u for u, degree in degrees.items() if degree != 2]
    seen = set()
    reduced = nx.Graph()
    reduced.add_nodes_from(branches)
    for root in branches:
        for target in witness[root]:
            edge = frozenset((root, target))
            if edge in seen:
                continue
            seen.add(edge)
            previous = root
            while degrees[target] == 2:
                following = next(v for v in witness[target] if v != previous)
                edge = frozenset((target, following))
                if edge in seen:
                    return False
                seen.add(edge)
                previous, target = target, following
            if root == target or reduced.has_edge(root, target):
                return False
            reduced.add_edge(root, target)
    if len(seen) != witness.number_of_edges():
        return False
    return nx.is_isomorphic(reduced, nx.complete_graph(5)) or nx.is_isomorphic(
        reduced, nx.complete_bipartite_graph(3, 3)
    )


def modularity(g, communities, resolution):
    if g.size(weight="weight") == 0:
        return 0.0
    return nx.community.modularity(g, communities, resolution=resolution)


def atlas():
    print(nx.__version__)
    for g in nx.graph_atlas_g():
        print(len(g), int(planar(g)), *(v for edge in g.edges for v in edge))


def main():
    if sys.argv[1] == "atlas":
        atlas()
        return
    repeats = int(sys.argv[1])
    tokens = iter(sys.stdin.buffer.read().split())
    print(nx.__version__)
    for _ in range(int(next(tokens))):
        mode = next(tokens)
        if mode == b"P":
            g = graph(tokens)
            elapsed, result = measure(lambda: planar(g), repeats)
            print(elapsed, int(result))
        elif mode == b"I":
            a, b = graph(tokens), graph(tokens)
            labels = bool(int(next(tokens)))
            if labels:
                for g in (a, b):
                    for u in g:
                        g.nodes[u]["label"] = int(next(tokens))
            elapsed, result = measure(lambda: iso(a, b, labels), repeats)
            print(elapsed, int(result))
        elif mode == b"S":
            host, pattern = graph(tokens), graph(tokens)
            if host.is_multigraph() or pattern.is_multigraph():
                graph_type = nx.MultiDiGraph if host.is_directed() else nx.MultiGraph
                host, pattern = graph_type(host), graph_type(pattern)
            labels = bool(int(next(tokens)))
            if labels:
                for g in (host, pattern):
                    for u in g:
                        g.nodes[u]["label"] = int(next(tokens))
            elapsed, result = measure(lambda: induced(host, pattern, labels), repeats)
            print(elapsed, int(result))
        elif mode == b"C":
            g = nx.Graph(graph(tokens))
            g.remove_edges_from(list(nx.selfloop_edges(g)))
            certificate = nx.Graph()
            for _ in range(int(next(tokens))):
                u, v = int(next(tokens)), int(next(tokens))
                assert g.has_edge(u, v)
                certificate.add_edge(u, v)
            if certificate:
                assert edge_subdivision(certificate)
            elapsed, (is_planar_graph, nx_certificate) = measure(
                lambda: nx.check_planarity(g, counterexample=True), repeats
            )
            assert bool(certificate) != is_planar_graph
            if not is_planar_graph:
                assert edge_subdivision(nx_certificate)
            print(elapsed, int(is_planar_graph), nx_certificate.number_of_edges() if not is_planar_graph else 0)
        elif mode == b"L":
            g = graph(tokens)
            resolution, threshold = float(next(tokens)), float(next(tokens))
            seed, limit = int(next(tokens)), int(next(tokens))
            assignment = {}
            for u in g:
                group = int(next(tokens))
                assignment.setdefault(group, set()).add(u)
            if g.size(weight="weight") == 0:
                print(0, 0.0, len(g), 0.0)
                continue
            options = dict(resolution=resolution, threshold=threshold, seed=seed)
            def solve():
                if limit:
                    return list(islice(nx.community.louvain_partitions(g, **options), limit))[-1]
                return nx.community.louvain_communities(g, **options)
            elapsed, result = measure(solve, repeats)
            print(
                elapsed, modularity(g, result, resolution), len(result),
                modularity(g, list(assignment.values()), resolution),
            )
        elif mode == b"G":
            g = graph(tokens)
            name = next(tokens).decode()
            f = getattr(nx, name)
            params = [int(next(tokens)) for _ in range(int(next(tokens)))]
            if name.endswith("_lattice_graph"):
                params[2] = bool(params[2])
            print(0, int(nx.is_isomorphic(g, f(*params))))
        else:
            raise ValueError(mode)


if __name__ == "__main__":
    main()
