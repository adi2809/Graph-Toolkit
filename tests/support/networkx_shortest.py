"""Distance, negative-cycle, and timing oracle for signed shortest paths."""

import statistics
import sys
import time

import networkx as nx


def bellman(graph, source):
    _, distances = nx.bellman_ford_predecessor_and_distance(graph, source)
    return distances


def johnson(graph, _source):
    return nx.johnson(graph)


def floyd(graph, _source):
    _, distances = nx.floyd_warshall_predecessor_and_distance(graph)
    return distances


def run(function, graph, source, repeats):
    samples = []
    result = None
    for _ in range(repeats):
        start = time.perf_counter_ns()
        try:
            result = function(graph, source)
        except nx.NetworkXUnbounded:
            result = None
        samples.append(time.perf_counter_ns() - start)
    return int(statistics.median(samples)), result


def encode(values):
    if values is None:
        return "cycle"
    return " ".join("x" if value is None or value == float("inf") else str(value) for value in values)


def main():
    repeats = int(sys.argv[1])
    tokens = iter(map(int, sys.stdin.buffer.read().split()))
    cases = next(tokens)
    print(nx.__version__)
    for _ in range(cases):
        n, m, source, directed = [next(tokens) for _ in range(4)]
        graph = nx.DiGraph() if directed else nx.Graph()
        graph.add_nodes_from(range(n))
        for _ in range(m):
            a, b, cost = [next(tokens) for _ in range(3)]
            if not graph.has_edge(a, b) or graph[a][b]["weight"] > cost:
                graph.add_edge(a, b, weight=cost)
        # The public predecessor API in some NetworkX releases checks
        # self-loops globally. Restrict it to vertices reachable from source.
        reachable = nx.descendants(graph, source) | {source}
        source_graph = graph if len(reachable) == n else graph.subgraph(reachable).copy()
        for function in (bellman, johnson, floyd):
            argument = source_graph if function is bellman else graph
            if repeats > 1:
                run(function, argument, source, 1)
            elapsed, result = run(function, argument, source, repeats)
            if result is None:
                values = None
            elif function is bellman:
                values = (result.get(node) for node in range(n))
            elif function is johnson:
                def distance(a, b):
                    path = result[a].get(b)
                    if path is None:
                        return None
                    return sum(graph[u][v]["weight"] for u, v in zip(path, path[1:]))
                values = (distance(a, b) for a in range(n) for b in range(n))
            elif any(result[node].get(node, 0) < 0 for node in range(n)):
                values = None
            else:
                values = (result[a].get(b) for a in range(n) for b in range(n))
            print(elapsed, encode(values))


if __name__ == "__main__":
    main()
