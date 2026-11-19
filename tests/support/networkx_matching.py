"""NetworkX objective and timing oracle. Read integer test cases from stdin."""

import statistics
import sys
import time

import networkx as nx


def solve(graph, count, operation, repetitions):
    elapsed = []
    answer = None
    for _ in range(repetitions):
        start = time.perf_counter_ns()
        if operation == "maximum":
            matching = nx.max_weight_matching(graph)
        elif operation == "cardinality":
            matching = nx.max_weight_matching(graph, maxcardinality=True)
        else:
            matching = nx.min_weight_matching(graph)
        elapsed.append(time.perf_counter_ns() - start)
        assert nx.is_matching(graph, matching)
        result = (len(matching), sum(graph[a][b]["weight"] for a, b in matching))
        assert answer is None or result == answer
        answer = result
    assert all(a < count and b < count for a, b in matching)
    return (*answer, int(statistics.median(elapsed)))


def main():
    repetitions = int(sys.argv[1])
    terms = iter(map(int, sys.stdin.buffer.read().split()))
    print(nx.__version__)
    for _ in range(next(terms)):
        vertices, edges = next(terms), next(terms)
        high, low = nx.Graph(), nx.Graph()
        high.add_nodes_from(range(vertices))
        low.add_nodes_from(range(vertices))
        for _ in range(edges):
            a, b, weight = next(terms), next(terms), next(terms)
            if a == b:
                continue
            if not high.has_edge(a, b) or weight > high[a][b]["weight"]:
                high.add_edge(a, b, weight=weight)
            if not low.has_edge(a, b) or weight < low[a][b]["weight"]:
                low.add_edge(a, b, weight=weight)
        for graph, operation in [(high, "maximum"), (high, "cardinality"), (low, "minimum")]:
            if repetitions > 1:
                solve(graph, vertices, operation, 1)
            print(*solve(graph, vertices, operation, repetitions))


if __name__ == "__main__":
    main()
