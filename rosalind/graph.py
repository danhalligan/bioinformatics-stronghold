from collections import defaultdict
from itertools import permutations


class Graph:
    """Store a graph in a dictionary where key is node and value is a list of
    connections"""

    def __init__(self, adjacency_list):
        self.graph = defaultdict(list)
        for x in adjacency_list:
            self.graph[x[0]].append(x[1])
            self.graph[x[1]].append(x[0])
        self.nodes = list(self.graph.keys())

    def count_distinct(self):
        visited = {}
        n_graphs = 0

        def visit_nodes(node):
            visited[node] = True
            for i in self.graph[node]:
                if i not in visited:
                    visit_nodes(i)

        for node in list(self.nodes):
            if node not in visited:
                visit_nodes(node)
                n_graphs += 1

        return n_graphs


def overlap_graph(seqs, n=3):
    """Build an overlap graph of sequences where prefix of A matches suffix of
    B"""
    for pair in permutations(seqs, 2):
        if pair[0].seq.endswith(pair[1].seq[:n]):
            yield (pair[0].id, pair[1].id)


def trie(v):
    def trier(v, r):
        if len(v) == 1 and len(v[0]) == 1:
            out.append([r, r + 1, v[0]])
            return r + 2
        elif len(v) == 1 and len(v[0]) == 0:
            return r + 1
        else:
            s = set([x[0] for x in v])
            n = r + 1
            for c in s:
                out.append([r, n, c])
                n = trier([x[1:] for x in v if x[0] == c], n)
            return n

    out = []
    trier(v, 1)
    return out
