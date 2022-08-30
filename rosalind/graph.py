from collections import defaultdict
from itertools import permutations, groupby


def overlap_graph(seqs, n=3):
    """Build an overlap graph of sequences where prefix of A matches suffix of
    B"""
    for pair in permutations(seqs, 2):
        if pair[0].seq.endswith(pair[1].seq[:n]):
            yield (pair[0].id, pair[1].id)


def trie(seqs):
    graph = {}
    if len(seqs):
        for base, nseqs in groupby(seqs, key=lambda s: s[0]):
            graph[base] = trie([seq[1:] for seq in nseqs if len(seq) > 1])
    return graph


def print_trie(graph, node=1):
    i = node
    for edge in sorted(graph):
        print(node, i + 1, edge)
        i = print_trie(graph[edge], i + 1)
    return i


def build_seq(node, rev, edges):
    seq = ""
    while node in rev:
        # print(edges[node], end = " ")
        seq = edges[node] + seq
        node = rev[node]
    # print()
    return seq


# traverse graph and find longest seq up to a node with a least k leaves?
def lrep(seq, k, graph):
    # parse data and build structures
    edges = defaultdict(list)
    rev = {}
    heads, tails = set(), set()
    for edge in graph:
        n1, n2, i, n = edge.split()
        tails.add(n2)
        heads.add(n1)
        rev[n2] = n1
        edges[n2] = seq[(int(i) - 1) : (int(i) + int(n) - 1)]

    # count the number of descendents (leaves) from each node
    descendents = defaultdict(int)
    for leaf in tails - heads:
        while leaf in rev:
            leaf = rev[leaf]
            descendents[leaf] += 1

    candidates = [x for x in descendents if descendents[x] >= k]
    seqs = [build_seq(cand, rev, edges) for cand in candidates]
    return max(seqs, key=len)
