from collections import defaultdict
from itertools import permutations, groupby
from os.path import commonprefix
from functools import cache


def overlap_graph(seqs, n=3):
    for pair in permutations(seqs, 2):
        if pair[0].seq.endswith(pair[1].seq[:n]):
            yield (pair[0].id, pair[1].id)


def trie(seqs):
    graph = {}
    if len(seqs):
        for base, nseqs in groupby(sorted(seqs), key=lambda s: s[0]):
            graph[base] = trie([seq[1:] for seq in nseqs if len(seq) > 1])
    return graph


def as_adjacency(graph, nodes=[]):
    node = len(nodes) + 1
    for edge in sorted(graph):
        i = len(nodes) + 1
        nodes.append((node, i + 1, edge))
        nodes + as_adjacency(graph[edge], nodes)
    return nodes


def build_seq(node, rev, edges):
    seq = ""
    while node in rev:
        seq = edges[node] + seq
        node = rev[node]
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


@cache
def suffix_tree(seq, starts):
    graph = {}
    bases = sorted(set([seq[start] for start in starts]))
    for base in bases:
        matching = [start for start in starts if seq[start] == base]
        seqs = [seq[s:] for s in matching]
        prefix = commonprefix(seqs)
        size = len(prefix)
        new_starts = [start + size for start in matching if start + size < len(seq)]
        graph[prefix] = suffix_tree(seq, tuple(new_starts))
    return graph


def suff(seq):
    return suffix_tree(seq, tuple(range(len(seq))))


def get_edges(graph):
    for k in graph.keys():
        yield k
        yield from get_edges(graph[k])


def ling(seq):
    s = sum(len(edge) for edge in get_edges(suff(seq)))
    m = sum(min(4 ** k, len(seq) - k + 1) for k in range(1, len(seq) + 1))
    return s / m
