from collections import defaultdict
from itertools import permutations, groupby
from rosalind.helpers import recursionlimit
from os.path import commonprefix
from math import log
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


def collapse_trie(graph):
    """Collapses single descendent nodes in a trie to yield a suffix tree"""
    new = {}
    for k in graph.keys():
        edge = k
        g = graph[k]
        while len(g) == 1:
            k, g = list(g.items())[0]
            edge += k
        new[edge] = collapse_trie(g)
    return new


def get_edges(graph):
    for k in graph.keys():
        yield k
        yield from get_edges(graph[k])


def suff(seq):
    seqs = [seq[i:] for i in range(len(seq))]
    return collapse_trie(trie(seqs))


# This is naive solution, but takes too long to run
def ling_old(seq):
    tree = suff(seq)
    s = sum([len(edge) for edge in get_edges(tree)])
    m = sum(min(4 ** k, len(seq) - k + 1) for k in range(1, len(seq) + 1))
    return s / m


# This is an "optimised" solution to ling that doesn't compute a suffix tree
# Instead it sums the lengths of the edges as it goes.
# Memoisation leads to a large speed up and takes ~25 secs to run.
@cache
def compute_sub(seq, starts):
    if len(starts) == 1:
        return len(seq) - starts[0]
    else:
        tot = 0
        bases = set([seq[start] for start in starts])
        for base in bases:
            matching = [start for start in starts if seq[start] == base]
            seqs = [seq[s:] for s in matching]
            prefix = commonprefix(seqs)
            size = len(prefix)
            nstarts = [start + size for start in matching if start + size < len(seq)]
            tot += compute_sub(seq, tuple(nstarts)) + size
        return tot


def ling(seq):
    n = len(seq)
    m = sum([n - k + 1 if k > log(n + 1) / log(4) else 4 ** k for k in range(1, n + 1)])
    with recursionlimit(10000):
        sub = compute_sub(seq, tuple(range(len(seq))))
    return sub / m
