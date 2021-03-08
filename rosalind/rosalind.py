import re
import numpy as np
import requests as r

from rosalind.helpers import read_fasta, Dna
from itertools import permutations
from math import comb
from functools import reduce
from io import StringIO
from collections import defaultdict


def max_gc(seqs):
    gc = [Dna(rec.seq).gc_content() for rec in seqs]
    m = gc.index(max(gc))
    return {"name": seqs[m].id, "value": round(gc[m] * 100, 5)}


def hamming(s1, s2):
    """Calculate Hamming distance"""
    return sum(xi != yi for xi, yi in zip(s1, s2))


def rabbits(n, k):
    """Divide and conquer solution to Fibonacci sequence"""
    a, b = 1, 1
    for i in range(2, n):
        a, b = b, k * a + b
    return b


def mendel1(k, m, n):
    tot = comb(k + m + n, 2)
    poss = comb(k, 2) + k * m + k * n + m * n / 2 + comb(m, 2) * 3 / 4
    return poss / tot


def genetic_code():
    return {
        "UUU": "F",
        "UCU": "S",
        "UAU": "Y",
        "UGU": "C",
        "UUC": "F",
        "UCC": "S",
        "UAC": "Y",
        "UGC": "C",
        "UUA": "L",
        "UCA": "S",
        "UAA": "*",
        "UGA": "*",
        "UUG": "L",
        "UCG": "S",
        "UAG": "*",
        "UGG": "W",
        "CUU": "L",
        "CCU": "P",
        "CAU": "H",
        "CGU": "R",
        "CUC": "L",
        "CCC": "P",
        "CAC": "H",
        "CGC": "R",
        "CUA": "L",
        "CCA": "P",
        "CAA": "Q",
        "CGA": "R",
        "CUG": "L",
        "CCG": "P",
        "CAG": "Q",
        "CGG": "R",
        "AUU": "I",
        "ACU": "T",
        "AAU": "N",
        "AGU": "S",
        "AUC": "I",
        "ACC": "T",
        "AAC": "N",
        "AGC": "S",
        "AUA": "I",
        "ACA": "T",
        "AAA": "K",
        "AGA": "R",
        "AUG": "M",
        "ACG": "T",
        "AAG": "K",
        "AGG": "R",
        "GUU": "V",
        "GCU": "A",
        "GAU": "D",
        "GGU": "G",
        "GUC": "V",
        "GCC": "A",
        "GAC": "D",
        "GGC": "G",
        "GUA": "V",
        "GCA": "A",
        "GAA": "E",
        "GGA": "G",
        "GUG": "V",
        "GCG": "A",
        "GAG": "E",
        "GGG": "G",
    }


def translate(seq):
    code = genetic_code()
    return "".join([code[seq[i : i + 3]] for i in range(0, len(seq), 3)])[:-1]


def find_motif(seq1, seq2):
    size = len(seq2)
    for i in range(len(seq1) - size + 1):
        if seq1[i : (i + size)] == seq2:
            yield i + 1


def profile_matrix(seqs):
    def count_bases(v):
        return [sum(v == b) for b in "ACGT"]

    x = np.array([list(x) for x in seqs])
    return np.array([count_bases(v) for v in x.T]).T


def consensus_sequence(mat):
    return "".join(["ACGT"[np.argmax(v)] for v in mat.T])


def mortal_rabbits(n, m):
    v = [1] + (m - 1) * [0]
    for i in range(2, n + 1):
        v = [sum(v[1:])] + v[:-1]
    return sum(v)


def expected_offspring(v):
    p = [1, 1, 1, 0.75, 0.5, 0]
    return sum([x[0] * x[1] * 2 for x in zip(v, p)])


def find_shared_motif(seqs):
    seqs = sorted(seqs, key=len)

    maxsub = ""
    s0 = seqs.pop()
    s = len(s0)

    for i in range(s):
        for j in range(i + len(maxsub), s):
            for seq in seqs:
                if s0[i:j] not in seq:
                    break
                else:
                    maxsub = s0[i:j]

    return maxsub


def mendel2(k, n):
    def dbinom(x, size, prob):
        return comb(size, x) * prob ** x * (1 - prob) ** (size - x)

    def pbinom(q, size, prob):
        return sum([dbinom(x, size, prob) for x in range(0, q + 1)])

    return 1 - pbinom(n - 1, 2 ** k, 0.25)


def get_uniprot(id):
    response = r.post(f"https://www.uniprot.org/uniprot/{id}.fasta")
    return list(read_fasta(StringIO(response.text)))[0]


def find_protein_motif(seq, pattern="N[^P][ST][^P]"):
    motif = re.compile("(?=(" + pattern + "))")
    return [m.start() + 1 for m in motif.finditer(seq)]


def count_rnas(seq, mod=1000000):

    codons = {
        "F": 2,
        "L": 6,
        "I": 3,
        "M": 1,
        "V": 4,
        "S": 6,
        "P": 4,
        "T": 4,
        "A": 4,
        "Y": 2,
        "H": 2,
        "Q": 2,
        "N": 2,
        "K": 2,
        "D": 2,
        "E": 2,
        "C": 2,
        "W": 1,
        "R": 6,
        "G": 4,
        "*": 3,
    }
    seq = seq + "*"
    return reduce(lambda p, c: p * codons[c] % mod, seq, 1)


def find_orfs(seq):
    for s in [seq, seq.revc()]:
        s = str(s.rna())
        for i in range(3):
            subseq = s[i : len(s) - (len(s) - i) % 3]
            for m in re.finditer("(?=(M[^\\*]*)\\*)", translate(subseq)):
                yield m.group(1)


def reverse_pallindromes(seq):
    comp = seq.translate(str.maketrans("ACGT", "TGCA"))
    n = len(seq)
    for i in range(n):
        for size in range(2, 7):
            lo = i - size
            hi = i + size
            if lo < 0 or hi > n or seq[lo:hi] != comp[lo:hi][::-1]:
                break
            else:
                yield [i - size + 1, size * 2]


def kmp_preprocess(seq):
    """KMP preprocessing algorithm"""
    j = -1
    b = [j]
    for i in range(len(seq)):
        while j >= 0 and seq[i] != seq[j]:
            j = b[j]
        j += 1
        b.append(j)
    return b[1:]


def overlap_graph(seqs, n=3):
    """Build an overlap graph of sequences where prefix of A matches suffix of
    B"""
    for pair in permutations(seqs, 2):
        if pair[0].seq.endswith(pair[1].seq[:n]):
            yield (pair[0].id, pair[1].id)


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
