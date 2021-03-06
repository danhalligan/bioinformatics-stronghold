import re
import numpy as np
from rosalind.helpers import *


def count_nucleotides(seq):
    """Count number of 'A', 'C', 'G', and 'T' occurrences"""
    return [seq.count(nuc) for nuc in "ACGT"]


def dna2rna(seq):
    """Convert DNA to RNA"""
    return re.sub("T", "U", seq)


def revcomp(seq):
    """Reverse complement"""
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def gc_content(seq):
    """Calculate GC content"""
    return sum([seq.count(base) for base in "GC"]) / len(seq)


def max_gc(seqs):
    gc = [gc_content(rec.seq) for rec in seqs]
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
    from math import comb

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
    from math import comb

    def dbinom(x, size, prob):
        return comb(size, x) * prob ** x * (1 - prob) ** (size - x)

    def pbinom(q, size, prob):
        return sum([dbinom(x, size, prob) for x in range(0, q + 1)])

    return 1 - pbinom(n - 1, 2 ** k, 0.25)


def get_uniprot(id):
    import requests as r
    from io import StringIO

    response = r.post(f"https://www.uniprot.org/uniprot/{id}.fasta")
    return list(read_fasta(StringIO(response.text)))[0]


def find_protein_motif(seq, pattern="N[^P][ST][^P]"):
    motif = re.compile("(?=(" + pattern + "))")
    return [m.start() + 1 for m in motif.finditer(seq)]


def count_rnas(seq, mod=1000000):
    from functools import reduce

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
    for s in [seq, revcomp(seq)]:
        s = dna2rna(s)
        for i in range(3):
            l = len(s)
            subseq = s[i : l - (l - i) % 3]
            for m in re.finditer("(?=(M[^\\*]*)\\*)", translate(subseq)):
                yield m.group(1)


def protein_mass(prot):
    mass = {
        "A": 71.03711,
        "C": 103.00919,
        "D": 115.02694,
        "E": 129.04259,
        "F": 147.06841,
        "G": 57.02146,
        "H": 137.05891,
        "I": 113.08406,
        "K": 128.09496,
        "L": 113.08406,
        "M": 131.04049,
        "N": 114.04293,
        "P": 97.05276,
        "Q": 128.05858,
        "R": 156.10111,
        "S": 87.03203,
        "T": 101.04768,
        "V": 99.06841,
        "W": 186.07931,
        "Y": 163.06333,
    }

    return sum([mass[x] for x in prot])


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
