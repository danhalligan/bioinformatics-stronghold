import re

from math import prod, factorial
from functools import reduce


def fib(n, k):
    """Rabbits and Recurrence Relations"""
    a, b = 1, 1
    for i in range(2, n):
        a, b = b, k * a + b
    return b


def fibd(n, m):
    """Mortal Fibonacci Rabbits"""
    v = [1] + (m - 1) * [0]
    for i in range(2, n + 1):
        v = [sum(v[1:])] + v[:-1]
    return sum(v)


def mrna(seq, mod=1000000):
    """Inferring mRNA from Protein"""
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


def orf(seq):
    """Open Reading Frames"""
    for x in [seq, seq.revc()]:
        x = str(x.rna())
        for i in range(3):
            subseq = x[i : len(x) - (len(x) - i) % 3]
            for m in re.finditer("(?=(M[^\\*]*)\\*)", translate(subseq)):
                yield m.group(1)


def pmch(seq):
    """Perfect Matchings and RNA Secondary Structures"""
    return prod([factorial(seq.count(x)) for x in "AG"])


def pper(n, k):
    """Partial Permutations"""
    return reduce(lambda p, i: int((p * i) % 1e6), range(n, n - k, -1))
