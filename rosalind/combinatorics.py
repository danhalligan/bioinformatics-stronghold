import re

from rosalind.helpers import codons, memoize
from math import prod, factorial
from functools import reduce
from itertools import permutations, product, combinations


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
    cod = codons()
    seq = seq + "*"
    return reduce(lambda p, c: p * cod[c] % mod, seq, 1)


def orf(seq):
    """Open Reading Frames"""
    for x in [seq, seq.revc()]:
        for i in range(3):
            subseq = x[i : len(x) - (len(x) - i) % 3]
            for m in re.finditer("(?=(M[^\\*]*)\\*)", subseq.translate().seq):
                yield m.group(1)


def pmch(seq):
    """Perfect Matchings and RNA Secondary Structures"""
    return prod([factorial(seq.count(x)) for x in "AG"])


def pper(n, k):
    """Partial Permutations"""
    return reduce(lambda p, i: int((p * i) % 1e6), range(n, n - k, -1))


def sign(n):
    """Enumerating Oriented Gene Orderings"""
    perm = list(permutations(range(1, n + 1)))
    sign = list(product([-1, 1], repeat=n))
    res = []
    for p in perm:
        for s in sign:
            res.append([x * y for x, y in zip(s, p)])
    return res


def valid_pair(x, y):
    pair = {"A": "U", "U": "A", "C": "G", "G": "C"}
    return x == pair[y]


@memoize
def cat(seq, mod=1000000):
    """Calculate total number of noncrossing perfect matchings"""
    if len(seq) in range(1):
        return 1
    else:
        tot = sum(
            cat(seq[1:m]) * cat(seq[m + 1 :])
            for m in range(1, len(seq), 2)
            if valid_pair(seq[0], seq[m])
        )
        return tot % mod


def mmch(seq):
    au = [seq.count(x) for x in "AU"]
    gc = [seq.count(x) for x in "GC"]
    return prod(au) * prod(gc)


def flip(x, i, j):
    """Flip a section of a sequence"""
    rev = list.copy(x)
    rev[i:j] = rev[i:j][::-1]
    return rev


def breaks(s, t):
    """Identify breaks between a sequence and target"""
    return [
        i + 1 for i in range(len(s) - 1) if abs(t.index(s[i]) - t.index(s[i + 1])) != 1
    ]


def min_breaks(seqs, t):
    rev = []
    for s in seqs:
        for i, j in combinations(breaks(s, t), 2):
            rev.append(flip(s, i, j))

    min_b = len(t)
    mr = []
    for r in rev:
        n = len(breaks(r, t))
        if n < min_b:
            min_b = n
            mr = [r]
        elif n == min_b:
            mr.append(r)
    return mr


# based on https://medium.com/@matthewwestmk/87c62d690eef
def get_distance(s, t):
    s = ["-"] + s + ["+"]
    t = ["-"] + t + ["+"]
    nr = 0
    c = [s]
    while t not in c:
        c = min_breaks(c, t)
        nr += 1
    return nr


@memoize
def motz(seq):
    if len(seq) in range(1):
        return 1
    else:
        tot = motz(seq[1:])
        tot += sum(
            motz(seq[1:m]) * motz(seq[m + 1 :])
            for m in range(1, len(seq))
            if valid_pair(seq[0], seq[m])
        )
        return tot % 10 ** 6
