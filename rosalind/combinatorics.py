import re

from rosalind.helpers import genetic_code, codons
from math import prod, factorial
from functools import reduce
from itertools import permutations, product


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


def memoize(f):
    cache = {}

    def wrapper(*args):
        if args not in cache:
            cache[args] = f(*args)
        return cache[args]

    return wrapper


@memoize
def cat(seq):
    if len(seq) in range(1):
        return 1
    else:
        tot = sum(
            cat(seq[1:m]) * cat(seq[m + 1 :])
            for m in range(1, len(seq), 2)
            if valid_pair(seq[0], seq[m])
        )
        return tot % 10 ** 6


def mmch(seq):
    au = [seq.count(x) for x in "AU"]
    gc = [seq.count(x) for x in "GC"]
    return prod(au) * prod(gc)
