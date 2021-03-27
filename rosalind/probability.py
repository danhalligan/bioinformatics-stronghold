from math import comb, log10


def dbinom(x, size, prob):
    """Binomial density"""
    return comb(size, x) * prob ** x * (1 - prob) ** (size - x)


def pbinom(q, size, prob):
    """Binomial distribution function"""
    return sum([dbinom(x, size, prob) for x in range(0, q + 1)])


def iprb(k, m, n):
    tot = comb(k + m + n, 2)
    poss = comb(k, 2) + k * m + k * n + m * n / 2 + comb(m, 2) * 3 / 4
    return poss / tot


def expected_offspring(v):
    p = [1, 1, 1, 0.75, 0.5, 0]
    return sum([x[0] * x[1] * 2 for x in zip(v, p)])


def mendel2(k, n):
    return 1 - pbinom(n - 1, 2 ** k, 0.25)


def eval(n, s, x):
    gc = sum([s.count(x) for x in "GC"])
    at = len(s) - gc
    p = ((1 - x) / 2) ** (at) * (x / 2) ** gc
    return p * (n - len(s) + 1)


def indc(n):
    res = [1 - pbinom(x, 2 * n, 0.5) for x in range(2 * n)]
    return [round(log10(x), 3) for x in res]
