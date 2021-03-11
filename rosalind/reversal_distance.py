from itertools import combinations


# based on https://medium.com/@matthewwestmk/87c62d690eef
def flip(x, i, j):
    rev = list.copy(x)
    rev[i:j] = rev[i:j][::-1]
    return rev


def breaks(s, t):
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


def get_distance(s, t):
    s = ["-"] + s + ["+"]
    t = ["-"] + t + ["+"]
    nr = 0
    c = [s]
    while t not in c:
        c = min_breaks(c, t)
        nr += 1
    return nr
