def hamm(s1, s2):
    """Calculate Hamming distance"""
    return sum(xi != yi for xi, yi in zip(s1, s2))


def ts(x, y):
    """Transition?"""
    ts = {"A": "G", "C": "T", "G": "A", "T": "C"}
    return x == ts[y]


def tran(seqs):
    """Transitions and Transversions"""
    mm, tr = 0, 0
    for x, y in zip(seqs[0], seqs[1]):
        if x != y:
            mm += 1
            tr += int(ts(x, y))

    return tr / (mm - tr)


def pdst(seqs):
    """Creating a Distance Matrix"""
    n = len(seqs)
    d = [[0.0 for x in range(n)] for y in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d[i][j] = d[j][i] = hamm(seqs[i], seqs[j]) / len(seqs[i])
    return d
