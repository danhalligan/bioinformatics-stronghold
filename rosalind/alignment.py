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


def edit(s1, s2):
    """Edit Distance"""
    # initialise
    m = {}
    for j in range(len(s2) + 1):
        m[j, 0] = j

    for i in range(len(s1) + 1):
        m[0, i] = i

    # fill matrices
    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i]
            else:
                m[j + 1, i + 1] = min([m[j + 1, i], m[j, i], m[j, i + 1]]) + 1

    return m[len(s2), len(s1)]


def edta(s1, s2):
    """Edit Distance Alignment"""
    # initialise
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = j
        p[j, 0] = [j - 1, 0]

    for i in range(len(s1) + 1):
        m[0, i] = i
        p[0, i] = [0, i - 1]

    # fill matrices
    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i]
                p[j + 1, i + 1] = [j, i]
            else:
                opt = [m[j + 1, i], m[j, i], m[j, i + 1]]
                m[j + 1, i + 1] = min(opt) + 1
                p[j + 1, i + 1] = [[j + 1, i], [j, i], [j, i + 1]][opt.index(min(opt))]

    # traceback
    a1, a2 = "", ""
    i, j = len(s1), len(s2)
    while i > 0 and j > 0:
        if p[j, i] == [j - 1, i - 1]:
            a1 += s1[i - 1]
            a2 += s2[j - 1]
        elif p[j, i] == [j, i - 1]:
            a1 += s1[i - 1]
            a2 += "-"
        elif p[j, i] == [j - 1, i]:
            a1 += "-"
            a2 += s2[j - 1]
        j, i = p[j, i]

    return {"dist": m[len(s2), len(s1)], "a1": a1[::-1], "a2": a2[::-1]}


def scsp(s1, s2):
    """Interleaving Two Motifs"""

    # We can calculate a shortest common supersequence by constructing a
    # 2D matrix as before, but disallowing mismatches and with a special
    # traceback.

    # Initialise
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = j
        p[j, 0] = [j - 1, 0]

    for i in range(len(s1) + 1):
        m[0, i] = i
        p[0, i] = [0, i - 1]

    # fill matrices
    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i]
                p[j + 1, i + 1] = [j, i]
            else:
                opt = [m[j + 1, i], m[j, i + 1]]
                m[j + 1, i + 1] = min(opt) + 1
                p[j + 1, i + 1] = [[j + 1, i], [j, i + 1]][opt.index(min(opt))]

    # traceback
    ss = ""
    i, j = len(s1), len(s2)
    while i > 0 or j > 0:
        if p[j, i] == [j - 1, i - 1]:
            ss += s1[i - 1]
        elif p[j, i] == [j, i - 1]:
            ss += s1[i - 1]
        elif p[j, i] == [j - 1, i]:
            ss += s2[j - 1]
        j, i = p[j, i]

    return {"dist": m[len(s2), len(s1)], "ss": ss[::-1]}
