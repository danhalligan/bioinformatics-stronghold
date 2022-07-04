from rosalind.helpers import blosum62


def printm(m, s1, s2):
    """print alignment dict (for debugging)"""
    print(" ", end="")
    s1 = " " + s1
    for i in range(len(s1)):
        print("%(number)4s" % {"number": s1[i]}, end="")
    print()
    s2 = " " + s2
    for j in range(len(s2)):
        print(s2[j], end="")
        for i in range(len(s1)):
            print("%(number)4s" % {"number": m[j, i]}, end="")
        print()
    print()


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


def lcsq(s1, s2):
    """Finding a Shared Spliced Motif"""
    # initialise
    m, p = {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = 0
        p[j, 0] = [j - 1, 0]

    for i in range(len(s1) + 1):
        m[0, i] = 0
        p[0, i] = [0, i - 1]

    # fill matrices
    for j in range(len(s2)):
        for i in range(len(s1)):
            if s1[i] == s2[j]:
                m[j + 1, i + 1] = m[j, i] + 1
                p[j + 1, i + 1] = [j, i]
            else:
                opt = [m[j + 1, i], m[j, i + 1]]
                m[j + 1, i + 1] = max(opt)
                p[j + 1, i + 1] = [[j + 1, i], [j, i + 1]][opt.index(max(opt))]

    # traceback
    subs = ""
    i, j = len(s1), len(s2)
    while i > 0 and j > 0:
        if p[j, i] == [j - 1, i - 1]:
            subs += s1[i - 1]
        j, i = p[j, i]

    return subs[::-1]


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


def ctea(s1, s2):
    """Counting Optimal Alignments"""

    score, routes = {}, {}
    for j in range(len(s2) + 1):
        score[j, 0] = j
        routes[j, 0] = 1

    for i in range(len(s1) + 1):
        score[0, i] = i
        routes[0, i] = 1

    for j in range(len(s2)):
        for i in range(len(s1)):
            pos = [(j + 1, i), (j, i), (j, i + 1)]
            cost = [1, int(s1[i] != s2[j]), 1]
            scores = [score[pos[x]] + cost[x] for x in range(3)]
            best = min(scores)
            new = (j + 1, i + 1)
            score[new] = best
            routes[new] = sum(routes[pos[x]] for x in range(3) if scores[x] == best)

    return routes[len(s2), len(s1)] % 134217727


def glob(s1, s2, penalty=-5):
    """Global Alignment with Scoring Matrix"""
    score = blosum62()
    m = {}
    for j in range(len(s2) + 1):
        m[j, 0] = penalty * j
    for i in range(len(s1) + 1):
        m[0, i] = penalty * i

    for j in range(len(s2)):
        for i in range(len(s1)):
            pos = [(j + 1, i), (j, i), (j, i + 1)]
            cost = [penalty, score[s1[i]][s2[j]], penalty]
            scores = [m[pos[x]] + cost[x] for x in range(3)]
            m[j + 1, i + 1] = max(scores)

    return m[len(s2), len(s1)]


def gcon(s1, s2, penalty=-5):
    """Global Alignment with Constant Gap Penalty"""
    # See Biological Sequence Analysis page 29
    # Its a simplification of the more general affine gap penalty
    # with an extension penalty of 0.
    # We now have to keep track of three matrices m, x and y

    score = blosum62()
    m, x, y = {}, {}, {}
    for j in range(len(s2) + 1):
        m[j, 0] = penalty
        y[j, 0] = -99
    for i in range(len(s1) + 1):
        m[0, i] = penalty
        x[0, i] = -99

    m[0, 0] = 0
    for j in range(len(s2)):
        for i in range(len(s1)):
            new = (j + 1, i + 1)
            x[new] = max([m[j, i + 1] + penalty, x[j, i + 1]])
            y[new] = max([m[j + 1, i] + penalty, y[j + 1, i]])
            m[new] = max([m[j, i] + score[s1[i]][s2[j]], x[new], y[new]])

    return m[len(s2), len(s1)]


# Quick function to insert indels.
def insert_indel(word, i):
    return word[:i] + "-" + word[i:]


def gaff(v, w, sigma=-11, epsilon=-1):
    """Global Alignment with Affine Gap Penalty"""
    # See Biological Sequence Analysis page 29
    # We now have to keep track of three matrices m, x and y

    score = blosum62()

    m = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    x = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    y = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    px = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    pm = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]
    py = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]

    # Initialize the edges with the given penalties.
    for i in range(1, len(v) + 1):
        x[i][0] = sigma + (i - 1) * epsilon
        m[i][0] = sigma + (i - 1) * epsilon
        y[i][0] = 10 * sigma
    for j in range(1, len(w) + 1):
        y[0][j] = sigma + (j - 1) * epsilon
        m[0][j] = sigma + (j - 1) * epsilon
        x[0][j] = 10 * sigma

    # Fill in the scores for the lower, middle, upper, and backtrack matrices.
    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            s = [x[i - 1][j] + epsilon, m[i - 1][j] + sigma]
            x[i][j] = max(s)
            px[i][j] = s.index(x[i][j])

            s = [y[i][j - 1] + epsilon, m[i][j - 1] + sigma]
            y[i][j] = max(s)
            py[i][j] = s.index(y[i][j])

            s = [x[i][j], m[i - 1][j - 1] + score[v[i - 1]][w[j - 1]], y[i][j]]
            m[i][j] = max(s)
            pm[i][j] = s.index(m[i][j])

    # Initialize the values of i, j and the aligned sequences.
    i, j = len(v), len(w)
    a1, a2 = v, w

    # Get the maximum score, and the corresponding backtrack starting position.
    scores = [x[i][j], m[i][j], y[i][j]]
    max_score = max(scores)
    s = scores.index(max_score)

    # Backtrack to the edge of the matrix starting bottom right.
    while i * j != 0:
        if s == 0:
            if px[i][j] == 1:
                s = 1
            i -= 1
            a2 = insert_indel(a2, j)
        elif s == 1:
            if pm[i][j] == 1:
                i -= 1
                j -= 1
            else:
                s = pm[i][j]
        else:
            if py[i][j] == 1:
                s = 1
            j -= 1
            a1 = insert_indel(a1, i)

    # Prepend the necessary preceeding indels to get to (0,0).
    for _ in range(i):
        a2 = insert_indel(a2, 0)
    for _ in range(j):
        a1 = insert_indel(a1, 0)

    return {"dist": str(max_score), "a1": a1, "a2": a2}
