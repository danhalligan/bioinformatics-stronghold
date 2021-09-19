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


# def gaff(s1, s2, penalty=-5, extension=0):
#     """Global Alignment with Affine Gap Penalty"""
#     # See Biological Sequence Analysis page 29
#     # We now have to keep track of three matrices m, x and y
#     score = blosum62()
#     m, x, y = {}, {}, {}
#     pm, px, py = {}, {}, {}

#     for j in range(len(s2) + 1):
#         m[j, 0] = penalty + extension * j
#         y[j, 0] = -99
#         pm[j, 0] = "↑"
#         px[j, 0] = "↑"
#         py[j, 0] = "↑"
#     for i in range(len(s1) + 1):
#         m[0, i] = penalty + extension * i
#         x[0, i] = -99
#         pm[0, i] = "←"
#         px[0, i] = "←"
#         py[0, i] = "←"

#     m[0, 0] = 0
#     for j in range(len(s2)):
#         for i in range(len(s1)):
#             new = (j + 1, i + 1)

#             scores = [m[j, i + 1] + penalty, x[j, i + 1] + extension]
#             x[new] = max(scores)
#             px[new] = ["↖", "←"][scores.index(x[new])]

#             scores = [m[j + 1, i] + penalty, y[j + 1, i] + extension]
#             y[new] = max(scores)
#             py[new] = ["↖", "↑"][scores.index(y[new])]

#             scores = [m[j, i] + score[s1[i]][s2[j]], x[new], y[new]]
#             m[new] = max(scores)
#             pm[new] = ["↖", "↑", "←"][scores.index(m[new])]

#     # traceback
#     i, j = len(s1), len(s2)
#     scores = [x[j, i], y[j, i], m[j, i]]
#     best = max(scores)
#     t = scores.index(best)
#     a1, a2 = "", ""
#     # p = [px, py, pm][scores.index(best)]

#     while i > 0 or j > 0:
#         print(i, j, t)
#         if t == 0:
#             if px[j, i] == "↖":
#                 t = 2
#             i -= 1
#             a1 += s1[i - 1]
#             a2 += "-"
#         elif t == 1:
#             if py[j, i] == "↖":
#                 t = 2
#             j -= 1
#             a1 += "-"
#             a2 += s2[j - 1]
#         elif t == 2:
#             if pm[j, i] == "←":
#                 t = 0
#             elif pm[j, i] == "↑":
#                 t = 1
#             else:
#                 a1 += s1[i - 1]
#                 a2 += s2[j - 1]
#                 i -= 1
#                 j -= 1

#     return {"dist": best, "a1": a1[::-1], "a2": a2[::-1]}


def gaff(s1, s2, penalty=-11, extension=-1):
    """Global Alignment with Affine Gap Penalty"""
    # See Biological Sequence Analysis page 29
    # We now have to keep track of three matrices m, x and y
    score = blosum62()
    m, x, y = {}, {}, {}
    pm, px, py = {}, {}, {}

    m[0, 0] = 0
    x[0, 0] = 0
    y[0, 0] = 0
    for i in range(1, len(s1) + 1):
        x[i, 0] = penalty + (i - 1) * extension
        m[i, 0] = penalty + (i - 1) * extension
        y[i, 0] = -99
    for j in range(1, len(s2) + 1):
        y[0, j] = penalty + (j - 1) * extension
        m[0, j] = penalty + (j - 1) * extension
        x[0, j] = -99

    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            scores = [m[i - 1, j] + penalty, x[i - 1, j] + extension]
            x[i, j] = max(scores)
            px[i, j] = ["↖", "←"][scores.index(x[i, j])]

            scores = [m[i, j - 1] + penalty, y[i, j - 1] + extension]
            y[i, j] = max(scores)
            py[i, j] = ["↖", "↑"][scores.index(y[i, j])]

            scores = [m[i - 1, j - 1] + score[s1[i - 1]][s2[j - 1]], x[i, j], y[i, j]]
            m[i, j] = max(scores)
            pm[i, j] = ["↖", "←", "↑"][scores.index(m[i, j])]

    # printm(m, s2, s1)
    # printm(x, s2, s1)
    # printm(y, s2, s1)
    i, j = len(s1), len(s2)
    scores = [m[i, j], x[i, j], y[i, j]]
    best = max(scores)
    p = ["↖", "←", "↑"][scores.index(best)]
    a1, a2 = "", ""
    pm[0, 0] = 0

    i, j = len(s1) - 1, len(s2) - 1
    while i >= 0 or j >= 0:
        if p == "←":
            p = px[i, j]
            i -= 1
            a1 += s1[i]
            a2 += "-"
        elif p == "↖":
            p = pm[i, j]
            a1 += s1[i]
            a2 += s2[j]
            i -= 1
            j -= 1
        else:
            p = py[i, j]
            j -= 1
            a1 += "-"
            a2 += s2[j]

    return {"dist": str(best), "a1": a1[::-1], "a2": a2[::-1]}


# def gaff(s1, s2, penalty = -11, extension = -1):
#     """Global Alignment with Affine Gap Penalty"""
#     # See Biological Sequence Analysis page 29
#     # We now have to keep track of three matrices m, x and y
#     score = blosum62()
#     m, x, y = {}, {}, {}
#     pm, px, py = {}, {}, {}

#     m[0, 0] = 0
#     x[0, 0] = 0
#     y[0, 0] = 0
#     for i in range(1, len(s1)+1):
#         x[i, 0] = penalty + (i-1)*extension
#         m[i, 0] = penalty + (i-1)*extension
#         y[i, 0] = -999
#     for j in range(1, len(s2)+1):
#         y[0, j] = penalty + (j-1)*extension
#         m[0, j] = penalty + (j-1)*extension
#         x[0, j] = -999

#     for i in range(1, len(s1)+1):
#         for j in range(1, len(s2)+1):
#             scores = [m[i-1, j] + penalty, x[i-1, j] + extension]
#             x[i, j] = max(scores)
#             px[i, j] = ['↖', '←'][scores.index(x[i, j])]

#             scores = [m[i, j-1] + penalty, y[i, j-1] + extension]
#             y[i, j] = max(scores)
#             py[i, j] = ['↖', '↑'][scores.index(y[i, j])]

#             scores = [m[i-1, j-1] + score[s1[i-1]][s2[j-1]], y[i, j], x[i, j]]
#             m[i, j] = max(scores)
#             pm[i, j] = ['↖', '↑', '←'][scores.index(m[i, j])]

#     i, j = len(s1), len(s2)
#     a1, a2 = s1, s2

#     matrix_scores = [x[i, j], m[i, j], y[i, j]]
#     max_score = max(matrix_scores)
#     backtrack_matrix = matrix_scores.index(max_score)

#     insert_indel = lambda word, i: word[:i] + '-' + word[i:]

#     while i > 0 or j > 0:
#         if backtrack_matrix == 0:
#             if px[i, j] == '↖':
#                 backtrack_matrix = 1
#             i -= 1
#             a2 = insert_indel(a2, j)
#         elif backtrack_matrix == 1:
#             if pm[i, j] == '←':
#                 backtrack_matrix = 0
#             elif pm[i, j] == '↑':
#                 backtrack_matrix = 2
#             else:
#                 i -= 1
#                 j -= 1
#         else:
#             if py[i, j] == '↖':
#                 backtrack_matrix = 1
#             j -= 1
#             a1 = insert_indel(a1, i)

#     return {"dist": str(max_score), "a1": a1, "a2": a2}


# Correct answers

# 348
# EYNGPMSRGRIRQYWYKWFDVRHIHDFAYYSNIYFAGA----KAKQPFDWLSYDHCADLGWDDSNPSYEMVRYYGNGLNM
# EYNGPMSRANIMQFWYKWFCEMHIHDFHYYSNIYFARADMNWKGAKAFDWNSYDHCADLGWDDSNPSYEMVRYYGNGLNM

# 326
# QVRPHADVNWMGVEPHMEYCFFTPWGIAYHRMSPASVQFLHAMRNFWASCGQSNIHEVMKMLACEHQCYLRHNLH--------CEM
# QVRPHADVNWMGVEPHMEYLFFTPWGIAYHVGA--------SKRNFWASCGQSNIHEVMKMLACEHQCYLRHELKVHWVQADGCEM

# 216
# DQQAEEY----NEMITRKFQ--QREMAHPTMNCTEHMMMFTTCNQKHVPGNQVGSIIPTRKSFDKCVNSRPHCYHER-DGL
# DQQAEQFIVQTEGMITGTFQFDHLEMAHPTMNCTTHMLSNLGLGQKVFYWLQVGSIIPTRKSFDKSVSSRPHCYHESYDGL

# 307
# FEAMNDVCPPINNCKHLFQPCEIYLPSPYDC--------DMKHHAHVYLKSWP-----------KHWPWAIYCVFDKRWLEHKFDDIMHHIHRCAFIERQCTSQLW
# -EAMNDVCPPI--------PC-IYLPSPYDTETDVVPCQDMKHHAHVYLKSWQAARCFMKYSMIKHWPWFIFCVF-----ENKFDDIMHHIHRCAFGERQCTSQLW
