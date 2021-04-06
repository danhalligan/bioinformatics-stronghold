import re
import numpy as np
import requests as r

from rosalind.helpers import read_fasta, Dna
from io import StringIO
from itertools import product


def max_gc(seqs):
    gc = [Dna(rec.seq).gc_content() for rec in seqs]
    m = gc.index(max(gc))
    return {"name": seqs[m].id, "value": round(gc[m] * 100, 5)}


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


def lcsm(seqs):
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


def uniprot_output(id):
    return r.post(f"https://www.uniprot.org/uniprot/{id}.fasta").text


def get_uniprot(id):
    txt = uniprot_output(id)
    return list(read_fasta(StringIO(txt)))[0]


def find_protein_motif(seq, pattern="N[^P][ST][^P]"):
    motif = re.compile("(?=(" + pattern + "))")
    return [m.start() + 1 for m in motif.finditer(seq)]


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


def lgis(x):
    """DP approach to longest increasing subsequence"""
    n = len(x)
    d = [1] * n  # length of the longest increasing subsequence ending at i
    p = [-1] * n  # pointer to previous part of subsequence

    for i in range(n):
        for j in range(i):
            if x[j] < x[i] and d[i] < d[j] + 1:
                d[i] = d[j] + 1
                p[i] = j

    ans = max(d)
    pos = d.index(ans)

    # traceback
    subseq = []
    while pos != -1:
        subseq.append(x[pos])
        pos = p[pos]

    subseq.reverse()
    return subseq


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


def lexv(s, n):
    s = ["_"] + s
    perm = list(product(s, repeat=n))
    perm = ["".join(x) for x in perm]
    perm = [re.sub("_+$", "", x) for x in perm]
    return list(filter(lambda x: "_" not in x, perm[1:]))
