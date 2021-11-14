from math import floor
from rosalind.helpers import Dna
import rosalind.alignment as aln
from itertools import chain


# Fast find match using `find`
# Find overlap of at least min_overlap of prefix of s2 in s1
# Stores an index to search for next match
def find_overlap(s1, s2, min_overlap=None):
    ix = 1
    if min_overlap is None:
        min_overlap = floor(len(s2) / 2)
    while ix < len(s1):
        ix = s1.find(s2[:min_overlap], ix)
        if ix == -1:
            break
        if s2.startswith(s1[ix:]):
            return len(s1) - ix
        ix += 1


def construct_assembly(seqs):
    # Build a forward and reverse overlap graph of sequences
    fmap, rmap, starts, ends = ({}, {}, {}, {})
    for p1 in seqs.keys():
        for p2 in seqs.keys():
            if p1 in starts or p2 in ends or p1 in p2:
                continue
            n = find_overlap(seqs[p1], seqs[p2])
            if n:
                fmap[p1] = {"overlap": n, "next": p2}
                rmap[p2] = p1
                starts[p1] = True
                ends[p2] = True
                break

    # Find starting key using rmap
    k = list(seqs.keys())[0]
    while k in rmap:
        k = rmap[k]

    # Initialise list with sequence 1, then add suffixes of matching sequences
    seq = [seqs[k]]
    while k in fmap:
        seq.append(seqs[fmap[k]["next"]][fmap[k]["overlap"] :])
        k = fmap[k]["next"]

    # Join all sequences
    return "".join(seq)


# Find sequences with errors as those that are unique
# The remainder (those that are not unique) are assumed to be correct
# To find how to fix our unique sequences, we find the non-unique version
# that is one hamming distance away. This ignores the possibility that two
# unique seqeunces may be from the same read with just one of them having an
# error. In this case though, we would not be able to determine which one is
# incorrect anyway...
def find_errors(seqs):
    rseqs = [Dna(x).revc().seq for x in seqs]
    unique = [x for x in seqs if seqs.count(x) + rseqs.count(x) == 1]
    correct = set(seqs + rseqs).difference(set(unique))

    for x in unique:
        for y in correct:
            if aln.hamm(x, y) == 1:
                yield x + "->" + y


def dbru(seqs):
    seqs = list(seqs)
    seqs = set(seqs).union([Dna(seq).revc().seq for seq in seqs])
    return [(x[:-1], x[1:]) for x in seqs]


def find_cycle(graph):
    s = sorted(list(graph.keys()))[0]
    visited = set()
    cycle = []
    while s not in visited:
        cycle += [s]
        visited.add(s)
        if s not in graph:
            return None
        s = graph[s]
    return cycle


def join_cycle(chain):
    return "".join(x[0] for x in chain)


def pcov(seqs):
    x = find_cycle(dict(dbru(seqs)))
    return join_cycle(sorted(x))


# Return all kmers of length k from sequence seq
def kmers(seq, k):
    return [seq[i : (i + k)] for i in range(len(seq) - k + 1)]


def extract_chain(graph):
    ch1 = find_cycle(graph)
    if ch1 is None:
        return None, None
    graph = dict(filter(lambda x: x[0] not in ch1, graph.items()))
    return ch1, graph


def gasm(seqs):
    s0 = seqs[0][:-1]
    for k in range(len(s0) + 1, 0, -1):
        subseqs = chain(*[kmers(x, k) for x in seqs])
        graph = dict(dbru(subseqs))
        ch1, graph = extract_chain(graph)
        if ch1 is None:
            continue

        ch2, graph = extract_chain(graph)
        if ch2 is None:
            continue
        if len(graph) == 0:
            return join_cycle(ch1)
