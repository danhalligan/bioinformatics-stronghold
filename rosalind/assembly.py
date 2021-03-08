from itertools import permutations
from math import floor

# This builds an overlap graph by considering all pairs of sequences and
# looking for suffixes matching prefixes. Then using this graph we reconstruct
# the sequence. This could be made faster I think by not considering all
# permutations and possibly constructing the sequence by collapsing pairs as
# we find them?


def overlap_graph(seqs):
    """Build an overlap graph of sequences where suffix of A matches prefix of
    B"""
    fmap = {}
    rmap = {}
    for pair in permutations(seqs, 2):
        for n in range(floor(len(pair[0]) / 2), len(pair[0])):
            if pair[0].seq.endswith(pair[1].seq[:n]):
                fmap[pair[0].id] = {"n": n, "l": pair[1].id}
                rmap[pair[1].id] = pair[0].id
    return {"map": fmap, "revmap": rmap}


def construct_assembly(seqs):
    maps = overlap_graph(seqs)

    # Find starting key
    seqs = {x.id: x.seq for x in seqs}
    k = list(seqs.keys())[0]
    while k in maps["revmap"]:
        k = maps["revmap"][k]

    # Initialise with first sequence, and then build sequence by appending tails
    # of next sequence in graph
    seq = [seqs[k]]
    g = maps["map"]
    while k in g:
        s = seqs[g[k]["l"]]
        seq.append(s[g[k]["n"] : len(s)])
        k = g[k]["l"]
    return "".join(seq)
