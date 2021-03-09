from math import floor

# This builds an overlap graph by considering all pairs of sequences and
# looking for suffixes matching prefixes. Then using this graph we reconstruct
# the sequence.


def find_match(s1, s2):
    for n in range(floor(len(s2) / 2), len(s1)):
        if s1.endswith(s2[:n]):
            return n


# Speed increases here can be obtained by not searching for suffix matches in
# sequences where we've already found a suffix (starts) and the same for
# prefixes (ends). Ending early like this uses ~1/3 of calls to find_match
def overlap_graph(seqs):
    """Build an overlap graph of sequences where suffix of A matches prefix of
    B"""
    fmap = {}
    rmap = {}
    ids = seqs.keys()
    starts = {}
    ends = {}
    for p1 in ids:
        if p1 in starts:
            continue
        for p2 in ids:
            if p2 in ends:
                continue
            n = find_match(seqs[p1], seqs[p2])
            if n:
                fmap[p1] = {"n": n, "l": p2}
                rmap[p2] = p1
                starts[p1] = True
                ends[p2] = True
                break
    return {"map": fmap, "revmap": rmap}


def construct_assembly(seqs):
    seqs = {x.id: x.seq for x in seqs}
    maps = overlap_graph(seqs)

    # Find starting key
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
