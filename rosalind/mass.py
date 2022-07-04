from rosalind.helpers import aa_mass
from collections import Counter
from collections import defaultdict
from math import isclose


def protein_mass(prot):
    aam = aa_mass()
    return sum([aam[x] for x in prot])


def match_mass(weight, rel_tol=1e-6):
    aam = aa_mass()
    matches = [k for k in aam if isclose(weight, aam[k], rel_tol=rel_tol)]
    return None if len(matches) == 0 else matches[0]


def conv(s1, s2):
    """Comparing Spectra with the Spectral Convolution"""
    x = sorted([round(i - j, 5) for i in s1 for j in s2])
    return Counter(x).most_common()[0]


def spectrum(seq):
    """Complete spectrum"""
    r = range(1, len(seq))
    spec = [seq] + [seq[:k] for k in r] + [seq[k:] for k in r]
    return [protein_mass(x) for x in spec]


def prsm(s, r):
    """Matching a Spectrum to a Protein"""
    mult = [conv(spectrum(x), r)[1] for x in s]
    return [max(mult), s[mult.index(max(mult))]]


# Return a dictionary that let's us look up complementary ion pairs
# e.g. "PRO" and "TEIN" for full peptide "PROTEIN"
# We can find pairs, because their weights approximately sum to that
# of the peptide.
def find_pairs(peptide, ions):
    pairs = {}
    for w in ions:
        for w2 in ions:
            if isclose(w2 + w, peptide):
                pairs[w] = w2
    return pairs


# This builds a graph of possibly adjacent ions.
# ions can be considered adjacent in a graph if the difference in
# their masses is close to the mass of a known amino acid.
def weight_graph(ions):
    graph = defaultdict(list)
    for i in range(len(ions)):
        for j in range(i + 1, len(ions)):
            aa = match_mass(ions[j] - ions[i])
            if aa:
                graph[ions[i]] += [[ions[j], aa]]
    return graph


def full(peptide, ions):
    """Inferring Peptide from Full Spectrum"""

    def infer_peptide(w, seq, rm):
        for w2, aa in graph[w]:
            if w2 in rm:
                continue
            if len(seq) + 1 == target_len:
                yield seq + aa
            else:
                yield from infer_peptide(w2, seq + aa, rm + [w, pairs[w]])

    graph = weight_graph(ions)
    pairs = find_pairs(peptide, ions)
    target_len = int(len(ions) / 2 - 1)
    w = min(ions)
    return list(infer_peptide(w, "", [w, pairs[w]]))


def sgra(ions):
    """Using the Spectrum Graph to Infer Peptides"""

    def infer_peptide(w, seq):
        for w2, aa in graph[w]:
            yield from infer_peptide(w2, seq + aa)
        yield seq

    graph = weight_graph(ions)
    w = min(ions)
    return max(list(infer_peptide(w, "")), key=len)
