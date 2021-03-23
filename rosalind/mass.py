from rosalind.helpers import aa_mass
from collections import Counter


def protein_mass(prot):
    aam = aa_mass()
    return sum([aam[x] for x in prot])


def match_mass(weight):
    aam = aa_mass()
    dist = [abs(weight - x) for x in aam.values()]
    return list(aam.keys())[dist.index(min(dist))]


def conv(s1, s2):
    """Comparing Spectra with the Spectral Convolution"""
    x = sorted([round(i - j, 5) for i in s1 for j in s2])
    return Counter(x).most_common()[0]
