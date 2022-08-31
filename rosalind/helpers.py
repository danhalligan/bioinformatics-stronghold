import re
import pkg_resources as pr
import yaml
import string
import sys


class Parser:
    """Parse problem data text files"""

    def __init__(self, file):
        self.file = file

    def dna(self):
        """Return the first line as a DNA string"""
        return Dna(self.line())

    def rna(self):
        """Return the first line as a RNA string"""
        return Rna(self.line())

    def line(self):
        """Return the first line"""
        return open(self.file).readline().rstrip()

    def lines(self):
        """Return lines as a chomped list"""
        return open(self.file).read().splitlines()

    def fastas(self):
        """Return fasta records as a list"""
        return list(read_fasta(open(self.file, "r")))

    def seqs(self):
        """Return sequences from fasta records as a list"""
        return [x.seq for x in self.fastas()]

    def ints(self):
        """Return space separated integers from first line"""
        return list(map(int, self.line().split()))

    def floats(self):
        """Return space separated floats from first line"""
        return list(map(float, self.line().split()))


class Rec:
    """A simple FASTA record"""

    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

    def __len__(self):
        return len(self.seq)


class Seq:
    """A sequence with a specified alphabet (either DNA, RNA or Protein)"""

    def __init__(self, seq: str):
        if 0 in [c in self.alphabet for c in seq]:
            raise TypeError("String contains invalid characters!")
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def __eq__(self, other):
        return self.seq == other

    def __getitem__(self, value):
        return type(self)(self.seq.__getitem__(value))

    def __repr__(self):
        return self.__class__.__name__ + "(" + self.seq + ")"

    def table(self):
        """Count number of each letter"""
        return {x: self.seq.count(x) for x in self.alphabet}

    @property
    def alphabet(self):
        return string.ascii_uppercase


class Dna(Seq):
    """A DNA sequence"""

    def rna(self):
        """Convert DNA to RNA"""
        return Rna(re.sub("T", "U", self.seq))

    def revc(self):
        """Reverse complement"""
        return Dna(self.seq[::-1].translate(str.maketrans("ACGT", "TGCA")))

    def gc_content(self):
        """Calculate GC content"""
        return sum([self.seq.count(base) for base in "GC"]) / len(self)

    def translate(self):
        """Translate DNA to protein sequence. The stop codon is removed
        automatically."""
        return self.rna().translate()

    @property
    def alphabet(self):
        return "ACGT"


class Rna(Seq):
    """An RNA sequence"""

    def __init__(self, seq: str):
        super().__init__(seq)
        self._code = self._genetic_code()

    @property
    def alphabet(self):
        return "ACGU"

    def translate(self):
        """Translate RNA to protein sequence. The stop codon is removed
        automatically."""
        x = self.seq
        prot = [self._code[x[i : i + 3]] for i in range(0, len(x), 3)]
        prot = "".join(prot)
        return Prot(re.sub("\\*$", "", prot))

    def _genetic_code(self):
        stream = pr.resource_stream(__name__, "data/genetic_code.yaml")
        return yaml.load(stream, Loader=yaml.FullLoader)


class Prot(Seq):
    """A protein sequence"""

    @property
    def alphabet(self):
        return "*ACDEFGHIKLMNPQRSTVWY"


def read_fasta(handle):
    header, sequence = "", []
    for line in handle:
        if line[0] == ">":
            if sequence:
                yield Rec(header, "".join(sequence))
            header, sequence = line[1:-1], []
        else:
            sequence.append(line.strip())
    yield Rec(header, "".join(sequence))


def codons():
    stream = pr.resource_stream(__name__, "data/codons.yaml")
    return yaml.load(stream, Loader=yaml.FullLoader)


def aa_mass():
    stream = pr.resource_stream(__name__, "data/aa_mass.yaml")
    return yaml.load(stream, Loader=yaml.FullLoader)


def memoize(f):
    cache = {}

    def wrapper(*args):
        if args not in cache:
            cache[args] = f(*args)
        return cache[args]

    return wrapper


@memoize
def blosum62():
    lines = pr.resource_string(__name__, "data/blosum62.txt").decode().split("\n")
    header = lines[0].split()
    return dict([x[0], dict(zip(header, map(int, x.split()[1:])))] for x in lines[1:])


@memoize
def pam250():
    lines = pr.resource_string(__name__, "data/pam250.txt").decode().split("\n")
    header = lines[0].split()
    return dict([x[0], dict(zip(header, map(int, x.split()[1:])))] for x in lines[1:])


class recursionlimit:
    def __init__(self, limit):
        self.limit = limit

    def __enter__(self):
        self.old_limit = sys.getrecursionlimit()
        sys.setrecursionlimit(self.limit)

    def __exit__(self, type, value, tb):
        sys.setrecursionlimit(self.old_limit)
