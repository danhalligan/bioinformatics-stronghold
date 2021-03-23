import re
import pkg_resources as pr
import yaml
import string


class Parser:
    """Parse problem data"""

    def __init__(self, file):
        self.file = file

    def dna(self):
        """Return the first line as a DNA string"""
        return Dna(self.line())

    def line(self):
        return open(self.file).readline().rstrip()

    def lines(self):
        return open(self.file).read().splitlines()

    def fastas(self):
        return list(read_fasta(open(self.file, "r")))

    def seqs(self):
        return [x.seq for x in self.fastas()]

    def ints(self):
        return list(map(int, self.line().split()))

    def floats(self):
        return list(map(float, self.line().split()))


class Rec:
    """A FASTA record"""

    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

    def __len__(self):
        return len(self.seq)


class Seq:
    """Sequence methods (either DNA, RNA or Protein)"""

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
        code = genetic_code()
        x = self.rna().seq
        prot = [code[x[i : i + 3]] for i in range(0, len(x), 3)]
        prot = "".join(prot)
        prot = Prot(re.sub("\\*$", "", prot))
        return prot

    @property
    def alphabet(self):
        return "ACGT"


class Rna(Seq):
    """An RNA sequence"""

    @property
    def alphabet(self):
        return "ACGU"


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


def genetic_code():
    stream = pr.resource_stream(__name__, "data/genetic_code.yaml")
    return yaml.load(stream, Loader=yaml.FullLoader)


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
