import re


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


class Rec:
    """A FASTA record"""

    def __init__(self, id, seq):
        self.id = id
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def seq(self):
        return self.seq


class Seq:
    """Sequence methods (either DNA, RNA or Protein)"""

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def __eq__(self, other):
        return self.seq == other

    def __getitem__(self, value):
        return self.seq.__getitem__(value)

    def __repr__(self):
        return self.__class__.__name__ + "(" + self.seq + ")"

    def table(self):
        """Count number of each letter"""
        return {x: self.seq.count(x) for x in self.alphabet}


class Dna(Seq):

    """A DNA sequence"""

    def __init__(self, seq: str):
        self.alphabet = "ACGT"
        self.seq = seq

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
        x = self.rna()
        prot = [code[x[i : i + 3]] for i in range(0, len(x), 3)]
        prot = "".join(prot)
        prot = re.sub("\\*$", "", prot)
        return prot


class Rna(Seq):

    """An RNA sequence"""

    def __init__(self, seq):
        self.alphabet = "ACGU"
        self.seq = seq


class Prot(Seq):

    """A protein sequence"""

    def __init__(self, seq):
        self.alphabet = "*ACDEFGHIKLMNPQRSTVWY"
        self.seq = seq


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
    return {
        "UUU": "F",
        "UCU": "S",
        "UAU": "Y",
        "UGU": "C",
        "UUC": "F",
        "UCC": "S",
        "UAC": "Y",
        "UGC": "C",
        "UUA": "L",
        "UCA": "S",
        "UAA": "*",
        "UGA": "*",
        "UUG": "L",
        "UCG": "S",
        "UAG": "*",
        "UGG": "W",
        "CUU": "L",
        "CCU": "P",
        "CAU": "H",
        "CGU": "R",
        "CUC": "L",
        "CCC": "P",
        "CAC": "H",
        "CGC": "R",
        "CUA": "L",
        "CCA": "P",
        "CAA": "Q",
        "CGA": "R",
        "CUG": "L",
        "CCG": "P",
        "CAG": "Q",
        "CGG": "R",
        "AUU": "I",
        "ACU": "T",
        "AAU": "N",
        "AGU": "S",
        "AUC": "I",
        "ACC": "T",
        "AAC": "N",
        "AGC": "S",
        "AUA": "I",
        "ACA": "T",
        "AAA": "K",
        "AGA": "R",
        "AUG": "M",
        "ACG": "T",
        "AAG": "K",
        "AGG": "R",
        "GUU": "V",
        "GCU": "A",
        "GAU": "D",
        "GGU": "G",
        "GUC": "V",
        "GCC": "A",
        "GAC": "D",
        "GGC": "G",
        "GUA": "V",
        "GCA": "A",
        "GAA": "E",
        "GGA": "G",
        "GUG": "V",
        "GCG": "A",
        "GAG": "E",
        "GGG": "G",
    }
