import inspect
from Bio import SeqIO

class Data:
    """Get problem data"""
    
    def __init__(self, fun, test=False):
        if test:
            self.file = "tests/data/test_" + fun + ".txt"
        else:
            self.file = "data/rosalind_" + fun + ".txt"
    
    def line(self):
        return open(self.file).readline().rstrip()

    def lines(self):
        return open(self.file).read().splitlines()

    def fastas(self):
        return list(SeqIO.parse(self.file, "fasta"))


class Parser:
    """Get problem data"""
    
    def __init__(self, file):
        self.file = file
    
    def line(self):
        return open(self.file).readline().rstrip()

    def lines(self):
        return open(self.file).read().splitlines()

    def fastas(self):
        return list(SeqIO.parse(self.file, "fasta"))
