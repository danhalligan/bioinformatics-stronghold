import inspect
from Bio import SeqIO

class Parser:
    """Parse problem data"""
    
    def __init__(self, file):
        self.file = file
    
    def line(self):
        return open(self.file).readline().rstrip()

    def lines(self):
        return open(self.file).read().splitlines()

    def fastas(self):
        return list(SeqIO.parse(self.file, "fasta"))
