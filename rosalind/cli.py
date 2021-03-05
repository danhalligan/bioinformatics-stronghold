import typer

from .rosalind import *
from .helpers import *

app = typer.Typer()

@app.command("dna")
def dna(file: str):
    """Counting DNA Nucleotides
    """
    print(*count_nucleotides(Parser(file).line()))


@app.command("rna")
def rna(file: str):
    """Transcribing DNA into RNA
    """
    print(dna2rna(Parser(file).line()))


@app.command("revc")
def revc(file: str):
    """Complementing a Strand of DNA
    """
    print(revcomp(Parser(file).line()))


@app.command("fib")
def fib(file: str):
    """Rabbits and Recurrence Relations
    """
    n, k = map(int, Parser(file).line())
    print(fibonacci(n, k))


@app.command("gc")
def gc(file: str):
    """Computing GC Content
    """
    x = Parser(file).fastas()
    res = max_gc(x)
    print(res['name'], res['value'], sep="\n")


@app.command("hamm")
def hamm(file: str):
    """Counting Point Mutations
    """
    s1, s2 = Parser(file).lines()
    print(hamming_distance(s1, s2))


@app.command("iprb")
def iprb(file: str):
    """Mendel's First Law
    """
    k, m, n = map(int, Parser(file).line().split())
    print(mendel1(k, m, n))


@app.command("prot")
def prot(file: str):
    """Translating RNA into Protein
    """
    seq = Parser(file).line()
    print(translate(seq))


@app.command("subs")
def subs(file: str):
    """Finding a Motif in DNA
    """
    s1, s2 = Parser(file).lines()
    print(*list(find_motif(s1, s2)))


@app.command("cons")
def cons(file: str):
    """Consensus and Profile
    """
    x = [x.seq for x in Parser(file).fastas()]
    mat = profile_matrix(x)
    cons = consensus_sequence(mat)
    print(cons)
    for i in range(0, 4):
        print('ACGT'[i] + ':', *mat[i])


@app.command("fibd")
def fibd(file: str):
    """Mortal Fibonacci Rabbits
    """
    n, m = map(int, Parser(file).line().split())
    print(mortal_rabbits(n, m))


@app.command("iev")
def iev(file: str):
    """Calculating Expected Offspring
    """
    v = map(int, Parser(file).line().split())
    print(expected_offspring(v))


@app.command("lcsm")
def lcsm(file: str):
    """Finding a Shared Motif
    """
    seqs = [x.seq for x in Parser(file).fastas()]
    print(find_shared_motif(seqs))


@app.command("lia")
def lia(file: str):
    """Independent Alleles
    """
    k, n = map(int, Parser(file).line().split())
    print(mendel2(k, n))

 
@app.command("mprt")
def mprt(file: str):
    """Finding a Protein Motif
    """
    for id in Parser(file).lines():
        seq = get_uniprot(id)
        matches = find_protein_motif(str(seq.seq))
        if len(matches):
            print(id)
            print(*matches)


@app.command("mrna")
def mrna(file: str):
    """Inferring mRNA from Protein
    """
    print(count_rnas(Parser(file).line()))


@app.command("orf")
def orf(file: str):
    """Open Reading Frames
    """

@app.command("kmp")
def kmp(file: str):
    """Shortening the Motif Search
    """
    seq = list(SeqIO.parse("data/rosalind_kmp.txt", "fasta"))[0].seq
    print(*kmp_preprocess(seq))


def main():
    app()
