import typer
import rosalind.rosalind as ros
from rosalind.helpers import Parser

app = typer.Typer()


@app.command("dna")
def dna(file: str):
    """Counting DNA Nucleotides"""
    print(Parser(file).dna().count_bases)


@app.command("rna")
def rna(file: str):
    """Transcribing DNA into RNA"""
    print(ros.dna2rna(Parser(file).line()))


@app.command("revc")
def revc(file: str):
    """Complementing a Strand of DNA"""
    print(ros.revcomp(Parser(file).line()))


@app.command("fib")
def fib(n: int, k: int):
    """Rabbits and Recurrence Relations"""
    print(ros.rabbits(n, k))


@app.command("gc")
def gc(file: str):
    """Computing GC Content"""
    x = Parser(file).fastas()
    res = ros.max_gc(x)
    print(res["name"], res["value"], sep="\n")


@app.command("hamm")
def hamm(file: str):
    """Counting Point Mutations"""
    s1, s2 = Parser(file).lines()
    print(ros.hamming(s1, s2))


@app.command("iprb")
def iprb(k: int, m: int, n: int):
    """Mendel's First Law"""
    print(ros.mendel1(k, m, n))


@app.command("prot")
def prot(file: str):
    """Translating RNA into Protein"""
    seq = Parser(file).line()
    print(ros.translate(seq))


@app.command("subs")
def subs(file: str):
    """Finding a Motif in DNA"""
    s1, s2 = Parser(file).lines()
    print(*list(ros.find_motif(s1, s2)))


@app.command("cons")
def cons(file: str):
    """Consensus and Profile"""
    x = [x.seq for x in Parser(file).fastas()]
    mat = ros.profile_matrix(x)
    cons = ros.consensus_sequence(mat)
    print(cons)
    for i in range(0, 4):
        print("ACGT"[i] + ":", *mat[i])


@app.command("fibd")
def fibd(n: int, m: int):
    """Mortal Fibonacci Rabbits"""
    print(ros.mortal_rabbits(n, m))


@app.command("iev")
def iev(file: str):
    """Calculating Expected Offspring"""
    v = map(int, Parser(file).line().split())
    print(ros.expected_offspring(v))


@app.command("lcsm")
def lcsm(file: str):
    """Finding a Shared Motif"""
    seqs = [x.seq for x in Parser(file).fastas()]
    print(ros.find_shared_motif(seqs))


@app.command("lia")
def lia(file: str):
    """Independent Alleles"""
    k, n = map(int, Parser(file).line().split())
    print(ros.mendel2(k, n))


@app.command("mprt")
def mprt(file: str):
    """Finding a Protein Motif"""
    for id in Parser(file).lines():
        seq = ros.get_uniprot(id)
        matches = ros.find_protein_motif(str(seq.seq))
        if len(matches):
            print(id)
            print(*matches)


@app.command("mrna")
def mrna(file: str):
    """Inferring mRNA from Protein"""
    print(ros.count_rnas(Parser(file).line()))


@app.command("orf")
def orf(file: str):
    """Open Reading Frames"""


@app.command("perm")
def perm(n: int):
    """Enumerating Gene Orders"""
    from itertools import permutations

    perm = list(permutations(range(1, n + 1)))
    print(len(perm))
    for i in perm:
        print(*i)


@app.command("prtm")
def prtm(file: str):
    """Calculating Protein Mass"""
    print(ros.protein_mass(Parser(file).line()))


@app.command("kmp")
def kmp(file: str):
    """Shortening the Motif Search"""
    seq = Parser(file).fastas()[0].seq
    print(*ros.kmp_preprocess(seq))


@app.command("revp")
def revp(file: str):
    """Locating Restriction Sites"""
    seq = Parser(file).fastas()[0].seq
    res = list(ros.reverse_pallindromes(seq))
    res.sort()
    for row in res:
        print(*row, sep=" ")


def main():
    app()
