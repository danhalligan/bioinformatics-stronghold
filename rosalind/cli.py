import typer
import rosalind.rosalind as ros
import rosalind.mass as mass
import rosalind.assembly as assembly
from rosalind.helpers import Parser

app = typer.Typer()


@app.command("dna")
def dna(file: str):
    """Counting DNA Nucleotides"""
    print(*list(Parser(file).dna().table().values()))


@app.command("rna")
def rna(file: str):
    """Transcribing DNA into RNA"""
    print(Parser(file).dna().rna())


@app.command("revc")
def revc(file: str):
    """Complementing a Strand of DNA"""
    print(Parser(file).dna().revc())


@app.command("fib")
def fib(file: str):
    """Rabbits and Recurrence Relations"""
    n, k = map(int, Parser(file).line().split())
    print(ros.rabbits(n, k))


@app.command("gc")
def gc(file: str):
    """Computing GC Content"""
    res = ros.max_gc(Parser(file).fastas())
    print(res["name"], res["value"], sep="\n")


@app.command("hamm")
def hamm(file: str):
    """Counting Point Mutations"""
    s1, s2 = Parser(file).lines()
    print(ros.hamming(s1, s2))


@app.command("iprb")
def iprb(file: str):
    """Mendel's First Law"""
    k, m, n = map(int, Parser(file).line().split())
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
def fibd(file: str):
    """Mortal Fibonacci Rabbits"""
    n, m = map(int, Parser(file).line().split())
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
    seq = ros.Dna(Parser(file).fastas()[0].seq)
    orfs = list(ros.find_orfs(seq))
    orfs = sorted(list(dict.fromkeys(orfs)))
    print("\n".join(orfs))


@app.command("perm")
def perm(file: str):
    """Enumerating Gene Orders"""
    from itertools import permutations

    n = int(Parser(file).line().split()[0])
    perm = list(permutations(range(1, n + 1)))
    print(len(perm))
    for i in perm:
        print(*i)


@app.command("prtm")
def prtm(file: str):
    """Calculating Protein Mass"""
    print(mass.protein_mass(Parser(file).line()))


@app.command("revp")
def revp(file: str):
    """Locating Restriction Sites"""
    seq = Parser(file).fastas()[0].seq
    res = list(ros.reverse_pallindromes(seq))
    res.sort()
    for row in res:
        print(*row, sep=" ")


@app.command("splc")
def splc(file: str):
    """RNA Splicing"""
    from functools import reduce

    seqs = [x.seq for x in Parser(file).fastas()]

    def trim(gene, intron):
        s = gene.find(intron)
        return gene[:s] + gene[s + len(intron) :]

    print(ros.Dna(reduce(trim, seqs)).translate())


@app.command("lexf")
def lexf(file: str):
    """Enumerating k-mers Lexicographically"""
    from itertools import product

    l1, l2 = Parser(file).lines()
    set = l1.split(" ")
    n = int(l2)
    perm = list(product(set, repeat=n))
    perm = ["".join(x) for x in perm]
    print("\n".join(sorted(perm)))


@app.command("pper")
def pper(file: str):
    """Partial Permutations"""
    from functools import reduce

    n, k = map(int, Parser(file).line().split())
    print(reduce(lambda p, i: int((p * i) % 1e6), range(n, n - k, -1)))


@app.command("prob")
def prob(file: str):
    """Introduction to Random Strings"""
    import math

    def logp(c, p):
        return math.log10({"G": p, "C": p, "A": 1 - p, "T": 1 - p}[c] / 2)

    seq, arr = Parser(file).lines()
    arr = [float(x) for x in arr.split(" ")]
    res = [round(sum([logp(c, p) for c in seq]), 3) for p in arr]
    print(*res)


@app.command("sseq")
def sseq(file: str):
    """Finding a Spliced Motif"""

    def matches(s1, s2):
        i, j = 0, 0
        while i < len(s1) and j < len(s2):
            if s2[j] == s1[i]:
                yield (i + 1)
                j += 1
            i += 1

    s1, s2 = [x.seq for x in Parser(file).fastas()]
    print(*list(matches(s1, s2)))


@app.command("tran")
def tran(file: str):
    """Transitions and Transversions"""

    def ts(x, y):
        return (x == "A" and y == "G") or (x == "C" and y == "T")

    seqs = [x.seq for x in Parser(file).fastas()]
    mm, tr = 0, 0
    for x, y in zip(seqs[0], seqs[1]):
        if x != y:
            mm += 1
            tr += int(ts(x, y) or ts(y, x))

    print(tr / (mm - tr))


def kmer_perm(k):
    from itertools import product

    set = ["A", "C", "G", "T"]
    perm = list(product(set, repeat=k))
    return ["".join(x) for x in perm]


@app.command("kmer")
def kmer(file: str):
    """k-Mer Composition"""

    seq = Parser(file).fastas()[0].seq

    # initialise hash with all possible 4-mers permutations
    d = {k: 0 for k in kmer_perm(4)}

    # Run through 4-mer slices of sequence and increment dictionary keys
    for i in range(len(seq) - 3):
        d[seq[i : (i + 4)]] += 1

    print(*d.values())


@app.command("kmp")
def kmp(file: str):
    """Speeding Up Motif Finding"""
    seq = Parser(file).fastas()[0].seq
    print(*ros.kmp_preprocess(seq))


@app.command("rstr")
def rstr(file: str):
    import math

    l1, seq = Parser(file).lines()
    n, x = map(float, l1.split(" "))
    gc = sum([seq.count(x) for x in "GC"])
    lam = ((1 - x) / 2) ** (len(seq) - gc) * (x / 2) ** gc * n
    print(1 - math.exp(-lam))


@app.command("spec")
def spec(file: str):
    """Inferring Protein from Spectrum"""

    weights = [float(x) for x in Parser(file).lines()]
    diff = [j - i for i, j in zip(weights[:-1], weights[1:])]
    print("".join([mass.match_mass(x) for x in diff]))


@app.command("grph")
def grph(file: str):
    """Overlap Graphs"""
    fa = Parser(file).fastas()
    out = list(ros.overlap_graph(fa))
    for i in out:
        print(*i)


@app.command("tree")
def tree(file: str):
    """Completing a Tree"""
    # Nb. A connected tree of n nodes will always contain n-1 edge
    data = Parser(file).lines()
    n_nodes = int(data[0])
    print(n_nodes - len(data[1:]) - 1)


@app.command("long")
def long(file: str):
    """Genome Assembly as Shortest Superstring"""
    seqs = Parser(file).fastas()
    print(assembly.construct_assembly(seqs))


def main():
    app()
