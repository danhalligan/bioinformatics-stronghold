import typer
import math
import builtins
import rosalind.rosalind as ros
import rosalind.mass as mass
import rosalind.assembly as assembly
import rosalind.alignment as aln
import rosalind.combinatorics as com
import rosalind.graph as graph
import rosalind.probability as pr

from rosalind.helpers import Parser
from itertools import permutations, product
from functools import reduce

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
    print(com.fib(*Parser(file).ints()))


@app.command("gc")
def gc(file: str):
    """Computing GC Content"""
    res = ros.max_gc(Parser(file).fastas())
    print(res["name"], res["value"], sep="\n")


@app.command("hamm")
def hamm(file: str):
    """Counting Point Mutations"""
    print(aln.hamm(*Parser(file).lines()))


@app.command("iprb")
def iprb(file: str):
    """Mendel's First Law"""
    print(pr.iprb(*Parser(file).ints()))


@app.command("prot")
def prot(file: str):
    """Translating RNA into Protein"""
    print(Parser(file).rna().translate())


@app.command("subs")
def subs(file: str):
    """Finding a Motif in DNA"""
    s1, s2 = Parser(file).lines()
    print(*list(ros.find_motif(s1, s2)))


@app.command("cons")
def cons(file: str):
    """Consensus and Profile"""
    x = Parser(file).seqs()
    mat = ros.profile_matrix(x)
    print(ros.consensus_sequence(mat))
    for i in range(0, 4):
        print("ACGT"[i] + ":", *mat[i])


@app.command("fibd")
def fibd(file: str):
    """Mortal Fibonacci Rabbits"""
    print(com.fibd(*Parser(file).ints()))


@app.command("grph")
def grph(file: str):
    """Overlap Graphs"""
    fa = Parser(file).fastas()
    for i in graph.overlap_graph(fa):
        print(*i)


@app.command("iev")
def iev(file: str):
    """Calculating Expected Offspring"""
    print(pr.iev(Parser(file).ints()))


@app.command("lcsm")
def lcsm(file: str):
    """Finding a Shared Motif"""
    print(ros.lcsm(Parser(file).seqs()))


@app.command("lia")
def lia(file: str):
    """Independent Alleles"""
    print(pr.lia(*Parser(file).ints()))


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
    print(com.mrna(Parser(file).line()))


@app.command("orf")
def orf(file: str):
    """Open Reading Frames"""
    seq = ros.Dna(Parser(file).fastas()[0].seq)
    orfs = list(com.orf(seq))
    orfs = sorted(list(dict.fromkeys(orfs)))
    print("\n".join(orfs))


@app.command("perm")
def perm(file: str):
    """Enumerating Gene Orders"""
    n = Parser(file).ints()[0]
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
    res = sorted(ros.reverse_pallindromes(seq))
    for row in res:
        print(*row)


@app.command("splc")
def splc(file: str):
    """RNA Splicing"""

    def trim(gene, intron):
        s = gene.find(intron)
        return gene[:s] + gene[s + len(intron) :]

    seqs = Parser(file).seqs()
    print(ros.Dna(reduce(trim, seqs)).translate())


@app.command("lexf")
def lexf(file: str):
    """Enumerating k-mers Lexicographically"""

    l1, l2 = Parser(file).lines()
    set = l1.split(" ")
    n = int(l2)
    perm = ["".join(x) for x in product(set, repeat=n)]
    print(*sorted(perm), sep="\n")


@app.command("lgis")
def lgis(file: str):
    """Longest Increasing Subsequence"""
    data = Parser(file).lines()[1]
    data = [int(x) for x in data.split(" ")]
    s1 = ros.lgis(data)
    print(*s1)

    s2 = ros.lgis([-x for x in data])
    s2 = [-x for x in s2]
    print(*s2)


@app.command("long")
def long(file: str):
    """Genome Assembly as Shortest Superstring"""
    seqs = Parser(file).fastas()
    seqs = dict([(x.id, x.seq) for x in seqs])
    print(assembly.construct_assembly(seqs))


@app.command("pper")
def pper(file: str):
    """Partial Permutations"""
    n, k = Parser(file).ints()
    print(reduce(lambda p, i: (p * i) % 1_000_000, range(n, n - k, -1)))


@app.command("prob")
def prob(file: str):
    """Introduction to Random Strings"""
    seq, arr = Parser(file).lines()
    arr = [float(x) for x in arr.split()]
    print(*pr.prob(seq, arr))


@app.command("sign")
def sign(file: str):
    """Enumerating Oriented Gene Orderings"""
    n = Parser(file).ints()[0]
    res = com.sign(n)
    print(len(res))
    for i in res:
        print(*i)


@app.command("sseq")
def sseq(file: str):
    """Finding a Spliced Motif"""

    def matches(s1, s2):
        i, j = 0, 0
        while j < len(s2):
            if s2[j] == s1[i]:
                yield i + 1
                j += 1
            i += 1

    s1, s2 = [x.seq for x in Parser(file).fastas()]
    print(*list(matches(s1, s2)))


@app.command("tran")
def tran(file: str):
    """Transitions and Transversions"""
    print(aln.tran(Parser(file).seqs()))


@app.command("tree")
def tree(file: str):
    """Completing a Tree"""
    # Nb. A connected tree of n nodes will always contain n-1 edges
    data = Parser(file).lines()
    n_nodes = int(data[0])
    print(n_nodes - len(data[1:]) - 1)


@app.command("corr")
def corr(file: str):
    """Error Correction in Reads"""
    seqs = Parser(file).fastas()
    seqs = [x.seq for x in seqs]
    print(*assembly.find_errors(seqs), sep="\n")


@app.command("kmer")
def kmer(file: str):
    """k-Mer Composition"""

    seq = Parser(file).fastas()[0].seq

    # initialise hash with all possible 4-mers permutations
    d = {k: 0 for k in ros.kmer_perm(4)}

    # Run through 4-mer slices of sequence and increment dictionary keys
    for i in range(len(seq) - 3):
        d[seq[i : (i + 4)]] += 1

    print(*d.values())


@app.command("dbru")
def dbru(file: str):
    """Constructing a De Bruijn Graph"""
    for pair in sorted(assembly.dbru(Parser(file).lines())):
        print("(", pair[0], ", ", pair[1], ")", sep="")


@app.command("pcov")
def pcov(file: str):
    """Genome Assembly with Perfect Coverage"""
    print(assembly.pcov(Parser(file).lines()))


@app.command("kmp")
def kmp(file: str):
    """Speeding Up Motif Finding"""
    seq = Parser(file).fastas()[0].seq
    print(*ros.kmp_preprocess(seq))


@app.command("rear")
def rear(file: str):
    """Reversal Distance"""
    data = open(file).read().strip().split("\n\n")
    data = [tuple([list(map(int, y.split())) for y in x.split("\n")]) for x in data]
    print(*[com.get_distance(s, t) for s, t in data])


@app.command("rstr")
def rstr(file: str):
    """Matching Random Motifs"""
    l1, seq = Parser(file).lines()
    n, x = map(float, l1.split(" "))
    gc = sum([seq.count(x) for x in "GC"])
    lam = ((1 - x) / 2) ** (len(seq) - gc) * (x / 2) ** gc * n
    print(1 - math.exp(-lam))


@app.command("sset")
def sset(file: str):
    """Counting Subsets"""
    n = Parser(file).ints()[0]
    print(2 ** n % 1000000)


@app.command("seto")
def seto(file: str):
    """Introduction to Set Operations"""
    n, s1, s2 = Parser(file).lines()
    n = int(n)
    s1 = builtins.eval(s1)
    s2 = builtins.eval(s2)
    s3 = set(range(1, n + 1))
    print(*[s1 | s2, s1 & s2, s1 - s2, s2 - s1, s3 - s1, s3 - s2], sep="\n")


@app.command("spec")
def spec(file: str):
    """Inferring Protein from Spectrum"""
    weights = [float(x) for x in Parser(file).lines()]
    diff = [j - i for i, j in zip(weights[:-1], weights[1:])]
    print("".join([mass.match_mass(x) for x in diff]))


@app.command("pdst")
def pdst(file: str):
    """Creating a Distance Matrix"""
    for r in aln.pdst(Parser(file).seqs()):
        print(*[round(x, 3) for x in r])


@app.command("aspc")
def aspc(file: str):
    """Introduction to Alternative Splicing"""
    n, k = Parser(file).ints()
    print(sum([math.comb(n, x) for x in range(k, n + 1)]) % 1000000)


@app.command("cat")
def cat(file: str):
    """Catalan Numbers and RNA Secondary Structures"""
    print(com.cat(*Parser(file).seqs()))


@app.command("inod")
def inod(file: str):
    """Counting Phylogenetic Ancestors"""
    print(Parser(file).ints()[0] - 2)


@app.command("mmch")
def mmch(file: str):
    """Maximum Matchings and RNA Secondary Structures"""
    print(com.mmch(Parser(file).seqs()[0]))


@app.command("afrq")
def afrq(file: str):
    """Counting Disease Carriers"""
    b = pr.afrq(Parser(file).floats())
    print(*[round(x, 3) for x in b])


@app.command("conv")
def conv(file: str):
    """Comparing Spectra with the Spectral Convolution"""
    l1, l2 = Parser(file).lines()
    s1 = list(map(float, l1.split()))
    s2 = list(map(float, l2.split()))
    res = mass.conv(s1, s2)
    print(res[1], res[0], sep="\n")


@app.command("ebin")
def ebin(file: str):
    """Wright-Fisher's Expected Behavior"""
    l1, l2 = Parser(file).lines()
    n = int(l1)
    s2 = list(map(float, l2.split()))
    print(*[round(x * n, 3) for x in s2])


@app.command("scsp")
def scsp(file: str):
    """Interleaving Two Motifs"""
    seqs = Parser(file).lines()
    out = aln.scsp(seqs[0], seqs[1])
    print(out["ss"], sep="\n")


@app.command("edit")
def edit(file: str):
    """Edit Distance"""
    seqs = Parser(file).seqs()
    print(aln.edit(seqs[0], seqs[1]))


@app.command("edta")
def edta(file: str):
    """Edit Distance Alignment"""
    seqs = Parser(file).seqs()
    out = aln.edta(seqs[0], seqs[1])
    print(out["dist"], out["a1"], out["a2"], sep="\n")


@app.command("ctea")
def ctea(file: str):
    """Counting Optimal Alignments"""
    seqs = Parser(file).seqs()
    print(aln.ctea(seqs[0], seqs[1]))


@app.command("glob")
def glob(file: str):
    """Global Alignment with Scoring Matrix"""
    seqs = Parser(file).seqs()
    print(aln.glob(seqs[0], seqs[1]))


@app.command("gcon")
def gcon(file: str):
    """Global Alignment with Constant Gap Penalty"""
    seqs = Parser(file).seqs()
    print(aln.gcon(seqs[0], seqs[1]))


@app.command("gaff")
def gaff(file: str):
    """Global Alignment with Scoring Matrix and Affine Gap Penalty"""
    seqs = Parser(file).seqs()
    res = aln.gaff(seqs[0], seqs[1], -11, -1)
    print(res["dist"])
    print(res["a1"])
    print(res["a2"])


@app.command("eval")
def eval(file: str):
    """Expected Number of Restriction Sites"""
    n, s, a = Parser(file).lines()
    n = int(n)
    a = map(float, a.split())
    print(*[round(pr.eval(n, s, x), 3) for x in a])


@app.command("indc")
def indc(file: str):
    """Independent Segregation of Chromosomes"""
    print(*pr.indc(*Parser(file).ints()))


@app.command("foun")
def foun(file: str):
    """The Founder Effect and Genetic Drift"""
    l1, l2 = Parser(file).lines()
    n, m = [int(x) for x in l1.split()]
    a = [int(x) for x in l2.split()]
    for x in pr.foun(n, m, a):
        print(*[round(f, 8) for f in x])


@app.command("sexl")
def sexl(file: str):
    """Sex-Linked Inheritance"""
    arr = Parser(file).floats()
    print(*[round(2 * v * (1 - v), 3) for v in arr])


@app.command("wfmd")
def wfmd(file: str):
    """The Wright-Fisher Model of Genetic Drift"""
    print(pr.wfmd(*Parser(file).ints()))


@app.command("lcsq")
def lcsq(file: str):
    """Locating Motifs Despite Introns"""
    print(aln.lcsq(*Parser(file).seqs()))


@app.command("motz")
def motz(file: str):
    """Motzkin Numbers and RNA Secondary Structures"""
    print(com.motz(Parser(file).seqs()[0]))


@app.command("pmch")
def pmch(file: str):
    """Perfect Matchings and RNA Secondary Structures"""
    print(com.pmch(Parser(file).seqs()[0]))


@app.command("lexv")
def lexv(file: str):
    """Ordering Strings of Varying Length Lexicographically"""
    l1, l2 = Parser(file).lines()
    print(*ros.lexv(l1.split(), int(l2)), sep="\n")


# @app.command("trie")
# def trie(file: str):
#     """Introduction to Pattern Matching"""
#     seqs = Parser(file).lines()
#     for line in graph.trie(seqs):
#         print(*line)


def main():
    app()
