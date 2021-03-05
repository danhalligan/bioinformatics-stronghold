from rosalind import __version__
from rosalind.rosalind import *
from rosalind.helpers import *


def testdata(fun=None):
    if fun is None:
        fun = inspect.stack()[1].function
    return Parser("tests/data/" + fun + ".txt")


def test_version():
    assert __version__ == "0.1.0"


def test_dna():
    dat = testdata().line()
    assert count_nucleotides(dat) == [20, 12, 17, 21]


def test_rna():
    inp = "GATGGAACTTGACTACGTAAATT"
    out = "GAUGGAACUUGACUACGUAAAUU"
    assert dna2rna(inp) == out


def test_revcomp():
    inp = "AAAACCCGGT"
    out = "ACCGGGTTTT"
    assert revcomp(inp) == out


def test_rabbits():
    assert rabbits(5, 3) == 19


def test_gc():
    dat = testdata().fastas()
    assert max_gc(dat) == {"name": "Rosalind_0808", "value": 60.919540}


def hamming_distance():
    dat = ["GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT"]
    assert hamming_distance(dat) == 7


def test_mendel1():
    assert round(mendel1(2, 2, 2), 5) == 0.78333


def test_translate():
    seq = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
    assert translate(seq) == "MAMAPRTEINSTRING"


def test_find_motif():
    s1, s2 = ["GATATATGCATATACTT", "ATAT"]
    assert list(find_motif(s1, s2)) == [2, 4, 10]


def test_profile_matrix():
    exp = [
        [5, 1, 0, 0, 5, 5, 0, 0],
        [0, 0, 1, 4, 2, 0, 6, 1],
        [1, 1, 6, 3, 0, 1, 0, 0],
        [1, 5, 0, 0, 0, 1, 1, 6],
    ]

    dat = testdata("test_cons").fastas()
    seqs = [x.seq for x in dat]
    out = profile_matrix(seqs)
    assert all([all(out[i] == exp[i]) for i in range(3)])


def test_consensus_sequence():
    dat = testdata("test_cons").fastas()
    seqs = [x.seq for x in dat]
    out = consensus_sequence(profile_matrix(seqs))
    assert out == "ATGCAACT"


def test_mortal_rabbits():
    assert mortal_rabbits(6, 3) == 4


def test_expected_offspring():
    assert expected_offspring([1, 0, 0, 1, 0, 1]) == 3.5


def test_find_shared_motif():
    seqs = [x.seq for x in testdata().fastas()]
    assert find_shared_motif(seqs) == "AC"


def test_mendel2():
    assert round(mendel2(2, 1), 3) == 0.684


def test_find_protein_motif():
    seq = "LDNFSDPLIDCKNCKANYSTDL"
    assert find_protein_motif(seq) == [3, 17]


def test_find_orfs():
    seq = str(testdata().fastas()[0].seq)
    exp = {"MLLGSFRLIPKETLIQVAGSSPCNLS", "M", "MGMTPRLGLESLLE", "MTPRLGLESLLE"}
    assert set(find_orfs(seq)) == exp


def test_reverse_pallindromes():
    exp = [[4, 6], [5, 4], [6, 6], [7, 4], [17, 4], [18, 4], [20, 6], [21, 4]]
    out = list(reverse_pallindromes("TCAATGCATGCGGGTCTATATGCAT"))
    assert sorted(out) == exp


def test_count_rnas():
    assert count_rnas("MA") == 12
