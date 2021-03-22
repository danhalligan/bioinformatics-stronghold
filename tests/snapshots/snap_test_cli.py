# -*- coding: utf-8 -*-
# snapshottest: v1 - https://goo.gl/zC4yUc
from __future__ import unicode_literals

from snapshottest import Snapshot


snapshots = Snapshot()

snapshots['test_cli_function[aspc] 1'] = '''42
'''

snapshots['test_cli_function[cat] 1'] = '''2
'''

snapshots['test_cli_function[cons] 1'] = '''ATGCAACT
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6
'''

snapshots['test_cli_function[corr] 1'] = '''TTCAT->TTGAT
GAGGA->GATGA
TTTCC->TTTCA
'''

snapshots['test_cli_function[dna] 1'] = '''20 12 17 21
'''

snapshots['test_cli_function[fib] 1'] = '''19
'''

snapshots['test_cli_function[fibd] 1'] = '''4
'''

snapshots['test_cli_function[gc] 1'] = '''Rosalind_0808
60.91954
'''

snapshots['test_cli_function[grph] 1'] = '''Rosalind_0498 Rosalind_2391
Rosalind_0498 Rosalind_0442
Rosalind_2391 Rosalind_2323
'''

snapshots['test_cli_function[hamm] 1'] = '''7
'''

snapshots['test_cli_function[iev] 1'] = '''3.5
'''

snapshots['test_cli_function[inod] 1'] = '''2
'''

snapshots['test_cli_function[iprb] 1'] = '''0.7833333333333333
'''

snapshots['test_cli_function[kmer] 1'] = '''4 1 4 3 0 1 1 5 1 3 1 2 2 1 2 0 1 1 3 1 2 1 3 1 1 1 1 2 2 5 1 3 0 2 2 1 1 1 1 3 1 0 0 1 5 5 1 5 0 2 0 2 1 2 1 1 1 2 0 1 0 0 1 1 3 2 1 0 3 2 3 0 0 2 0 8 0 0 1 0 2 1 3 0 0 0 1 4 3 2 1 1 3 1 2 1 3 1 2 1 2 1 1 1 2 3 2 1 1 0 1 1 3 2 1 2 6 2 1 1 1 2 3 3 3 2 3 0 3 2 1 1 0 0 1 4 3 0 1 5 0 2 0 1 2 1 3 0 1 2 2 1 1 0 3 0 0 4 5 0 3 0 2 1 1 3 0 3 2 2 1 1 0 2 1 0 2 2 1 2 0 2 2 5 2 2 1 1 2 1 2 2 2 2 1 1 3 4 0 2 1 1 0 1 2 2 1 1 1 5 2 0 3 2 1 1 2 2 3 0 3 0 1 3 1 2 3 0 2 1 2 2 1 2 3 0 1 2 3 1 1 3 1 0 1 1 3 0 2 1 2 2 0 2 1 1
'''

snapshots['test_cli_function[kmp] 1'] = '''0 0 0 1 2 0 0 0 0 0 0 1 2 1 2 3 4 5 3 0 0
'''

snapshots['test_cli_function[lcsm] 1'] = '''AC
'''

snapshots['test_cli_function[lexf] 1'] = '''AA
AC
AG
AT
CA
CC
CG
CT
GA
GC
GG
GT
TA
TC
TG
TT
'''

snapshots['test_cli_function[lgis] 1'] = '''1 2 3
5 4 2
'''

snapshots['test_cli_function[lia] 1'] = '''0.68359375
'''

snapshots['test_cli_function[long] 1'] = '''ATTAGACCTGCCGGAATAC
'''

snapshots['test_cli_function[mmch] 1'] = '''6
'''

snapshots['test_cli_function[mprt] 1'] = '''B5ZC00
85 118 142 306 395
P07204_TRBM_HUMAN
47 115 116 382 409
P20840_SAG1_YEAST
79 109 135 248 306 348 364 402 485 501 614
'''

snapshots['test_cli_function[mrna] 1'] = '''12
'''

snapshots['test_cli_function[orf] 1'] = '''M
MGMTPRLGLESLLE
MLLGSFRLIPKETLIQVAGSSPCNLS
MTPRLGLESLLE
'''

snapshots['test_cli_function[pdst] 1'] = '''0.0 0.4 0.1 0.1
0.4 0.0 0.4 0.3
0.1 0.4 0.0 0.2
0.1 0.3 0.2 0.0
'''

snapshots['test_cli_function[perm] 1'] = '''6
1 2 3
1 3 2
2 1 3
2 3 1
3 1 2
3 2 1
'''

snapshots['test_cli_function[pper] 1'] = '''51200
'''

snapshots['test_cli_function[prob] 1'] = '''-5.737 -5.217 -5.263 -5.36 -5.958 -6.628 -7.009
'''

snapshots['test_cli_function[prot] 1'] = '''MAMAPRTEINSTRING
'''

snapshots['test_cli_function[prtm] 1'] = '''821.3919199999999
'''

snapshots['test_cli_function[rear] 1'] = '''9 4 5 7 0
'''

snapshots['test_cli_function[revc] 1'] = '''ACCGGGTTTT
'''

snapshots['test_cli_function[revp] 1'] = '''4 6
5 4
6 6
7 4
17 4
18 4
20 6
21 4
'''

snapshots['test_cli_function[rna] 1'] = '''GAUGGAACUUGACUACGUAAAUU
'''

snapshots['test_cli_function[rstr] 1'] = '''0.6885137241525929
'''

snapshots['test_cli_function[sign] 1'] = '''8
-1 -2
-1 2
1 -2
1 2
-2 -1
-2 1
2 -1
2 1
'''

snapshots['test_cli_function[spec] 1'] = '''WMQS
'''

snapshots['test_cli_function[splc] 1'] = '''MVYIADKQHVASREAYGHMFKVCA
'''

snapshots['test_cli_function[sseq] 1'] = '''3 4 5
'''

snapshots['test_cli_function[sset] 1'] = '''8
'''

snapshots['test_cli_function[subs] 1'] = '''2 4 10
'''

snapshots['test_cli_function[tran] 1'] = '''1.2142857142857142
'''

snapshots['test_cli_function[tree] 1'] = '''3
'''
