import os
import sys

from Bio.PopGen.GenePop.Controller import GenePopController

from igrat.genetics.popgen.ldne.Controller import LDNeController
from igrat.genetics.popgen import ldne

if len(sys.argv) != 2:
    print "%s <thres>"
    exit(-1)

thres = sys.argv[1]  # float but we really prefer string here

os.chdir(thres)

ldnec = LDNeController()
gp_ctrl = GenePopController()

fname = "ld" + str(os.getpid()) + ".ldne"
out = open(fname, "w")


l = sys.stdin.readline()
cnt = 0
while l != "":
    out.write(l)
    l = sys.stdin.readline()
    cnt += 1
if cnt == 0:  # Sample size above indivs
    sys.exit(0)
out.close()
ldnec.run_ldne(fname, fname + ".out")
ldout = open(fname + '.out')
ldres = ldne.RecordParser().parse(ldout)
mNes = []
mOr2s = []
mNesPow = []
mNesCI = []
for id, fcases in ldres.populations:
    hm, ic, or2, er2, (ne, (ne95, ne05), (j95, j05)) = fcases[0]
    mNes.append(ne)
    mOr2s.append(or2)
    mNesCI.append((ne95, ne05))
ldout.close()

popi, loci = gp_ctrl.calc_allele_genotype_freqs(fname)
poli = []
for pop in popi:
    name, locs = pop
    genos_len = [len(x[0]) for x in locs.values()]
    poli.append(len(genos_len) - genos_len.count(1))
os.remove(fname)
os.remove(fname + ".out")
print mNes
print mNesCI
print >>sys.stderr, mOr2s
print >>sys.stderr, poli
sys.stdout.flush()

os.chdir("..")
