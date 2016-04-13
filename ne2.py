from __future__ import print_function
import os
import sys

from genomics.popgen.ne2.Controller import NeEstimator2Controller
from genomics.popgen import ne2

if len(sys.argv) != 2:
    print("%s <thres>")
    exit(-1)

thres = float(sys.argv[1])

ne2c = NeEstimator2Controller()

fname = 'ld' + str(os.getpid()) + '.ne2'
out = open(fname, 'w')
l = sys.stdin.readline()
cnt = 0
while l != "":
    out.write(l)
    l = sys.stdin.readline()
    cnt += 1
if cnt == 0:  # Sample size above indivs
    sys.exit(0)
out.close()
ne2c.run_neestimator2('.', fname, '.', fname + '.out', crits=[thres])
ldout = open(fname + '.out')
ldres = ne2.parse(ldout)

mNes = []
mOr2s = []
mSmpr2s = []
mNesPow = []
mNesCI = []
mIndep = []
mHMean = []
for fcases in ldres.ld:
    case = fcases[0]
    ne = case['EstNe']
    or2 = case['OvRSquare']
    sr2 = case['ExpRSquareSample']
    indep = case['IndepComp']
    hmean = case['HMean']
    ne05, ne975 = tuple(case['ParaNe'])
    mNes.append(ne)
    mOr2s.append(or2)
    mSmpr2s.append(sr2)
    mNesCI.append((ne975, ne05))
    mIndep.append(indep)
    mHMean.append(hmean)
ldout.close()

os.remove(fname)
os.remove(fname + ".out")
print(str(mNes))
print(str(mNesCI))
print(str(mOr2s), file=sys.stderr)
print(str(mSmpr2s), file=sys.stderr)
print(str(mIndep), file=sys.stderr)
print(str(mHMean), file=sys.stderr)
sys.stdout.flush()
