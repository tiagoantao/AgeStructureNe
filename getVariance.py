from __future__ import print_function

import os
import sys

import myUtils

if len(sys.argv) not in [2]:
    print("python %s varConfFile" % sys.argv[0])
    sys.exit(-1)

varConfFile = sys.argv[1]

myDir = varConfFile.split("/")[0]

N0, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = myUtils.getVarConf(varConfFile)

models = list(N0.keys())
models.sort()

outdir = 'output/variance'

try:
    os.mkdir(outdir)
except OSError:
    pass

for model in models:
    Ns = N0[model]
    Ns.sort()
    try:
        for N in Ns:
            w = open('%s/%s-%d.txt' % (outdir, model, N), 'w')
            cfg = myUtils.getConfig("%s/%d%s.conf" % (myDir, N, model))
            startGen = cfg.gens - numGens
            for rep in range(reps):
                for rec in myUtils.getDemo(myDir, model, N, rep):
                    cycle = rec['cycle']
                    pyramid = rec['pyramid']
                    if cycle >= startGen:
                        w.write('%d\t' % sum(pyramid))
                w.write('\n')
            w.write('\n')
            for rep in range(reps):
                for rec in myUtils.getVk(myDir, model, N, rep):
                    cycle = rec['cycle']
                    nb = rec['nb']
                    if cycle >= startGen:
                        w.write('%d\t' % nb)
                w.write('\n')
            w.close()
    except IOError:
        os.remove('%s/%s-%d.txt' % (outdir, model, N))
