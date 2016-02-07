import os
import sys

import myUtils

if len(sys.argv) not in [2]:
    print "python {0!s} varConfFile".format(sys.argv[0])
    sys.exit(-1)

varConfFile = sys.argv[1]

myDir = varConfFile.split("/")[0]

N0, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = myUtils.getVarConf(varConfFile)

models = N0.keys()
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
            w = open('{0!s}/{1!s}-{2:d}.txt'.format(outdir, model, N), 'w')
            cfg = myUtils.getConfig("{0!s}/{1:d}{2!s}.conf".format(myDir, N, model))
            startGen = cfg.gens - numGens
            for rep in range(reps):
                for rec in myUtils.getDemo(myDir, model, N, rep):
                    cycle = rec['cycle']
                    pyramid = rec['pyramid']
                    if cycle >= startGen:
                        w.write('{0:d}\t'.format(sum(pyramid)))
                w.write('\n')
            w.write('\n')
            for rep in range(reps):
                for rec in myUtils.getVk(myDir, model, N, rep):
                    cycle = rec['cycle']
                    nb = rec['nb']
                    if cycle >= startGen:
                        w.write('{0:d}\t'.format(nb))
                w.write('\n')
            w.close()
    except IOError:
        os.remove('{0!s}/{1!s}-{2:d}.txt'.format(outdir, model, N))
