import sys
import os

import myUtils
import Executor

if len(sys.argv) not in [2]:
    print "python %s varConfFile" % sys.argv[0]
    sys.exit(-1)

varConfFile = sys.argv[1]

myDir = varConfFile.split("/")[0]

N0, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = myUtils.getVarConf(varConfFile)



models = N0.keys()
models.sort()

for model in models:
    Ns = N0[model]
    Ns.sort()
    for N in Ns:
        cfg = myUtils.getConfig("%s/%d%s.conf" % (myDir, N, model))
        startGen = cfg.gens - numGens
