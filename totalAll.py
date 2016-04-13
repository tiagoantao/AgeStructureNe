from __future__ import print_function
import multiprocessing
import sys
import os

import myUtils
import Executor

if len(sys.argv) not in [2, 3]:
    print("python %s varConfFile dogen" % sys.argv[0])
    sys.exit(-1)

varConfFile = sys.argv[1]
if len(sys.argv) == 3:
    doGen = True
else:
    doGen = False

myDir = varConfFile.split("/")[0]

lexec = Executor.Local(-multiprocessing.cpu_count())
#lexec = Executor.Local(-10)

N0, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = myUtils.getVarConf(varConfFile)


def nongen(model, N, ageM, ageF, reps, startGen):
    age = ageM  # XXX simplification
    for rep in range(reps):
        lexec.err = "stderr"
        lexec.submit("bash", "nongenAll.sh %d %d %d%s %s %d %d %d" %
                     (rep, rep + 1, N, model, myDir, startGen,
                      startGen + numGens, age))


def gen(model, N, ageM, ageF, reps):
    print(1, reps)
    os.chdir(myDir)
    cond = "(sex==1 and age>%d) or (sex==2 and age>%d)" % (ageM, ageF)
    for rep in range(reps):
        if os.path.isfile("bothTop.sh"):
            lexec.submit("bash", 'bothTop.sh %d %d %d%s "%s"' %
                         (rep, rep + 1, N, model, cond))
        else:
            lexec.submit("bash", 'topGo.sh %d %d %d%s "%s"' %
                         (rep, rep + 1, N, model, cond))
            print(("bash", 'topGo.sh %d %d %d%s "%s"' %
                  (rep, rep + 1, N, model, cond)))
    os.chdir("..")
    #bash topGo.sh $REPS $REPE ${a[$i]}  ${age[$i]} ;

models = list(N0.keys())
models.sort()

for model in models:
    Ns = N0[model]
    Ns.sort()
    for N in Ns:
        print(model, N)
        cfg = myUtils.getConfig("%s/%d%s.conf" % (myDir, N, model))
        startGen = cfg.gens - numGens
        ageM, ageF = myUtils.getAgeFecund(cfg)
        if doGen:
            gen(model, N, ageM, ageF, reps)
        else:
            nongen(model, N, ageM, ageF, reps, startGen)
lexec.wait(True)
