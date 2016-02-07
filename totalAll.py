import multiprocessing
import sys
import os

import myUtils
import Executor

if len(sys.argv) not in [2, 3]:
    print "python {0!s} varConfFile dogen".format(sys.argv[0])
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
        lexec.submit("bash", "nongenAll.sh {0:d} {1:d} {2:d}{3!s} {4!s} {5:d} {6:d} {7:d}".format(rep, rep + 1, N, model, myDir, startGen,
                      startGen + numGens, age))


def gen(model, N, ageM, ageF, reps):
    print 1, reps
    os.chdir(myDir)
    cond = "(sex==1 and age>{0:d}) or (sex==2 and age>{1:d})".format(ageM, ageF)
    for rep in range(reps):
        if os.path.isfile("bothTop.sh"):
            lexec.submit("bash", 'bothTop.sh {0:d} {1:d} {2:d}{3!s} "{4!s}"'.format(rep, rep + 1, N, model, cond))
        else:
            lexec.submit("bash", 'topGo.sh {0:d} {1:d} {2:d}{3!s} "{4!s}"'.format(rep, rep + 1, N, model, cond))
            print("bash", 'topGo.sh {0:d} {1:d} {2:d}{3!s} "{4!s}"'.format(rep, rep + 1, N, model, cond))
    os.chdir("..")
    #bash topGo.sh $REPS $REPE ${a[$i]}  ${age[$i]} ;

models = N0.keys()
models.sort()

for model in models:
    Ns = N0[model]
    Ns.sort()
    for N in Ns:
        print model, N
        cfg = myUtils.getConfig("{0!s}/{1:d}{2!s}.conf".format(myDir, N, model))
        startGen = cfg.gens - numGens
        ageM, ageF = myUtils.getAgeFecund(cfg)
        if doGen:
            gen(model, N, ageM, ageF, reps)
        else:
            nongen(model, N, ageM, ageF, reps, startGen)
lexec.wait(True)
