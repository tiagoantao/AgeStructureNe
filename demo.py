from __future__ import division, print_function
import sys

from scipy import stats
import numpy

import myUtils

if len(sys.argv) != 2:
    print("python %s varConfFile" % (sys.argv[0], ))
    sys.exit(-1)

confFile = sys.argv[1]

N0, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = \
  myUtils.getVarConf(confFile)

models = list(N0.keys())
models.sort()


def doModels(fun):
    for model in models:
        Ns = N0[model]
        Ns.sort()
        for N in Ns:
            cfg = myUtils.getConfig(dataDir + "/" + str(N) + model + ".conf")
            startGen = cfg.gens - numGens - 1
            # sys.stdout.write("startGen: %d\n" % startGen)
            sys.stdout.write(model + "\t")
            sys.stdout.write(str(N) + "\t")
            try:
                ret = fun(model, N, startGen)
                sys.stdout.write("\t".join([str(x) for x in ret]))
            except IOError:
                sys.stdout.write("not done")
            sys.stdout.write("\n")


def cleanSpaces(ls, minVal):
    vals = []
    for l in ls:
        toks = [float(y) for y in
                [x for x in l.rstrip().split(" ") if x != ""]]
        if toks[0] < minVal:
            continue
        vals.append(toks)
    return vals


def acu(fnamer, harit, startGen):
    ret = []
    acu = {}
    cnt = 0.0
    for rep in range(reps):
        fname = fnamer % rep
        f = open(fname)
        ls = f.readlines()
        f.close()
        vals = cleanSpaces(ls, startGen)
        for val in vals:
            cnt += 1
            for i in range(1, len(val)):
                if i in harit:
                    acu[i-1] = acu.get(i-1, 0) + 1.0/val[i]
                else:
                    acu[i-1] = acu.get(i-1, 0) + val[i]
    ks = list(acu.keys())
    ks.sort()
    for k in ks:
        if k+1 in harit:
            ret.append(cnt/acu[k])
        else:
            ret.append(acu[k]/cnt)
    return ret


def pyr(model, N, startGen):
    fname = "data/%s/%d%s%%d.demo" % (dataDir, N, model)
    myPyr = acu(fname, [], startGen)[1:]
    return [sum(myPyr)] + myPyr


def vk(model, N, startGen):
    if doMature:
        fname = "data/%s/%d%s%%d.vk.mature" % (dataDir, N, model)
    else:
        fname = "data/%s/%d%s%%d.vk" % (dataDir, N, model)
    return acu(fname, [3, 4], startGen)


def doOfs():
    for model in models:
        Ns = N0[model]
        Ns.sort()
        for N in Ns:
            print("%s\t%d" % (model, N))
            cfg = myUtils.getConfig(dataDir + "/" + str(N) + model + ".conf")
            startGen = cfg.gens - numGens - 1
            nowMale = {}
            totMale = {}
            nowFemale = {}
            totFemale = {}
            maxV = 0
            cnt = 0.0
            for rep in range(reps):
                fname = "data/%s/%d%s%d.ofs" % (dataDir, N, model, rep)
                f = open(fname)
                ls = f.readlines()
                f.close()
                vals = cleanSpaces(ls, startGen)
                for val in vals:
                    cnt += 1
                    sex = val[3]
                    topOfs = val[4]
                    nowOfs = val[5]
                    if int(topOfs) > maxV:
                        maxV = int(topOfs)
                    if sex == 1:
                        nowMale[nowOfs] = nowMale.get(nowOfs, 0) + 1
                        totMale[topOfs] = totMale.get(topOfs, 0) + 1
                    else:
                        nowFemale[nowOfs] = nowFemale.get(nowOfs, 0) + 1
                        totFemale[topOfs] = totFemale.get(topOfs, 0) + 1
            for i in range(maxV):
                if i > 20 and i < maxV-1:
                    continue
                print("\t\t%d\t%f\t%f\t%f\t%f" % (i, nowMale.get(i, 0.0)/cnt,
                      nowFemale.get(i, 0)/cnt, totMale.get(i, 0)/cnt,
                      totFemale.get(i, 0)/cnt))

print()
print("All individuals")
print("model\tN1\tk\tvk\t (k*popSizes[gen] - 2) / (k-1+ (Vk/k))\tHillNe")
doMature = False
doModels(vk)

print()
print("MATURE")
print("model\tN1\tk\tvk\t (k*popSizes[gen] - 2) / (k-1+ (Vk/k))\tHillNe")
try:
    doMature = True
    doModels(vk)
except IOError:
    pass

try:
    print("Model\tN1\tcnt\tnowMaleOfs\tnowFemaleOfs\ttotMaleOfs\ttotFemaleOfs")
    doOfs()
except IOError:
    pass

print("model\tN1\tNc\tN1\tN2\t...")
doModels(pyr)
