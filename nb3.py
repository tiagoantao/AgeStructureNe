from __future__ import print_function
import bz2
import sys

import scipy
from scipy import stats

from myUtils import getVarConf, getConfig

cf = sys.argv[1]

N0s, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = getVarConf(cf)

models = list(N0s.keys())
models.sort()


def dumpOut(gen, born, curr, prev, parents, gender, cfg):
    cntPossParents = 0
    cntPossDads = 0
    cntPossMoms = 0
    cntOfs = 0
    ofsCnt = {}
    ofsCntDads = {}
    ofsCntMoms = {}
    for ci in curr:
        age = gen - born[ci]
        if age == 1:
            cntOfs += 1
    for pi in prev:

        age = (gen - 1) - born[pi]
        cntPossParents += 1

        if gender[pi] == 1:
            # fec = cfg.fecundityMale
            cntPossDads += 1
            ofsCntDads[pi] = 0
        else:
            # fec = cfg.fecundityFemale
            cntPossMoms += 1
            ofsCntMoms[pi] = 0
        ofsCnt[pi] = 0
    for parent in parents:
        # print prev.index(parent)
        ofsCnt[parent] += 1
        if gender[parent] == 1:
            ofsCntDads[parent] += 1
        else:
            ofsCntMoms[parent] += 1
    kbar = 2.0 * cntOfs / cntPossParents
    a = scipy.array(list(ofsCnt.values()))
    vk = a.var()

    kbarm = 1.0 * cntOfs / cntPossDads
    a = scipy.array(list(ofsCntDads.values()))
    vkm = a.var()

    kbarf = 1.0 * cntOfs / cntPossMoms
    a = scipy.array(list(ofsCntMoms.values()))
    vkf = a.var()
    nb = (kbar * cntPossParents - 2) / (kbar - 1 + vk / kbar)
    return cntOfs, cntPossParents, \
        len(set(parents)), kbar, vk, kbarm, vkm, kbarf, vkf, nb


def doModel(model, N0, rep, cfg):
    f = bz2.BZ2File("data/%s/%d%s%d.sim.bz2" % (dataDir, N0, model, rep))
    myStart = 50
    currGen = 0
    gender = {}
    curr = []
    prev = []
    parents = []
    born = {}
    res = {}
    for l in f:
        toks = l.rstrip().split(" ")
        gen = int(toks[0])
        if gen > currGen:
            if currGen >= myStart:
                res[gen] = dumpOut(gen, born, curr, prev, parents, gender, cfg)
            prev = curr
            curr = []
            parents = []
            currGen = gen
        id = int(float(toks[2]))
        age = int(float(toks[6]))
        sex = int(toks[3])
        father = int(float(toks[4]))
        mother = int(float(toks[5]))
        curr.append(id)
        if age == 1:
            gender[id] = sex
            parents.extend([father, mother])
            born[id] = gen
    f.close()
    return res


def dumpMeans(allRes):
    keys = list(allRes[0].keys())
    numVals = len(allRes[0][keys[0]])
    keys.sort()
    vals = []
    for i in range(numVals):
        vals.append([])
    for key in keys:
        acu = []
        for i in range(numVals):
            acu.append(0.0)
        for i in range(len(allRes)):
            myRes = allRes[i][key]
            for j in range(numVals):
                acu[j] += myRes[j]
        print(key, end=' ')
        for i in range(numVals):
            v = acu[i] / len(allRes)
            print(v, end=' ')
            vals[i].append(v)
        print()
    return vals

for model in models:
    for N0 in N0s[model]:
        allRes = []
        print(model, N0)
        print("gen N1 Nall Npar kbar vk kbarm vkm kbarf vkf nb")
        cfg = getConfig(dataDir + "/" + str(N0) + model + ".conf")
        startGen = cfg.gens - numGens - 1
        try:
            for rep in range(reps):
                repCase = doModel(model, N0, rep, cfg)
                print(model, N0, rep, end=' ', file=sys.stderr)
                myGens = list(repCase.keys())
                myGens.sort()
                print(myGens[0], end=' ', file=sys.stderr)
                for myGen in myGens:
                    print(repCase[myGen][-1], end=' ', file=sys.stderr)
                print(file=sys.stderr)
                # print >>sys.stderr, model, N0, rep, repCase
                miniRep = {}
                for gen, vals in list(repCase.items()):
                    if gen >= startGen:
                        miniRep[gen] = vals
                allRes.append(miniRep)
        except IOError:
            print("NoData", rep)
            continue
        vals = dumpMeans(allRes)
        print("", end=' ')
        for i in range(len(vals) - 1):
            print(scipy.mean(vals[i]), end=' ')
        try:
            hm = stats.hmean(vals[-1])
            print(hm)
        except ValueError:
            print("Negs", end=' ')
            print(stats.hmean([x for x in vals[-1] if x > 0]))
