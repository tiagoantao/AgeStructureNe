import bz2
import sys

import scipy
from scipy import stats

from myUtils import getVarConf, getConfig

cf = sys.argv[1]

N0s, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = getVarConf(cf)

models = N0s.keys()
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

        age = (gen-1) - born[pi]
        cntPossParents +=1

        if gender[pi]==1:
            fec = cfg.fecundityMale
            cntPossDads +=1
            ofsCntDads[pi] = 0
        else:
            fec = cfg.fecundityFemale
            cntPossMoms +=1
            ofsCntMoms[pi] = 0
        ofsCnt[pi] = 0
    for parent in parents:
        #print prev.index(parent)
        ofsCnt[parent] +=1
        if gender[parent] == 1:
            ofsCntDads[parent] += 1
        else:
            ofsCntMoms[parent] += 1
    kbar = 2.0*cntOfs/cntPossParents
    a = scipy.array(ofsCnt.values())
    vk = a.var()

    kbarm = 1.0*cntOfs/cntPossDads
    a = scipy.array(ofsCntDads.values())
    vkm = a.var()

    kbarf = 1.0*cntOfs/cntPossMoms
    a = scipy.array(ofsCntMoms.values())
    vkf = a.var()
    nb = (kbar*cntPossParents-2)/(kbar-1+vk/kbar)
    return cntOfs, cntPossParents, len(set(parents)), kbar, vk, kbarm, vkm, kbarf, vkf, nb


def doModel(model, N0, rep, cfg):
    f = bz2.BZ2File("data/%s/%d%s%d.sim.bz2" %(dataDir, N0, model, rep))
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
        if gen>currGen:
            if currGen>=myStart:
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
        if age==1:
            gender[id] = sex
            parents.extend([father, mother])
            born[id] = gen
    f.close()
    return res

def dumpMeans(allRes):
    keys = allRes[0].keys()
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
        print key,
        for i in range(numVals):
            v = acu[i]/len(allRes)
            print v,
            vals[i].append(v)
        print
    return vals

for model in models:
    for N0 in N0s[model]:
        allRes = []
        print model, N0
        print "gen N1 Nall Npar kbar vk kbarm vkm kbarf vkf nb"
        cfg = getConfig(dataDir + "/" + str(N0) + model + ".conf")
        startGen = cfg.gens - numGens - 1
        try:
            for rep in range(reps):
                repCase = doModel(model, N0, rep, cfg)
                print >>sys.stderr, model, N0, rep,
                myGens = repCase.keys()
                myGens.sort()
                print >>sys.stderr, myGens[0],
                for myGen in myGens:
                    print >>sys.stderr, repCase[myGen][-1],
                print >>sys.stderr
                #print >>sys.stderr, model, N0, rep, repCase
                miniRep = {}
                for gen, vals in repCase.items():
                    if gen >= startGen:
                        miniRep[gen] = vals
                allRes.append(miniRep)
        except IOError:
            print "NoData", rep
            continue
        vals = dumpMeans(allRes)
        print "",
        for i in range(len(vals)-1):
            print scipy.mean(vals[i]),
        try:
            hm = stats.hmean(vals[-1])
            print hm
        except ValueError:
            print "Negs",
            print stats.hmean(filter(lambda x: x > 0,vals[-1]))
