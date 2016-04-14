from __future__ import print_function
import sys

import myUtils

if len(sys.argv) not in [3]:
    print("Do not forget to run plotHz")
    print(("python %s varConfFile prefix" % sys.argv[0]))
    sys.exit(-1)

confFile = sys.argv[1]
pref = sys.argv[2]


def doCase(w, age, indivs, loci, isSNP, isRel, startGens, numGens,
           pcrit, case):
    for model in models:
        Ns = N0[model]
        Ns.sort()
        for N in Ns:
            print("%s\t%s\t%s\t%d" % (age, case, model, N), end=' ', file=w)
            if type(startGens) == int:
                startGen = startGens  # hack
            else:
                try:
                    startGen = startGens[str(N) + model]
                except KeyError:
                    print("err", N, model, file=sys.stderr)
                    print("\t", end=' ', file=w)
                    print(file=w)
                    continue
            if pcrit is None:
                name = ""
            else:
                name = "{pcrit}-".format(pcrit=pcrit)
            name += str(N) + model + age + str(indivs) + str(loci)
            name += "-snp-" if isSNP else "-"
            name += "rel-" if isRel else ""
            ldne = []
            ldneCI = []
            r2 = []
            sr2 = []
            j = []
            ssize = []
            try:
                for rep in range(reps):
                    mName = name + str(rep)
                    f = open(dataDir + "/ldout/" + mName)
                    l = f.readline().rstrip()
                    if l == "":
                        raise IOError("Not complete")
                    l = l.replace("inf", "1000000")
                    le = eval(l)
                    ldne.extend(le[startGen:startGen + numGens])
                    ldne = [x if x > 0 else 1000000 for x in ldne]
                    l = f.readline().rstrip()
                    if l == "":
                        raise IOError("Not complete")
                    l = l.replace("inf", "1000000")
                    le = eval(l)
                    ldneCI.extend(le[startGen:startGen + numGens])
                    f.close()

                    f = open(dataDir + '/ldout/' + mName + '.r2')
                    l = f.readline()
                    le = eval(l)
                    r2.extend(le[startGen:startGen + numGens])
                    l = f.readline()
                    le = eval(l)
                    sr2.extend(le[startGen:startGen + numGens])
                    l = f.readline()
                    le = eval(l)
                    j.extend(le[startGen:startGen + numGens])
                    l = f.readline()
                    le = eval(l)
                    ssize.extend(le[startGen:startGen + numGens])
                    f.close()

                if len(ldne) > 0:
                    print("\t" + str(ldne), end=' ', file=w)
                    print("\t" + str(ldneCI), end=' ', file=w)
                    print("\t" + str(r2), end=' ', file=w)
                    print("\t" + str(sr2), end=' ', file=w)
                    print("\t" + str(j), end=' ', file=w)
                    print("\t" + str(ssize), end=' ', file=w)
            except IOError:
                print("err", mName, model, indivs, loci, rep, name,
                      file=sys.stderr)
            print(file=w)

N0, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = \
    myUtils.getVarConf(confFile)

models = list(N0.keys())
models.sort()


def doHz(w, startGens):
    for model in models:
        Ns = N0[model]
        Ns.sort()
        for N in Ns:
            cfg = myUtils.getConfig(dataDir + "/" + str(N) + model + ".conf")
            if len(startGens) == 0:  # Use last
                startGen = cfg.gens - numGens - 1
            else:
                startGen = startGens

    ages = list(sampCohort.keys())
    for age in ages:
        for indivs, loci in sampSize:
            case = '\t%d\t%d\tMSAT' % (indivs, loci)
            doCase(w, age, indivs, loci, False, False, startGen, numGens,
                   None, case)
            case = '\t%d\t%d\tMSAT-rel' % (indivs, loci)
            doCase(w, age, indivs, loci, False, True, startGen, numGens,
                   None, case)
            for pcrit in [0.021, 0.035, 0.05, 0.1]:
                case = '%s\t%d\t%d\tMSAT' % ("{pcrit}".format(pcrit=pcrit),
                                             indivs, loci)
                doCase(w, age, indivs, loci, False, False, startGen, numGens,
                       pcrit, case)
        for indivs, loci in sampSNP:
            case = '\t%d\t%d\tSNP' % (indivs, loci)
            doCase(w, age, indivs, loci, True, False, startGen, numGens,
                   None, case)
            for pcrit in [0.021, 0.035, 0.05, 0.1]:
                case = '%s\t%d\t%d\tSNP' % ("{pcrit}".format(pcrit=pcrit),
                                            indivs, loci)
                doCase(w, age, indivs, loci, True, False, startGen, numGens,
                       pcrit, case)
            case = '\t%d\t%d\tSNP-rel' % (indivs, loci)
            doCase(w, age, indivs, loci, True, True, startGen, numGens,
                   None, case)

w = open("output/%s.txt" % pref, "w")
doHz(w, {})
w.close()
allCuts = {}
f = open("output/hz-cut.txt")
for l in f:
    toks = l.rstrip().split("\t")
    model, cut, gen = toks[0], float(toks[1]), int(toks[2])
    # Taking the cut from SNP
    # Should get a different one for MSATs
    model, marker = model.split("-")
    if marker == "0":  # msat
        continue
    modelCuts = allCuts.setdefault(cut, {})
    modelCuts[model] = gen
f.close()
for cut, modelCuts in list(allCuts.items()):
    print(cut)
    print((list(modelCuts.keys())))
    w = open("output/%s-%d.txt" % (pref, cut * 100), "w")
    doHz(w, modelCuts)
    w.close()
