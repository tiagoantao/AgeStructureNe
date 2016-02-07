import sys

import myUtils

if len(sys.argv) not in [3]:
    print("Do not forget to run plotHz")
    print("python {0!s} varConfFile prefix".format(sys.argv[0]))
    sys.exit(-1)

confFile = sys.argv[1]
pref = sys.argv[2]


def doCase(w, age, indivs, loci, isSNP, isRel, startGens, numGens,
           pcrit, case):
    for model in models:
        Ns = N0[model]
        Ns.sort()
        for N in Ns:
            print >>w, "{0!s}\t{1!s}\t{2!s}\t{3:d}".format(age, case, model, N),
            if type(startGens) == int:
                startGen = startGens  # hack
            else:
                try:
                    startGen = startGens[str(N) + model]
                except KeyError:
                    print >>sys.stderr, "err", N, model
                    print >>w, "\t",
                    print >>w
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
                    ldne = map(lambda x: x if x > 0 else 1000000, ldne)
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
                    print >>w, "\t" + str(ldne),
                    print >>w, "\t" + str(ldneCI),
                    print >>w, "\t" + str(r2),
                    print >>w, "\t" + str(sr2),
                    print >>w, "\t" + str(j),
                    print >>w, "\t" + str(ssize),
            except IOError:
                print >>sys.stderr, "err", mName, model, indivs, loci, rep, name
            print >>w

N0, sampCohort, sampSize, sampSNP, numGens, reps, dataDir = myUtils.getVarConf(confFile)

models = N0.keys()
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

    ages = sampCohort.keys()
    for age in ages:
        for indivs, loci in sampSize:
            case = '\t{0:d}\t{1:d}\tMSAT'.format(indivs, loci)
            doCase(w, age, indivs, loci, False, False, startGen, numGens,
                   None, case)
            case = '\t{0:d}\t{1:d}\tMSAT-rel'.format(indivs, loci)
            doCase(w, age, indivs, loci, False, True, startGen, numGens,
                   None, case)
            for pcrit in [0.021, 0.035, 0.05, 0.1]:
                case = '{0!s}\t{1:d}\t{2:d}\tMSAT'.format("{pcrit}".format(pcrit=pcrit),
                                             indivs, loci)
                doCase(w, age, indivs, loci, False, False, startGen, numGens,
                       pcrit, case)
        for indivs, loci in sampSNP:
            case = '\t{0:d}\t{1:d}\tSNP'.format(indivs, loci)
            doCase(w, age, indivs, loci, True, False, startGen, numGens,
                   None, case)
            for pcrit in [0.021, 0.035, 0.05, 0.1]:
                case = '{0!s}\t{1:d}\t{2:d}\tSNP'.format("{pcrit}".format(pcrit=pcrit),
                                            indivs, loci)
                doCase(w, age, indivs, loci, True, False, startGen, numGens,
                       pcrit, case)
            case = '\t{0:d}\t{1:d}\tSNP-rel'.format(indivs, loci)
            doCase(w, age, indivs, loci, True, True, startGen, numGens,
                   None, case)

w = open("output/{0!s}.txt".format(pref), "w")
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
for cut, modelCuts in allCuts.items():
    print(cut)
    print(modelCuts.keys())
    w = open("output/{0!s}-{1:d}.txt".format(pref, cut * 100), "w")
    doHz(w, modelCuts)
    w.close()
