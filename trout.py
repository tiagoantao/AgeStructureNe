from __future__ import division

import math

dataDir = "trout"
nindivs = [15, 25, 50, 100]
nlocis = [15]
nsnps = [100, 200, 400]
N0s = [90, 180, 361, 722]  # , 1805]
Nbs = {('bulltrout', 90): 25, ('bulltrout', 180): 50, ('bulltrout', 361): 100,
       ('bulltrout', 722): 200,  # 1805: 500,
       ('bullpred', 193): 50, ('bullpred', 387): 100, ('bullpred', 775): 200,
       ('bullt2', 305): 5, ('bullt2', 610): 10, ('bullt2', 915): 15,
       ('bullt2', 1220): 20,
       ('bullt2', 1525): 25, ('bullt2', 1830): 30, ('bullt2', 2440): 40,
       ('bullt2', 3050): 50, ('bullt2', 4575): 75, ('bullt2', 6100): 100,
       # "btrout-1619": 50, "btrout-6476": 200,
       ('shepard', 518): 50, ('shepard', 1036): 100,
       ('fraley', 641): 50, ('fraley', 1282): 100}
#       "lake-18": 50, "lake-72": 200}
NbNames = [("bulltrout", "BuTrout"), ("bullpred", "BuPred"), ("bullt2", "BuLong"),
           ("shepard", "WCT-S"), ("fraley", "WCT-F")]
Nes = {90: 32.3, 180: 64.7, 361: 129.4, 722: 258.8}
cohorts = ["All",  "Newb", "c2c", "c3c"]
cuts = [0.45, 0.4, 0.35, 0.3, 0.25]

realNbs = {}

corrs = {"BuTrout": 0.77,
         "BuLong": 0.78,
         "BuPred": 0.65,
         "WCT-S": 0.69,
         "WCT-F": 0.71}

pcrits = [None, 0.021, 0.035, 0.05, 0.1]


def load_file(pref, cut=None):
    if cut is None:
        f = open("output/%s.txt" % pref)
    else:
        f = open("output/%s-%d.txt" % (pref, cut))
    case = {}
    f.readline()
    f.readline()
    for l in f:
        toks = l.rstrip().split("\t")
        cohort = toks[0]
        pcrit = float(toks[1]) if toks[1] != '' else None
        nindivs = int(toks[2])
        nloci = int(toks[3])
        my_type = toks[4]
        model = toks[5]
        N0 = int(toks[6])
        curr_case = case.setdefault(cohort, {})
        case_id = model, N0

        def get_vals(tok):
            ts = tok.rstrip()
            if ts == "":
                return [], []
            else:
                try:
                    point, ci = ts.split("!")
                    return eval(point), eval(ci)
                except ValueError:
                    return [], []
        if len(toks) > 7:
            modeln0 = curr_case.setdefault(case_id, {})
            modeln0[(pcrit, nindivs, nloci, my_type)] = eval(toks[7]), eval(toks[8])
    return case


def load_nb(pref):
    f = open("output/" + pref + "-nb.txt")
    for l in f:
        toks = l.rstrip().split(" ")
        model = toks[0]
        N0 = int(toks[1])
        rep = int(toks[2])
        startGen = int(toks[3])
        vals = [float(x) for x in toks[4:]]
        realNbs.setdefault((model, N0), {})
        realNbs[model, N0][rep] = startGen, vals
    f.close()


def correct_ci(bname, vals, ci, fixed=None):
    cvals = []
    cci = []
    diffs = []
    if fixed is None:
        my_corr = corrs[bname]
    else:
        my_corr = fixed
    for v in vals:
        cvals.append(v * abs(my_corr))
        if my_corr < 0:
            diffs.append(v * (1 + my_corr))
    for i in range(len(ci)):
        bot, top = ci[i]
        if my_corr < 0:
            cbot = (bot - diffs[i]) * (- my_corr)
            ctop = (top - diffs[i]) * (1 + (1 + my_corr))
            cci.append((cbot, ctop))
        else:
            cci.append((bot * my_corr, top * my_corr))
    return cvals, cci


def correct_logquad(bname, vals, ci, abc):
    a, b, c = abc
    cvals = []
    cci = []

    for i in range(len(ci)):
        val = vals[i]
        bot, top = ci[i]
        lval = math.log10(val)
        corr = a * lval * lval + b * lval + c
        corr = 1 - corr
        #print val, val*corr, lval, corr, "X"
        diff = val - val * corr
        cvals.append(val * corr)
        cci.append((bot - diff, top - diff))
    return cvals, cci


def get_corrs(bname, vals, ci):
    return [("None", (vals, ci)),
            ("Nb/Ne", correct_ci(bname, vals, ci)),
            ("0.9", correct_ci(bname, vals, ci, 0.9)),
            ("Int0.9", correct_ci(bname, vals, ci, -0.9)),
            #("LogQuad", correct_logquad(bname, vals, ci,
            #                            [-0.17599607, 0.75649721, 0.07641839]))]
            ("LogQuad", correct_logquad(bname, vals, ci,
                                        [0.15458222, -0.671991958, 0.799127]))]


def get_bname(model):
    for name, bname in NbNames:
        if name == model:
            return bname
