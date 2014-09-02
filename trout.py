from __future__ import division

import math

from scipy.stats import chi2

dataDir = "trout"
nindivs = [15, 25, 50, 100]
nlocis = [15, 100]
nsnps = [100, 200, 400]
N0s = [90, 180, 361, 722]  # , 1805]
Nbs = {('bulltrout', 90): 25, ('bulltrout', 180): 50, ('bulltrout', 361): 100,
       ('bulltrout', 722): 200,  # 1805: 500,
       ('restricted', 90): 25, ('restricted', 180): 50,
       ('restricted', 361): 100, ('restricted', 722): 200,
       ('bullpred', 193): 50, ('bullpred', 387): 100, ('bullpred', 775): 200,
       ('bullt2', 305): 5, ('bullt2', 610): 10, ('bullt2', 915): 15,
       ('bullt2', 1220): 20,
       ('bullt2', 1525): 25, ('bullt2', 1830): 30, ('bullt2', 2440): 40,
       ('bullt2', 3050): 50, ('bullt2', 4575): 75, ('bullt2', 6100): 100,
       ('btrout', 1619): 50, ('btrout', 6476): 200,
       ('shepard', 518): 50, ('shepard', 1036): 100,
       ('fraley', 641): 50, ('fraley', 1282): 100,
       ('wfrog', 600): 100, ('wfrog', 300): 50, ('wfrog', 150): 25,
       ('mosquito', 110): 25, ('mosquito', 176): 40,
       ('mosquito', 218): 50, ('mosquito', 436): 100,
       ('seaweed', 70): 100, ('seaweed', 35): 50,
       ('synseaweed', 89): 20, ('synseaweed', 178): 40,
       ('synseaweed', 219): 50,
       ('sagegrouse', 28): 100, ('sagegrouse', 14): 50,
       ('lake', 18): 50, ('lake', 36): 100
       }
NbNames = [("bulltrout", "BuTrout"), ("bullpred", "BuPred"), ("bullt2", "BuLong"),
           ("shepard", "WCT-S"), ("fraley", "WCT-F"), ('restricted', "Restr"),
           ('mosquito', 'Mosq'), ('grizzly', 'Griz'), ('wfrog', 'WFrog'),
           ('sagegrouse', 'SGrouse'), ('seaweed', 'Sweed'),
           ('synseaweed', 'SynSweed')]
Nes = {90: 32.3, 180: 64.7, 361: 129.4, 722: 258.8}
cohorts = ["All",  "Newb", "c2c", "c3c"]
cuts = [0.45, 0.4, 0.35, 0.3, 0.25]

realNbs = {}

nb_corrs = {"BuTrout": 0.77,
            "Restricted": 0.77,
            "BuLong": 0.78,
            "BuPred": 0.65,
            "WCT-S": 0.69,
            "WCT-F": 0.71,
            'WFrog': 0.599,
            'Griz': 1.153,
            'Mosq': 0.267,
            'Sweed': 1.261,
            'SynSweed': 1.261,
            'SGrouse': 1.694}

log_a = 0.15458222
log_b = -0.671991958
log_c = 0.799127

corrs = {}
for sp, nbne in nb_corrs.items():
    corrs[sp] = 1 / (1.26 - 0.323 * nbne)

pcrits = [None, 0.021, 0.035, 0.05, 0.1]


def load_file(pref, cut=None, mydir='.'):
    if cut is None:
        f = open("%s/output/%s.txt" % (mydir, pref))
    else:
        f = open("%s/output/%s-%d.txt" % (mydir, pref, cut))
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

        if len(toks) > 7:
            modeln0 = curr_case.setdefault(case_id, {})
            ldne = eval(toks[7])
            ci = eval(toks[8])
            r2 = eval(toks[9])
            sr2 = eval(toks[10])
            j = eval(toks[11])
            ssize = eval(toks[12])
            modeln0[(pcrit, nindivs, nloci, my_type)] = (ldne, ci, r2, sr2,
                                                         j, ssize)
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


def get_ldne(nindivs, sr2, r2):
    if nindivs < 30:
        #r2l = 0.0018 + 0.907 / nindivs + 4.44 / nindivs ** 2
        myr2 = r2 - sr2
        try:
            ldne = (.308 + math.sqrt(.308 ** 2 - 2.08 * myr2)) / (2 * myr2)
        except ValueError:
            return None
    else:
        #r2l = 1 / nindivs + 3.19 / nindivs ** 2
        myr2 = r2 - sr2
        try:
            ldne = (1 / 3 + math.sqrt(1 / 9 - 2.76 * myr2)) / (2 * myr2)
        except ValueError:
            return None
    return ldne if ldne > 0 else 100000


def patch_ci(nindivs, r2, sr2, j, cut=0.025):
    b = chi2.isf(cut, j)
    t = chi2.isf(1 - cut, j)
    r2b, r2t = j * r2 / b, j * r2 / t
    return get_ldne(nindivs, sr2, r2b), get_ldne(nindivs, sr2, r2t)


def correct_ci(bname, nindivs, vals, ci, r2=None, fixed=None,
               sr2=None, j=None, jcorr=1):
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
    if r2 is not None:
        # We are going to patch the CIs with r2
        ci = []
        for i in range(len(r2)):
            ci.append(patch_ci(nindivs, r2[i], sr2[i], j[i] * jcorr))
    for i in range(len(ci)):
        bot, top = ci[i]
        if my_corr < 0:
            if bot is None:
                cbot = None
            else:
                cbot = (bot - diffs[i]) * (- my_corr)
            if top is None:
                ctop = None
            else:
                ctop = (top - diffs[i]) * (1 + (1 + my_corr))
            cci.append((cbot, ctop))
        else:
            bcorr = None if bot is None else bot * my_corr
            tcorr = None if top is None else top * my_corr
            cci.append((bcorr, tcorr))
    return cvals, cci


def correct_logquad(bname, nindivs, vals, ci, abc,
                    r2=None, sr2=None, j=None, jcorr=None, cut=0.025):
    a, b, c = abc
    cvals = []
    cci = []

    for i in range(len(ci)):
        val = vals[i]
        if r2 is None:
            bot, top = ci[i]
        else:
            bot, top = patch_ci(r2[i], sr2[i], j[i] * jcorr, cut)
        if val >= 100000:
            # Not correcting
            cvals.append(val)
            cci.append((bot, top))
            continue
        lval = math.log10(val)
        corr = a * lval * lval + b * lval + c
        corr = 1 - corr
        diff = val - val * corr
        cvals.append(val * corr)
        bdiff = bot - diff if bot is not None else None
        tdiff = top - diff if top is not None else None
        cci.append((bdiff, tdiff))
    return cvals, cci


def get_corrs(bname, nindivs, vals, ci, r2, sr2, j):
    return [("None", (vals, ci)),
            #("NbNe", correct_ci(bname, nindivs, vals, ci)),
            ("NbNe", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                  sr2=sr2, j=j)),
            ("NbNe0.01", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                     sr2=sr2, j=j, jcorr=0.01)),
            ("NbNe0.05", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                     sr2=sr2, j=j, jcorr=0.05)),
            ("NbNe0.1", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                    sr2=sr2, j=j, jcorr=0.1)),
            ("NbNe0.2", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                    sr2=sr2, j=j, jcorr=0.2)),
            ("NbNe0.3", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                    sr2=sr2, j=j, jcorr=0.3)),
            ("NbNe0.4", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                    sr2=sr2, j=j, jcorr=0.4)),
            ("NbNe0.5", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                    sr2=sr2, j=j, jcorr=0.5)),
            ("NbNe0.9", correct_ci(bname, nindivs, vals, ci, r2=r2,
                                    sr2=sr2, j=j, jcorr=0.9)),
            #("0.9", correct_ci(bname, nindivs, vals, ci, fixed=0.9)),
            #("Int0.9", correct_ci(bname, nindivs, vals, ci, fixed=-0.9)),
            #("LogQuad", correct_logquad(bname, nindivs, vals, ci,
            #                            [-0.17599607, 0.75649721, 0.07641839]))]
            #("LogQuad2", correct_logquad(bname, nindivs, vals, ci,
            #                             [0.144624, -0.654859, 0.8009])),
            #("Log0.001", correct_logquad(bname, nindivs, vals, ci,
            #                            [log_a, log_b,
            #                             log_c], r2=r2, sr2=sr2,
            #                            j=j, jcorr=0.001)),
            #("Log0.01", correct_logquad(bname, nindivs, vals, ci,
            #                            [log_a, log_b,
            #                             log_c], r2=r2, sr2=sr2,
            #                            j=j, jcorr=0.01)),
            #("Log0.05", correct_logquad(bname, nindivs, vals, ci,
            #                            [log_a, log_b,
            #                             log_c], r2=r2, sr2=sr2,
            #                            j=j, jcorr=0.05)),
            #("Log0.1", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.1)),
            #("Log0.2", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.2)),
            #("Log0.2-90", correct_logquad(bname, nindivs, vals, ci,
            #                              [log_a, log_b,
            #                               log_c], r2=r2, sr2=sr2,
            #                              j=j, jcorr=0.2, cut=0.05)),
            #("Log0.3", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.3)),
            #("Log0.4", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.4)),
            #("Log0.5", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.5)),
            #("Log0.6", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.6)),
            #("Log0.7", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.7)),
            #("Log0.8", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.8)),
            #("Log0.9", correct_logquad(bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.9)),
            ("LogQuad", correct_logquad(bname, nindivs, vals, ci,
                                        [log_a, log_b,
                                         log_c])),
            ]


def get_bname(model):
    for name, bname in NbNames:
        if name == model:
            return bname
