from __future__ import division

import math

from scipy.stats import chi2

dataDir = "trout"
nindivs = [15, 25, 50, 100]
nlocis = [15, 100]
nsnps = [100, 200, 400]
#N0s = [90, 180, 361, 722]  # , 1805]
Nbs = {('bulltrout', 86): 25, ('bulltrout', 174): 50,
       ('bulltrout', 353): 100, ('bulltrout', 713): 200,
       ('bullpred', 189): 50, ('bullpred', 381): 100, ('bullpred', 765): 200,
       ('bullt2', 606): 10, ('bullt2', 1214): 20,
       ('bullt2', 1822): 30, ('bullt2', 2433): 40,
       ('bullt2', 3040): 50, ('bullt2', 4565): 75, ('bullt2', 6090): 100,
       ('shepard', 513): 50, ('shepard', 1030): 100,
       ('fraley', 637): 50, ('fraley', 1278): 100,
       ('wfrog', 597): 100, ('wfrog', 297): 50, ('wfrog', 148): 25,
       ('mosquito', 102): 25, ('mosquito', 167): 40,
       ('mosquito', 210): 50, ('mosquito', 428): 100,
       ('synseaweed', 81): 20, ('synseaweed', 164): 40,
       ('synseaweed', 207): 50, ('synseaweed', 422): 100
       #('sagegrouse', 28): 100, ('sagegrouse', 14): 50,
       }
NbNames = [("bulltrout", "BT-Std"), ("bullpred", "BT-Pred"),
           ("bullt2", "BT-Long"),
           ("shepard", "WCT-S"), ("fraley", "WCT-F"), ('restricted', "Restr"),
           ('mosquito', 'Mosq'), ('grizzly', 'Griz'), ('wfrog', 'WFrog'),
           ('sagegrouse', 'SGrouse'), ('seaweed', 'Sweed'),
           ('synseaweed', 'SynSweed')]
#Nes = {90: 32.8, 180: 65.2, 361: 130.3, 722: 259.8}
cohorts = ["All",  "Newb", "c2c", "c3c"]
cuts = [0.45, 0.4, 0.35, 0.3, 0.25]

realNbs = {}

nb_corrs = {"BT-Std": {86: 0.7993630573, 174: 0.7908082409,
                         353: 0.7849293564, 713: 0.7798129384},
            "BT-Long": {606: 0.78125, 1214: 0.78125, 1822: 0.78125,
                        2433: 0.7797270955, 3040: 0.7800312012,
                        4565: 0.7796257796, 6090: 0.7788161994},
            "BT-Pred": {189: 0.6622516556, 381: 0.6602902375,
                        765: 0.6585446164},
            "WCT-S": {513: 0.702247191, 1030: 0.700280112},
            "WCT-F": {637: 0.7132667618, 1278: 0.7117437722},
            'WFrog': {148: 0.5966587112, 297: 0.596889952, 597: 0.5976119403},
            'Mosq': {102: 0.2835990888, 167: 0.2835990888,
                     210: 0.2773763202, 428: 0.2735229759},
            'SynSweed': {81: 1.3673469388, 164: 1.3468013468, 207: 1.336,
                         422: 1.309},
            #'SGrouse': 1.694,
            }

log_a = 0.15458222
log_b = -0.671991958
log_c = 0.799127

corrs = {}
for sp, nbne in nb_corrs.items():
    corrs[sp] = {}
    for N1, mcorr in nbne.items():
        corrs[sp][N1] = 1 / (1.26 - 0.323 * mcorr)

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


def correct_ci(N0, bname, nindivs, vals, ci, r2=None, fixed=None,
               sr2=None, j=None, jcorr=1):
    cvals = []
    cci = []
    diffs = []
    if fixed is None:
        my_corr = corrs[bname][N0]
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


def correct_logquad(N0, bname, nindivs, vals, ci, abc,
                    r2=None, sr2=None, j=None, jcorr=None, cut=0.025):
    a, b, c = abc
    cvals = []
    cci = []

    for i in range(len(ci)):
        val = vals[i]
        if r2 is None:
            bot, top = ci[i]
        else:
            bot, top = patch_ci(nindivs, r2[i], sr2[i], j[i] * jcorr, cut)
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


def get_corrs(N0, bname, nindivs, vals, ci, r2, sr2, j):
    return [("None", (vals, ci)),
            ("NbNe", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                sr2=sr2, j=j)),
            ("NbNe0.9", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                   sr2=sr2, j=j, jcorr=0.9)),
            ("NbNe0.5", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                   sr2=sr2, j=j, jcorr=0.5)),
            ("NbNe0.4", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                   sr2=sr2, j=j, jcorr=0.4)),
            ("NbNe0.3", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                   sr2=sr2, j=j, jcorr=0.3)),
            ("NbNe0.2", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                   sr2=sr2, j=j, jcorr=0.2)),
            ("NbNe0.1", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                   sr2=sr2, j=j, jcorr=0.1)),
            ("NbNe0.05", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                    sr2=sr2, j=j, jcorr=0.05)),
            ("NbNe0.01", correct_ci(N0, bname, nindivs, vals, ci, r2=r2,
                                    sr2=sr2, j=j, jcorr=0.01)),
            #("0.9", correct_ci(bname, nindivs, vals, ci, fixed=0.9)),
            #("Int0.9", correct_ci(bname, nindivs, vals, ci, fixed=-0.9)),
            #("Log0.01", correct_logquad(N0, bname, nindivs, vals, ci,
            #                            [log_a, log_b,
            #                             log_c], r2=r2, sr2=sr2,
            #                            j=j, jcorr=0.01)),
            #("Log0.05", correct_logquad(N0, bname, nindivs, vals, ci,
            #                            [log_a, log_b,
            #                             log_c], r2=r2, sr2=sr2,
            #                            j=j, jcorr=0.05)),
            #("Log0.1", correct_logquad(N0, bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.1)),
            #("Log0.2", correct_logquad(N0, bname, nindivs, vals, ci,
            #                           [log_a, log_b,
            #                            log_c], r2=r2, sr2=sr2,
            #                           j=j, jcorr=0.2)),
            #("LogQuad", correct_logquad(N0, bname, nindivs, vals, ci,
            #                            [log_a, log_b,
            #                             log_c])),
            ]


def get_bname(model):
    for name, bname in NbNames:
        if name == model:
            return bname
