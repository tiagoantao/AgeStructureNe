from __future__ import division

import functools
import math
import os
import shutil
import sys

import numpy as np
from scipy.stats import hmean

import matplotlib
matplotlib.use('AGG')
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

import seaborn as sns
import matplotlib.pyplot as plt

from trout import Nbs, nindivs, nlocis, cohorts, NbNames, cuts, load_file
from trout import pcrits, dataDir, realNbs, get_corrs, correct_ci, nsnps
from trout import get_bname
Nes = None  # Wrong for now

pref = sys.argv[1]


def do_nb(case, cohort, N0, nsnps, pref, corr_name):
    fig, ax = plt.subplots()
    model = 'bulltrout'
    bname = get_bname(model)
    fig.suptitle("Nb: %d (N1: %d) - cohort %s - %s" % (Nbs[model, N0],
                                                       N0, cohort, corr_name))
    box_vals = []
    tops = []
    bottoms = []
    labels = []
    hmeans = []
    for nindiv in nindivs:
        for nloci in nlocis:
            try:
                vals, ci, r2, sr2, j, ssize = \
                    case[cohort][(model, N0)][(None, nindiv, nloci, "MSAT")]
                for cname, corrections in get_corrs(N0, bname, nindiv, vals,
                                                    ci, r2, sr2, j):
                    if cname != corr_name:
                        continue
                    cvals, cci = corrections
                    vals = cvals
                    ci = cci
                    break
                if len(ci) > 0:
                    bottom, top = zip(*ci)
                    tops.append(np.percentile([x if x is not None else
                                               100000 for x in top], 90))
                    bottoms.append(np.percentile([x if x is not None
                                                  else 0.1 for x
                                                  in bottom], 10))
                else:
                    tops.append(None)
                    bottoms.append(None)
            except KeyError:
                vals = []
                tops.append(None)
                bottoms.append(None)
            box_vals.append(vals)
            hmeans.append(hmean([x for x in vals if x > 0]))
            labels.append("%dMS" % nloci)
        for nsnp in nsnps:
            try:
                vals, ci, r2, sr2, j, ssize = case[cohort][(model, N0)][(None, nindiv, nsnp, "SNP")]
                for cname, corrections in get_corrs(N0, bname, nindiv, vals,
                                                    ci, r2, sr2, j):
                    if cname != corr_name:
                        continue
                    cvals, cci = corrections
                    vals = cvals
                    ci = cci
                    break
                if len(ci) > 0:
                    bottom, top = zip(*ci)
                    tops.append(np.percentile([x if x is not None else
                                               100000 for x in top], 90))
                    bottoms.append(np.percentile([x if x is not None
                                                  else 0.1 for x in
                                                  bottom], 10))
                else:
                    tops.append(None)
                    bottoms.append(None)

            except KeyError:
                vals = []
                tops.append(None)
                bottoms.append(None)
            box_vals.append(vals)
            hmeans.append(hmean([x for x in vals if x > 0]))
            #Check above
            labels.append("%d" % nsnp)
        pos = len(labels)
        ax.axvline(pos + 0.5, color="k", lw=0.2)
        ax.text(pos - 0.5, 0, "%d Indivs" % nindiv,
                ha="center", va="bottom", size="small",
                rotation="horizontal")
    ax.set_ylabel("$\hat{N}_{e}$")
    ax.set_ylim(0, Nbs[(model, N0)] * 3)
    ax.axhline(Nbs[(model, N0)], color="k", lw=0.3)
    sns.boxplot(box_vals, notch=0, sym="")
    ax.set_xticks(1 + np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation="vertical")
    #ax.plot([1 + x for x in range(len(tops))], tops, "r.", ms=20)
    #ax.plot([1 + x for x in range(len(bottoms))], bottoms, "r.", ms=20)
    #ax.plot([1 + x for x in range(len(hmeans))], hmeans, "r.", ms=20)
    #fig.savefig("output/%s%s-%s-%d.png" % (pref, cohort, corr_name, N0))
    return fig


def do_cohort(case, model, N0, nindiv, corr_name):
    last = 0.5
    fig, ax = plt.subplots()
    fig.suptitle("Nb: %d (N1: %d) - different cohorts - 100 SNPs -%s" %
                (Nbs[(model, N0)], N0, corr_name))
    box_vals = []
    labels = []
    tops = []
    bottoms = []
    hmeans = []
    bname = get_bname(model)

    for cohort in cohorts:
        print(cohort, model, N0)
        vals, ci, r2, sr2, j, ssize = \
            case[cohort][(model, N0)][(None, nindiv, 100, "SNP")]
        for cname, corrections in get_corrs(N0, bname, nindiv, vals,
                                            ci, r2, sr2, j):
            if cname != corr_name:
                continue
            cvals, cci = corrections
            vals = cvals
            ci = cci
            break
        box_vals.append(vals)
        hmeans.append(hmean(vals))
        bottom, top = zip(*ci)
        top = [100000 if x is None else x for x in top]
        bottom = [100000 if x is None else x for x in bottom]
        tops.append(np.percentile(top, 90))
        bottoms.append(np.percentile(bottom, 10))
        labels.append("%s" % cohort)
        if cohort == cohorts[-1]:
            pos = len(labels) + 0.5
            ax.axvline(pos, color="k", lw=0.2)
            ax.text(last + (pos - last) / 2, 0, "%d Indivs" % nindiv,
                    ha="center", va="bottom", size="small",
                    rotation="horizontal")
            last = pos
    ax.set_ylim(0, Nbs[(model, N0)] * 3)
    ax.set_ylabel("$\hat{N}_{e}$")
    ax.axhline(Nbs[(model, N0)], color="k", lw=0.3)
    sns.boxplot(box_vals, notch=0, sym="")
    ax.set_xticks(1 + np.arange(len(labels)))
    ax.set_xticklabels(labels)
    ax.plot([1 + x for x in range(len(tops))], tops, "rx")
    ax.plot([1 + x for x in range(len(bottoms))], bottoms, "rx")
    ax.plot([1 + x for x in range(len(hmeans))], hmeans, "k+")
    #fig.savefig("output/cohort-%s-%s-%d.png" % (model, corr_name, N0))
    return fig


def do_rel(case, model, N0, nindiv, corr_name):
    fig, ax = plt.subplots()
    cohort = "Newb"
    fig.suptitle("Nb: %d (N1: %d) - 20pc related individuals - cohort %s - %s"
                 % (Nbs[(model, N0)], N0, cohort, corr_name))
    box_vals = []
    labels = []
    bname = get_bname(model)

    # vals, ci, r2, sr2, j, ssize = case[cohort][N0][(None, nindiv, 15, "MSAT-rel")]
    # box_vals.append(vals)
    # labels.append("MSAT-15-rel")
    # vals, ci, r2, sr2, j, ssize = case[cohort][N0][(None, nindiv, 15, "MSAT")]
    # box_vals.append(vals)
    # labels.append("MSAT-15")
    vals, ci, r2, sr2, j, ssize = case[cohort][(model, N0)][(None, nindiv, 100, "SNP-rel")]
    for cname, corrections in get_corrs(N0, bname, nindiv, vals, ci, r2, sr2, j):
        if cname != corr_name:
            continue
        cvals, cci = corrections
        vals = cvals
        ci = cci
        break
    box_vals.append(vals)
    labels.append("SNP-100-rel")
    vals, ci, r2, sr2, j, ssize = case[cohort][(model, N0)][(None, nindiv, 100, "SNP")]
    for cname, corrections in get_corrs(N0, bname, nindiv, vals, ci, r2, sr2, j):
        if cname != corr_name:
            continue
        cvals, cci = corrections
        vals = cvals
        ci = cci
        break
    box_vals.append(vals)
    labels.append("SNP-100")
    sns.boxplot(box_vals, notch=0, sym="")
    ax.set_ylim(0, Nbs[(model, N0)] * 3)
    ax.set_ylabel("$\hat{N}_{e}$")
    ax.axhline(Nbs[(model, N0)], color="k", lw=0.3)
    ax.set_xticks(1 + np.arange(len(labels)))
    ax.set_xticklabels(labels)
    #fig.savefig("output/rel-%s-%d.png" % (model, N0))
    return fig


def do_nb_comp(case):
    fig, ax = plt.subplots()
    box_vals = []
    labels = []
    n0s = Nbs.keys()
    n0s.sort()
    fig.suptitle("Nb comparison - 100 SNPs - 50 indivs")
    for model, n0 in n0s:
        if type(n0) != int:
            continue
        nb = Nbs[(model, n0)]
        vals, ci, r2, sr2, j, ssize = \
            case["Newb"][(model, n0)][(None, 50, 100, "SNP")]
        labels.append(str(nb))
        box_vals.append(vals)
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylabel("$\hat{N}_{e}$")
    sns.boxplot(box_vals, notch=0, sym="")
    #fig.savefig("output/nb-comp.png")


def do_lt_comp(case, nb, strat, corr_name):
    fig, ax = plt.subplots()
    nindiv = 50
    tops = []
    box_vals = []
    bottoms = []
    hmeans = []
    labels = []
    n0s = list(Nbs.keys())
    n0s.sort()
    ax.set_title("Nb: %d - %s sampled - 100 SNPs - 50 indivs - %s corr" % (
        nb, strat, corr_name))
    for pref, name in NbNames:
        if name == "Restr":
            continue
        for k, nb2 in Nbs.items():
            model, n0 = k
            if nb2 != nb or n0 < 50:
                continue
            elif pref != model:
                continue
            vals, ci, r2, sr2, j, ssize = \
                case[strat][(model, n0)][(None, 50, 100, "SNP")]
            for cname, corrections in get_corrs(n0, name, nindiv, vals,
                                                ci, r2, sr2, j):
                if cname != corr_name:
                    continue
                cvals, cci = corrections
                vals = cvals
                ci = cci
                break
            hmeans.append(hmean(vals))
            try:
                top, bottom = zip(*ci)
            except ValueError:
                top, bottom = [], []
            top = [100000 if x is None else x for x in top]
            bottom = [100000 if x is None else x for x in bottom]
            labels.append(name)
            try:
                tops.append(np.percentile(top, 90))
                box_vals.append(vals)
                bottoms.append(np.percentile(bottom, 10))
            except ValueError:
                tops.append(None)
                box_vals.append([])
                bottoms.append(None)
    sns.violinplot(box_vals, notch=0, sym="", ax=ax, alpha=0.9)
    #ax.plot([1 + x for x in range(len(tops))], tops, 'r.', ms=20)
    #ax.plot([1 + x for x in range(len(bottoms))], bottoms, "r.", ms=20)
    #ax.plot([1 + x for x in range(len(hmeans))], hmeans, "k.", ms=20)
    ymin, ymax = ax.get_ylim()
    plt.ylabel("$\hat{N}_{e}$")
    plt.ylim(ymin, min([ymax, 3 * nb]))
    plt.axhline(nb, color="k", lw=0.3)
    ax.set_xticks(range(1, 1 + len(labels)))
    ax.set_xticklabels(labels)
    #plt.savefig("output/lt-comp-%s-%s-%d.png" % (corr_name, strat, nb))
    return fig


def do_bias(case):
    n0s = Nbs.keys()
    n0s.sort()
    sampling = "Newb"
    table = [["Nb", "Sampling", "Model", "N1", "NeEst",
              "NeEst/Nb", "Above", "Below"]]
    for model, n0 in n0s:
        nb = Nbs[(model, n0)]
        vals, ci, r2, sr2, j, ssize = \
            case[sampling][(model, n0)][(None, 50, 100, "SNP")]
        for p, bname in NbNames:
            if p == model:
                break
        cnt = 0
        below = 0
        above = 0
        for i in range(len(vals)):
            mi, ma = ci[i]
            if nb < mi:
                below += 1
            if nb > ma:
                above += 1
            cnt += 1
        if cnt > 0:
            table.append([nb, sampling, bname, n0, np.median(vals),
                          np.median(vals) / nb,
                          round(100 * above / cnt), round(100 * below / cnt)])
    w = open("output/bias.txt", "w")
    for row in table:
        w.write("\t".join([str(x) for x in row]) + "\n")
    w.close()


def do_hz(model, ltype, loc, N1s):
    title = "%s %s" % (model, ltype)
    print('python plotHz.py "%s" data/trout %s' %
          (title, " ".join([str(N1) + model + "-" + str(loc)
                            for N1 in N1s])))
    os.system('python plotHz.py "%s" data/trout %s' %
              (title, " ".join([str(N1) + model + "-" + str(loc)
                                for N1 in N1s])))
    shutil.move("hz.png", "output/hz-%s-%s.png" % (model, ltype))
    shutil.move("hhz.png", "output/hhz-%s-%s.png" % (model, ltype))
    for N1 in N1s:
        shutil.move("ahz-%d%s-%d.png" % (N1, model, loc),
                    "output/ahz-%s-%d-%s.png" % (model, N1, ltype))
        os.remove("ahz-%d%s-%d.eps" % (N1, model, loc))


def do_hz_comp(pref, mydir, model, N0):
    snps = [100]  # , 200, 400]
    cohort = "Newb"
    cutCase = {}
    for cut in cuts:
        cutCase[cut] = {}
        case = load_file(pref, cut * 100, mydir)
        for nsnp in snps:
            vals, ci, r2, sr2, j, ssize = case[cohort][
                (model, N0)][(None, 50, nsnp, "SNP")]
            cutCase[cut][nsnp] = vals, ci, r2, sr2, j, ssize
    fig, ax = plt.subplots()
    fig.suptitle("Hz comparison: %s - %d - Newb - SNPs - 50 indivs " % (model, N0))
    box = []
    cnt = 0
    labels = []
    tops = []
    bottoms = []
    hmeans = []

    for cut in cuts:
        # ax.text(cnt, 0, str(cut), rotation="vertical")
        cases = cutCase[cut]
        for nsnp in snps:
            vals, ci, r2, sr2, j, ssize = cases[nsnp]
            hmeans.append(hmean(vals))
            bottom, top = zip(*ci)
            top = [100000 if x is None else x for x in top]
            bottom = [100000 if x is None else x for x in bottom]
            tops.append(np.percentile(top, 90))
            bottoms.append(np.percentile(bottom, 10))
            # labels.append(str(nsnp))
            labels.append(str(cut))
            box.append(vals)
            cnt += 1
    sns.boxplot(box, notch=0, sym="",)
    ax.set_ylabel("$\hat{N}_{e}$")
    ax.plot([1 + x for x in range(len(tops))], tops, "rx")
    ax.plot([1 + x for x in range(len(bottoms))], bottoms, "rx")
    ax.plot([1 + x for x in range(len(hmeans))], hmeans, "k+")
    ax.axhline(Nbs[(model, N0)], color="k", lw=0.3)
    ax.set_ylim(0, Nbs[(model, N0)] * 3)
    ax.set_xticks(1 + np.arange(len(labels)))
    ax.set_xticklabels(labels)
    #fig.savefig("output/hz-comp-%s-%d.png" % (model, N0))


def do_pcrit(case, model, N0, isSNP):
    if isSNP:
        markers = [100, 200, 400]
        marker_name = 'SNP'
    else:
        markers = [15, 50, 100]
        marker_name = 'MSAT'
    cohort = "Newb"
    sampCase = {}
    for nmarkers in markers:
        sampCase[nmarkers] = {}
        for pcrit in pcrits:
            #print case[cohort][N0].keys()
            try:
                vals, ci, r2, sr2, j, ssize = case[cohort][
                    (model, N0)][(pcrit, 50, nmarkers, marker_name)]
            except KeyError:
                vals, ci, r2, sr2, j, ssize = [], [], [], [], [], []
            sampCase[nmarkers][pcrit] = vals, ci, r2, sr2, j, ssize
    fig, ax = plt.subplots()
    fig.suptitle("pcrit : %s - %d - Newb - %ss - 50 individuals" % (
        model, N0, marker_name))
    box = []
    cnt = 0
    labels = []
    tops = []
    bottoms = []
    hmeans = []

    for nmarkers in markers:
        ax.text(1 + cnt, 0, str(nmarkers) + " %ss" % marker_name, rotation="vertical", ha="left", va="bottom")
        critCases = sampCase[nmarkers]
        for pcrit in pcrits:
            vals, ci, r2, sr2, j, ssize = critCases[pcrit]
            if len(vals) > 0:
                hmeans.append(hmean(vals))
                bottom, top = zip(*ci)
                tops.append(np.percentile([x if x is not None else 100000
                                           for x in top], 90))
                bottoms.append(np.percentile([x if x is not None else 0.1
                                              for x in bottom], 10))
            else:
                hmeans.append(None)
                tops.append(None)
                bottoms.append(None)
            labels.append(str(pcrit) if pcrit is not None else "std")
            box.append(vals)
            cnt += 1
    sns.boxplot(box, notch=0, sym="",)
    ax.plot([1 + x for x in range(len(tops))], tops, "r.")
    ax.plot([1 + x for x in range(len(bottoms))], bottoms, "r.")
    ax.plot([1 + x for x in range(len(hmeans))], hmeans, "k.")
    ax.set_ylabel("$\hat{N}_{e}$")
    ax.set_ylim(0, Nbs[(model, N0)] * 2)
    ax.axhline(Nbs[(model, N0)], color="k", lw=2)
    ax.set_xticks(1 + np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation="vertical")
    #fig.savefig("output/pcrit-%s-%d.png" % (marker_name, N0))


def _do_window(lst):
    win = []
    for i in range(len(lst)):
        win.append(np.median(lst[max([0, i - 20]): min([len(lst), i + 20])]))
    return win


def do_ld_progress(model, N0s):
    ninds = 50
    strat = "Newb"

    for N0 in N0s:
        fig, ax = plt.subplots()
        fig.suptitle("LDNe estimates %s %d (%d)" %
                     (model, N0, Nbs[(model, N0)]))
        for nsnp in nsnps:
            vals = []
            rep = 0
            try:
                while True:
                    fname = "%s/ldout/%d%s%s%d%d-snp-%d" % (
                        dataDir, N0, model, strat, ninds, nsnp, rep)
                    f = open(fname)
                    myEsts = eval(f.readline().rstrip())
                    for i in range(len(myEsts)):
                        if i >= len(vals):
                            vals.append(myEsts[i])
                        else:
                            vals[i] += myEsts[i]
                    f.close()
                    rep += 1
            except IOError:
                pass
                # We are done
            if rep > 0:
                ax.plot(_do_window(map(lambda x: x / rep, vals)),
                        label=str(nsnp))
        ax.set_ylim(1.0 * Nbs[(model, N0)], 1.5 * Nbs[(model, N0)])
        ax.legend()
        #fig.savefig("output/ldne-%d.png" % N0)


def fetch_nes(model, N0, rep, strat):
    ninds = 50
    nsnp = 100
    fname = "%s/ldout/%d%s%s%d%d-snp-%d" % (dataDir, N0,
                                            model, strat, ninds, nsnp, rep)
    f = open(fname)
    myEsts = eval(f.readline().rstrip())
    l = f.readline().rstrip()
    l = l.replace("inf", "100000")
    cis = eval(l)
    f.close()
    return myEsts, cis


def do_nb_ne(model, N0, rep):
    points = 20
    in_naive, in_corr = 0, 0
    start_year, vals = realNbs[model, N0][rep]
    nes, cis = fetch_nes(model, N0, rep, "Newb")
    fig, ax = plt.subplots()
    fig.suptitle("%s N1: %d Nb: %d Ne: %.2f" % (model, N0,
                                                Nbs[(model, N0)], Nes[N0]))
    ax.plot(vals[start_year:start_year + points], "+")
    shift = 1
    myNes = nes[start_year + shift:start_year + shift + 100]
    myCis = cis[start_year + shift:start_year + shift + 100]
    nes20 = myNes[:points]
    cis20 = myCis[:points]
    cvals, cci = correct_ci(model, 50, nes20, cis20, r2=None, fixed=-0.9)
    errs = []
    cerrs = []
    for i in range(len(cis20)):
        errs.append((nes20[i] - cis20[i][0], cis20[i][1] - nes20[i]))
        if vals[start_year + i] >= cci[i][0] and \
                vals[start_year + i] <= cci[i][1] and cci[i][1] < 100000:
            in_corr += 1
        cerrs.append((cvals[i] - cci[i][0], cci[i][1] - cvals[i]))
    ax.errorbar(np.array(range(len(cis20))) + 0.1, nes20, fmt="+",
                yerr=zip(*errs), mec="red", color="red")
    ax.errorbar(range(len(cis20)), cvals, fmt="+", yerr=zip(*cerrs),
                mec="green", color="green")
    ax.axhline(Nbs[(model, N0)])
    ax.axhline(Nes[N0])
    #fig.savefig("output/nb-ne-%s-%d-%d.png" % (model, N0, rep))

    fig, ax = plt.subplots()
    median = np.median(myNes)
    median = Nbs[(model, N0)]
    ne = Nes[N0]
    basics = []
    nbcomps = []
    necomps = []
    harmcomps = []
    for i in range(len(myNes)):
        basics.append(abs(myNes[i] - median))
        nbcomps.append(abs(myNes[i] - vals[i]))
        hm = 2.0 / (1 / vals[i] + 1 / ne)
        harmcomps.append(abs(myNes[i] - hm))
        necomps.append(abs(myNes[i] - ne))
    ax.plot(basics, "g")
    ax.plot(nbcomps, "r")
    ax.plot(necomps, "b")
    ax.plot(harmcomps, "k")
    #fig.savefig("output/nb-ne-diff-%s-%d-%d.png" % (model, N0, rep))
    #print model, N0, rep, sum(basics) / len(basics), sum(nbcomps) / len(nbcomps), len(basics), len(nbcomps)
    return median, basics, nbcomps, necomps, harmcomps, vals, in_naive, in_corr


def compare_correction_ci(case, model, N0, all_snps, all_indivs):
    cohort = "Newb"
    f, axs = plt.subplots(len(all_snps), len(all_indivs), squeeze=False,
                          sharex=True, sharey=True, figsize=(30, 20))
    bname = get_bname(model)
    Nb = Nbs[(model, N0)]
    top_flex_nb = Nb + 2 * math.sqrt(Nb / 2)
    bottom_flex_nb = Nb - 2 * math.sqrt(Nb / 2)

    def plot_case(ax, nindivs, nsnps):
        case[cohort][(model, N0)][(None, nindivs, nsnps, "SNP")]
        vals, ci, r2, sr2, j, ssize = \
            case[cohort][(model, N0)][(None, nindivs, nsnps, "SNP")]
        bottom_box_vals = []
        top_box_vals = []
        corr_names = []
        ax.axhline(Nb)
        ax.axhline(top_flex_nb)
        ax.axhline(bottom_flex_nb)
        i = 0
        top_y = 3 * Nb
        for corr_name, corrections in get_corrs(N0, bname, nindivs, vals, ci, r2,
                                                sr2, j):
            cvals, cci = corrections
            corr_names.append(corr_name)
            tops, bottoms = zip(*cci)
            top_box_vals.append(tops)
            bottom_box_vals.append(bottoms)
            aboveTop = len([x for x in tops if x is not None and x < Nb])
            belowBottom = len([x for x in bottoms if x is None or x > Nb])
            ax.text(i + 1.05, top_y, "%.1f" % (100 * aboveTop / len(tops)),
                    va='top', ha='left', rotation='vertical',
                    backgroundcolor='white', size=24)
            ax.text(i + 1.05, 0, "%.1f" % (100 * belowBottom / len(bottoms)),
                    va='bottom', ha='left', rotation='vertical',
                    backgroundcolor='white', size=24)
            aboveFlexTop = len([x for x in tops if x is not None and x <
                                bottom_flex_nb])
            belowFlexBottom = len([x for x in bottoms if x is None or x >
                                   top_flex_nb])
            ax.text(i + 0.95, top_y, "%.1f" % (100 * aboveFlexTop / len(tops)),
                    va='top', ha='right', rotation='vertical',
                    backgroundcolor='white', size=24)
            ax.text(i + 0.95, 0, "%.1f" % (100 * belowFlexBottom / len(bottoms)),
                    va='bottom', ha='right', rotation='vertical',
                    backgroundcolor='white', size=24)
            i += 1
        sns.boxplot(top_box_vals, ax=ax)
        sns.boxplot(bottom_box_vals, ax=ax)
        ind = np.arange(len(corr_names))
        ax.set_xticks(ind + 1)
        ax.set_xticklabels(corr_names, rotation="vertical", size=24)
        ax.set_ylim(0, top_y)
        yticks = [0, Nb / 2, Nb, 2 * Nb, 3 * Nb]
        ax.set_yticks(yticks)
        ax.set_yticklabels([str(x) for x in yticks], size=24)
        ax.text(1, Nb, 'inds=%d snps=%d' % (nindivs, nsnps),
                ha='right', va='center', size=36, rotation='vertical')
    for i, n_snps in enumerate(all_snps):
        for j, n_indivs in enumerate(all_indivs):
            try:
                plot_case(axs[i, j], n_indivs, n_snps)
            except KeyError:
                pass  # Might be ok for Nb < sample size
    f.tight_layout()
    #fig.savefig("output/compare-correction-%s-%d-%s.png" % (model, N0, suff))
    return f


def do_table_ci(Nbs, case, dir_pref, modelN0s, nsnps,
                nindivs, ci_percentile=50.0, flex_nb=False):
    cohort = "Newb"
    w = open(dir_pref + "/table-ci-%d-%d-%.1f.txt" % (
        nsnps, nindivs, ci_percentile), "w")
    w.write('Standard deviations problematic because of infinites\n')
    w.write("Corr Model Nb N1 J median mean stdDev percTopCI meanTopCI stdDevTopCI aboveTop medianTopErr percBotCI meanBotCI stdDevBotCI belowBot medianBotErr\n")
    medians = {}
    for bname, model, N0 in modelN0s:
        nb_ = Nbs[(model, N0)]
        if flex_nb:
            top_nb = nb_ + 2 * math.sqrt(nb_ / 2)
            bottom_nb = nb_ - 2 * math.sqrt(nb_ / 2)
        else:
            top_nb = nb_
            bottom_nb = nb_
        vals, ci, r2, sr2, j, ssize = \
            case[cohort][(model, N0)][(None, nindivs, nsnps, "SNP")]
        print(bname, model, N0, top_nb, nb_, bottom_nb)
        for has_corr, corrections in get_corrs(N0, bname, nindivs,
                                               vals, ci, r2, sr2, j):
            cvals, cci = corrections
            topErr = [0, 0.0]
            botErr = [0, 0.0]
            if len(ci) > 0:
                tops, bottoms = zip(*cci)
                topProb = botProb = 0
                for bottom in bottoms:
                    if bottom is None or bottom > top_nb:
                        botProb += 1
                        botErr[0] += 1
                        botErr[1] += bottom - top_nb if bottom is not None else top_nb
                for top in tops:
                    if top is None or top < bottom_nb or top > 100000:
                        topProb += 1
                        if top is None:
                            top = 100000
                        if top < bottom_nb:
                            topErr[0] += 1
                            topErr[1] += bottom_nb - top
                topMean = np.mean([x for x in tops if x is not None])
                botMean = np.mean([x for x in bottoms if x is not None])
                topMedian = np.percentile([x if x is not None else 100000
                                           for x in tops],
                                          100 - ci_percentile)
                botMedian = np.percentile([x if x is not None else 0.1
                                           for x in bottoms],
                                          ci_percentile)
                topStd = np.std([x for x in tops if x is not None])
                botStd = np.std([x for x in bottoms if x is not None])
                topProb /= len(tops)
                botProb /= len(bottoms)
                topProb *= 100
                botProb *= 100
            else:
                topMedian = botMedian = topProb = botProb = topMean = botMean = topStd = botStd = "NA"
            w.write(' '.join([str(x) for x in [has_corr, bname, nb_, N0, np.median(j),
                              np.median(cvals), np.mean(cvals), np.std(cvals),
                              topMedian, topMean, topStd, topProb,
                              np.median(topErr[1]), botMedian, botMean,
                              botStd, botProb, np.median(botErr[1])]]))
            w.write('\n')
            if has_corr == 'NbNe':
                medians[(model, N0)] = nb_, np.median(cvals)
    w.close()
    return medians


def do_all_nb_ne():
    allBasics = {}
    allNbComps = {}
    allNbs = {}
    #in_naives = {}
    in_corrs = {}
    for rep in range(20):
        for N0 in [180, 361, 722]:
            median, basics, nbcomps, necomps, harmcomps, nbs, in_naive, in_corr = do_nb_ne("bulltrout", N0, rep)
            in_corrs[N0] = in_corrs.get(N0, 0) + in_corr
            allBasics.setdefault(N0, []).append(basics)
            allNbComps.setdefault(N0, []).append(nbcomps)
            allNbs.setdefault(N0, []).append(nbs)
    print("naive", in_naive)
    print("corr", in_corrs)
    basicMeans = []
    basicStds = []
    compMeans = []
    compStds = []
    N0s = [180, 361, 722]
    for N0 in N0s:
        print("NbSTD", N0, np.std(allNbs[N0]))
        meanBasic = np.mean(allBasics[N0])
        meanComp = np.mean(allNbComps[N0])
        stdBasic = np.std(allBasics[N0])
        stdComp = np.std(allNbComps[N0])
        basicMeans.append(meanBasic)
        basicStds.append(stdBasic)
        compMeans.append(meanComp)
        compStds.append(stdComp)
    fig, ax = plt.subplots()
    width = 0.35
    ind = np.arange(len(N0s))
    basicRects = ax.bar(ind, basicMeans, width, color="r", yerr=basicStds)
    compRects = ax.bar(ind + width, compMeans, width, color="y", yerr=compStds)
    ax.legend((basicRects, compRects),
              ("$\hat{N}_{e} - \\bar{N}_{b}$", "$\hat{N}_{e} - N_{b}$"))
    ax.set_xticks(ind + width)
    ax.set_xticklabels([str(x) for x in N0s])
    #fig.savefig("output/nb-err.png")


def do_robin_nb_ne():
    model = 'bulltrout'
    allBasics = {}
    allNbComps = {}
    allNeComps = {}
    allHarmComps = {}
    allNbs = {}
    #allNbNes = {}
    for rep in range(20):
        for N0 in [180, 361, 722]:
            median, basics, nbcomps, necomps, harmcomps, nbs, in_naive, in_corr = do_nb_ne("bulltrout", N0, rep)
            allBasics.setdefault(N0, []).append(basics)
            allNbComps.setdefault(N0, []).append(nbcomps)
            allNeComps.setdefault(N0, []).append(necomps)
            allHarmComps.setdefault(N0, []).append(harmcomps)
            allNbs.setdefault(N0, []).append(nbs)
    # basicMeans = []
    # basicStds = []
    # compMeans = []
    # compStds = []
    N0s = [180, 361, 722]
    for N0 in N0s:
        meanBasics = [x for sl in allBasics[N0] for x in sl]
        meanComps = [x for sl in allNbComps[N0] for x in sl]
        neComps = [x for sl in allNeComps[N0] for x in sl]
        harmComps = [x for sl in allHarmComps[N0] for x in sl]

        fig = plt.figure()
        tops = {"X": 0, "Y": 0}
        stats = []

        def plotCase(case, tops):
            ax.hist(case, 50)
            ymin, ymax = ax.get_ylim()
            xmin, xmax = ax.get_xlim()
            if tops["Y"] < ymax:
                tops["Y"] = ymax
            if tops["X"] < xmax:
                tops["X"] = xmax
            return np.mean(case), np.std(case)

        ax = fig.subplot(2, 2, 1)
        stats.append(plotCase(meanBasics, tops))

        ax = fig.subplot(2, 2, 2)
        stats.append(plotCase(meanComps, tops))

        ax = fig.subplot(2, 2, 3)
        stats.append(plotCase(neComps, tops))

        ax = fig.subplot(2, 2, 4)
        stats.append(plotCase(harmComps, tops))

        for i in range(4):
            ax = fig.subplot(2, 2, i + 1)
            ax.set_xlim(0, tops["X"])
            ax.set_ylim(0, tops["Y"])
            if i == 0:
                text = "$\hat{N}_{e} - \\bar{N}_{b}$"
            elif i == 1:
                text = "$\hat{N}_{e} - N_b$"
            elif i == 2:
                text = "$\hat{N}_{e} - N_e$"
            elif i == 3:
                text = "$\hat{N}_{e} - HM(\hat{N}_{e}, N_b)$"

            ax.text(tops["X"], tops["Y"], text, va="bottom", ha="right")
            ax.text(tops["X"], tops["Y"],
                    "%.2f (%.2f)" % tuple(stats[i]), va="top", ha="right")
        fig.text(0.02, 0.5, "bulltrout N1: %d Nb: %d Ne: %.2f" %
                 (N0, Nbs[(model, N0)], Nes[N0]),
                 rotation="vertical", ha="left", va="center")
        #fig.savefig("output/nb-robin-%d.png" % N0)


def do_nb_linear(case, models, name, fun):
    fig, ax = plt.subplots()
    nbs = []
    pes = []
    tops = []
    bottoms = []
    nindivs = 50
    for model, n0 in models:
        nb = Nbs[(model, n0)]
        vals, ci, r2, sr2, j, ssize = \
            case["Newb"][(model, n0)][(None, nindivs, 100, "SNP")]
        vals, ci = fun(n0, get_bname(model), nindivs, vals, ci,
                       r2=r2, sr2=sr2, j=j)
        if len(vals) == 0:
            continue
        bottom, top = zip(*ci)
        nbs.append(nb)
        tops.append(top)
        pes.append(vals)
        bottoms.append(bottom)
    # pylab.yscale('log')
    sns.boxplot(tops, notch=0, sym="")
    sns.boxplot(bottoms)
    ax.set_xticks(1 + np.arange(len(nbs)))
    ax.set_xticklabels([str(nb) for nb in nbs])
    ax.set_ylim(0, max(nbs))
    #fig.savefig("output/nb-linear-%s.png" % name)


def old_main():
    linear = [("bullt2", 305), ("bullt2", 610), ("bullt2", 915),
              ("bullt2", 1220), ("bullt2", 1525), ("bullt2", 1830),
              ("bullt2", 2440), ("bullt2", 3050), ("bullt2", 4575),
              ("bullt2", 6100)]
    do_nb_linear(linear, "None",
                 functools.partial(lambda model, nindivs, x, y: (x, y)))
    #do_nb_linear(linear, "LogQuad",
    #             functools.partial(correct_logquad, abc=[log_a, log_b, log_c]))
    do_nb_linear(linear, "NbNe", functools.partial(correct_ci))
