from __future__ import division

import functools
import os
import shutil
import sys

import numpy
from scipy.stats import hmean

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

import pylab

from trout import *

pref = sys.argv[1]


def do_nb(cohort, nsnps, pref):
    model = 'bulltrout'
    for N0 in N0s:
        pylab.clf()
        pylab.title("Nb: %d (N1: %d) - cohort %s" % (Nbs[model, N0], N0, cohort))
        box_vals = []
        tops = []
        bottoms = []
        labels = []
        hmeans = []
        for nindiv in nindivs:
            for nloci in nlocis:
                try:
                    vals, ci, r2, poli = case[cohort][(model, N0)][(None, nindiv, nloci, "MSAT")]
                    if len(ci) > 0:
                        bottom, top = zip(*ci)
                        tops.append(numpy.percentile(top, 90))
                        bottoms.append(numpy.percentile(bottom, 10))
                    else:
                        tops.append(None)
                        bottoms.append(None)
                except KeyError:
                    vals = []
                    tops.append(None)
                    bottoms.append(None)
                box_vals.append(vals)
                hmeans.append(hmean(vals))
                labels.append("%dMS" % nloci)
            for nsnp in nsnps:
                try:
                    vals, ci, r2, poli = case[cohort][(model, N0)][(None, nindiv, nsnp, "SNP")]
                    if len(ci) > 0:
                        bottom, top = zip(*ci)
                        tops.append(numpy.percentile(top, 90))
                        bottoms.append(numpy.percentile(bottom, 10))
                    else:
                        tops.append(None)
                        bottoms.append(None)

                except KeyError:
                    vals = []
                    tops.append(None)
                    bottoms.append(None)
                box_vals.append(vals)
                hmeans.append(hmean(vals))
                labels.append("%d" % nsnp)
            pos = len(labels)
            pylab.axvline(pos + 0.5, color="k", lw=0.2)
            pylab.text(pos - 0.5, 0, "%d Indivs" % nindiv,
                       ha="center", va="bottom", size="small",
                       rotation="horizontal")
        pylab.ylabel("$\hat{N}_{e}$")
        pylab.ylim(0, Nbs[(model, N0)] * 3)
        pylab.axhline(Nbs[(model, N0)], color="k", lw=0.3)
        pylab.xticks(range(len(labels)), labels, rotation="vertical")
        pylab.boxplot(box_vals, notch=0, sym="")
        pylab.plot([1 + x for x in range(len(tops))], tops, "rx", ms=20)
        pylab.plot([1 + x for x in range(len(bottoms))], bottoms, "rx", ms=20)
        pylab.plot([1 + x for x in range(len(hmeans))], hmeans, "k+", ms=20)
        pylab.savefig("output/%s%s-%d.png" % (pref, cohort, N0))


def do_cohort(model, N0, nindiv):
    last = 0.5
    pylab.clf()
    pylab.title("Nb: %d (N1: %d) - different cohorts - 100 SNPs" %
                (Nbs[(model, N0)], N0))
    box_vals = []
    labels = []
    tops = []
    bottoms = []
    hmeans = []

    for cohort in cohorts:
        vals, ci, r2, poli = case[cohort][(model, N0)][(None, nindiv, 100, "SNP")]
        box_vals.append(vals)
        hmeans.append(hmean(vals))
        bottom, top = zip(*ci)
        tops.append(numpy.percentile(top, 90))
        bottoms.append(numpy.percentile(bottom, 10))
        labels.append("%s" % cohort)
        if cohort == cohorts[-1]:
            pos = len(labels) + 0.5
            pylab.axvline(pos, color="k", lw=0.2)
            pylab.text(last + (pos - last) / 2, 0, "%d Indivs" % nindiv,
                       ha="center", va="bottom", size="small",
                       rotation="horizontal")
            last = pos
    pylab.ylim(0, Nbs[(model, N0)] * 3)
    pylab.ylabel("$\hat{N}_{e}$")
    pylab.axhline(Nbs[(model, N0)], color="k", lw=0.3)
    pylab.xticks(range(len(labels)), labels)
    pylab.boxplot(box_vals, notch=0, sym="")
    pylab.plot([1 + x for x in range(len(tops))], tops, "rx")
    pylab.plot([1 + x for x in range(len(bottoms))], bottoms, "rx")
    pylab.plot([1 + x for x in range(len(hmeans))], hmeans, "k+")
    pylab.savefig("output/cohort-%s-%d.png" % (model, N0))


def do_rel(model, N0, nindiv):
    #last = 0.5
    pylab.clf()
    cohort = "All"
    pylab.title("Nb: %d (N1: %d) - 20pc related individuals - cohort %s" %
                (Nbs[(model, N0)], N0, cohort))
    box_vals = []
    labels = []
    #vals, ci, r2, poli = case[cohort][N0][(None, nindiv, 15, "MSAT-rel")]
    #box_vals.append(vals)
    #labels.append("MSAT-15-rel")
    #vals, ci, r2, poli = case[cohort][N0][(None, nindiv, 15, "MSAT")]
    #box_vals.append(vals)
    #labels.append("MSAT-15")
    vals, ci, r2, poli = case[cohort][(model, N0)][(None, nindiv, 100, "SNP-rel")]
    box_vals.append(vals)
    labels.append("SNP-100-rel")
    vals, ci, r2, poli = case[cohort][(model, N0)][(None, nindiv, 100, "SNP")]
    box_vals.append(vals)
    labels.append("SNP-100")
    pylab.ylim(0, Nbs[(model, N0)] * 3)
    pylab.ylabel("$\hat{N}_{e}$")
    pylab.axhline(Nbs[(model, N0)], color="k", lw=0.3)
    pylab.xticks(range(len(labels)), labels)
    pylab.boxplot(box_vals, notch=0, sym="")
    pylab.savefig("output/rel-%s-%d.png" % (model, N0))


def do_nb_comp():
    pylab.clf()
    box_vals = []
    labels = []
    n0s = Nbs.keys()
    n0s.sort()
    pylab.title("Nb comparison - 100 SNPs - 50 indivs")
    for model, n0 in n0s:
        if type(n0) != int:
            continue
        nb = Nbs[(model, n0)]
        vals, ci, r2, poli = case["Newb"][(model, n0)][(None, 50, 100, "SNP")]
        labels.append(str(nb))
        box_vals.append(vals)
    pylab.xticks(range(len(labels)), labels)
    pylab.ylabel("$\hat{N}_{e}$")
    pylab.boxplot(box_vals, notch=0, sym="")
    pylab.savefig("output/nb-comp.png")


def do_lt_comp(nb, strat):
    pylab.clf()
    tops = []
    box_vals = []
    bottoms = []
    hmeans = []
    labels = []
    n0s = Nbs.keys()
    n0s.sort()
    pylab.title("Nb: %d (Life table comparison)" % nb)
    pylab.suptitle("%s sampled - 100 SNPs - 50 individuals" % strat)
    for pref, name in NbNames:
        for k, nb2 in Nbs.items():
            model, n0 = k
            if nb2 != nb:
                continue
            elif pref != model:
                continue
            vals, ci, r2, poli = case[strat][(model, n0)][(None, 50, 100, "SNP")]
            hmeans.append(hmean(vals))
            try:
                bottom, top = zip(*ci)
            except ValueError:
                top, bottom = [], []
            labels.append(name)
            try:
                tops.append(numpy.percentile(top, 90))
                box_vals.append(vals)
                bottoms.append(numpy.percentile(bottom, 10))
            except ValueError:
                tops.append(None)
                box_vals.append([])
                bottoms.append(None)
    pylab.xticks(range(len(labels)), labels)
    pylab.plot([1 + x for x in range(len(tops))], tops, "rx", ms=20)
    pylab.plot([1 + x for x in range(len(bottoms))], bottoms, "rx", ms=20)
    pylab.plot([1 + x for x in range(len(hmeans))], hmeans, "k+", ms=20)
    bp = pylab.boxplot(box_vals, notch=0, sym="")
    ymin, ymax = pylab.ylim()
    pylab.ylabel("$\hat{N}_{e}$")
    pylab.ylim(ymin, min([ymax, 3 * nb]))
    pylab.axhline(nb, color="k", lw=0.3)
    pylab.setp(bp['boxes'], linewidth=1.0, color="k")
    pylab.setp(bp['caps'], linewidth=1.0, color="k")
    pylab.setp(bp['fliers'], linewidth=1.0, color="k")
    pylab.savefig("output/lt-comp-%s-%d.png" % (strat, nb))


def do_bias():
    #box_vals = []
    #labels = []
    n0s = Nbs.keys()
    n0s.sort()
    sampling = "Newb"
    table = [["Nb", "Sampling", "Model", "N1", "NeEst", "NeEst/Nb", "Above", "Below"]]
    for model, n0 in n0s:
        nb = Nbs[(model, n0)]
        vals, ci, r2, poli = case[sampling][(model, n0)][(None, 50, 100, "SNP")]
        for p, bname in NbNames:
            if p == model:
                break
        cnt = 0
        below = 0
        above = 0
        for i in range(len(vals)):
            #v = vals[i]
            mi, ma = ci[i]
            if nb < mi:
                below += 1
            if nb > ma:
                above += 1
            cnt += 1
        if cnt > 0:
            table.append([nb, sampling, bname, n0, numpy.median(vals),
                          numpy.median(vals) / nb,
                          round(100 * above / cnt), round(100 * below / cnt)])
    w = open("output/bias.txt", "w")
    for row in table:
        w.write("\t".join([str(x) for x in row]) + "\n")
    w.close()


def do_hz(model, ltype, loc, N1s):
    title = "%s %s" % (model, ltype)
    os.system('python plotHz.py "%s" data/trout %s' %
              (title, " ".join([str(N1) + model + "-" + str(loc) for N1 in N1s])))
    shutil.move("hz.png", "output/hz-%s-%s.png" % (model, ltype))
    shutil.move("hhz.png", "output/hhz-%s-%s.png" % (model, ltype))
    for N1 in N1s:
        shutil.move("ahz-%d%s-%d.png" % (N1, model, loc), "output/ahz-%s-%d-%s.png" % (model, N1, ltype))
        os.remove("ahz-%d%s-%d.eps" % (N1, model, loc))


def do_hz_comp(model, N0):
    snps = [100]  # , 200, 400]
    cohort = "Newb"
    cutCase = {}
    for cut in cuts:
        cutCase[cut] = {}
        case = load_file(pref, cut * 100)
        for nsnps in snps:
            vals, ci, r2, poli = case[cohort][(model, N0)][(None, 50, nsnps, "SNP")]
            cutCase[cut][nsnps] = vals, ci, r2, poli
    pylab.clf()
    pylab.title("Hz comparison: %s - %d " % (model, N0))
    pylab.suptitle("Newb sampled - SNPs - 50 individuals")
    box = []
    cnt = 0
    labels = []
    tops = []
    bottoms = []
    hmeans = []

    for cut in cuts:
        #pylab.text(cnt, 0, str(cut), rotation="vertical")
        cases = cutCase[cut]
        for nsnps in snps:
            vals, ci, r2, poli = cases[nsnps]
            hmeans.append(hmean(vals))
            bottom, top = zip(*ci)
            tops.append(numpy.percentile(top, 90))
            bottoms.append(numpy.percentile(bottom, 10))
            #labels.append(str(nsnps))
            labels.append(str(cut))
            box.append(vals)
            cnt += 1
    pylab.xticks(range(len(labels)), labels)
    pylab.boxplot(box, notch=0, sym="",)
    pylab.ylabel("$\hat{N}_{e}$")
    pylab.plot([1 + x for x in range(len(tops))], tops, "rx")
    pylab.plot([1 + x for x in range(len(bottoms))], bottoms, "rx")
    pylab.plot([1 + x for x in range(len(hmeans))], hmeans, "k+")
    pylab.axhline(Nbs[(model, N0)], color="k", lw=0.3)
    pylab.ylim(0, Nbs[(model, N0)] * 3)
    pylab.savefig("output/hz-comp.png")


def do_pcrit(model, N0, isSNP):
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
                vals, ci, r2, poli = case[cohort][(model, N0)][(pcrit,
                                                                50, nmarkers,
                                                                marker_name)]
            except KeyError:
                vals, ci, r2, poli = [], [], [], []
            sampCase[nmarkers][pcrit] = vals, ci, r2, poli
    pylab.clf()
    pylab.title("pcrit comparison: %s - %d " % (model, N0))
    pylab.suptitle("Newb sampled - %ss - 50 individuals" % marker_name)
    box = []
    cnt = 0
    labels = []
    tops = []
    bottoms = []
    hmeans = []

    for nmarkers in markers:
        pylab.text(cnt + 1, 0, str(nmarkers) + " %ss" % marker_name, rotation="vertical", ha="left", va="bottom")
        critCases = sampCase[nmarkers]
        for pcrit in pcrits:
            vals, ci, r2, poli = critCases[pcrit]
            if len(vals) > 0:
                hmeans.append(hmean(vals))
                bottom, top = zip(*ci)
                tops.append(numpy.percentile(top, 90))
                bottoms.append(numpy.percentile(bottom, 10))
            else:
                hmeans.append(None)
                tops.append(None)
                bottoms.append(None)
            labels.append(str(pcrit) if pcrit is not None else "std")
            box.append(vals)
            cnt += 1
    pylab.axhline(Nbs[(model, N0)], color="k", lw=1)
    pylab.xticks(range(len(labels)), labels, rotation="vertical")
    pylab.boxplot(box, notch=0, sym="",)
    pylab.plot([1 + x for x in range(len(tops))], tops, "rx")
    pylab.plot([1 + x for x in range(len(bottoms))], bottoms, "rx")
    pylab.plot([1 + x for x in range(len(hmeans))], hmeans, "k+")
    pylab.ylabel("$\hat{N}_{e}$")
    pylab.ylim(0, Nbs[(model, N0)] * 3)
    pylab.savefig("output/pcrit-%s-%d.png" % (marker_name, N0))


def _do_window(lst):
    win = []
    for i in range(len(lst)):
        win.append(numpy.median(lst[max([0, i - 20]): min([len(lst), i + 20])]))
    return win


def do_ld_progress(model, N0s):
    ninds = 50
    strat = "Newb"

    for N0 in N0s:
        pylab.clf()
        pylab.title("LDNe estimates %s %d (%d)" % (model, N0, Nbs[(model, N0)]))

        for nsnp in nsnps:
            vals = []
            rep = 0
            try:
                while True:
                    fname = "%s/ldout/%d%s%s%d%d-snp-%d" % (dataDir, N0, model, strat,
                                                            ninds, nsnp, rep)
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
                pylab.plot(_do_window(map(lambda x: x / rep, vals)), label=str(nsnp))
        pylab.ylim(1.0 * Nbs[(model, N0)], 1.5 * Nbs[(model, N0)])
        pylab.legend()
        pylab.savefig("output/ldne-%d.png" % N0)


def fetch_nes(model, N0, rep, strat):
    ninds = 50
    nsnp = 100
    fname = "%s/ldout/%d%s%s%d%d-snp-%d" % (dataDir, N0, model, strat, ninds, nsnp, rep)
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
    pylab.clf()
    pylab.title("%s N1: %d Nb: %d Ne: %.2f" % (model, N0, Nbs[(model, N0)], Nes[N0]))
    pylab.plot(vals[start_year:start_year + points], "+")
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
        if vals[start_year + i] >= cci[i][0] and vals[start_year + i] <= cci[i][1] and cci[i][1] < 10000:
            in_corr += 1
        cerrs.append((cvals[i] - cci[i][0], cci[i][1] - cvals[i]))
    pylab.errorbar(numpy.array(range(len(cis20))) + 0.1, nes20, fmt="+", yerr=zip(*errs), mec="red", color="red")
    pylab.errorbar(range(len(cis20)), cvals, fmt="+", yerr=zip(*cerrs), mec="green", color="green")
    pylab.axhline(Nbs[(model, N0)])
    pylab.axhline(Nes[N0])
    pylab.savefig("output/nb-ne-%s-%d-%d.png" % (model, N0, rep))

    pylab.clf()
    median = numpy.median(myNes)
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
    pylab.plot(basics, "g")
    pylab.plot(nbcomps, "r")
    pylab.plot(necomps, "b")
    pylab.plot(harmcomps, "k")
    pylab.savefig("output/nb-ne-diff-%s-%d-%d.png" % (model, N0, rep))
    #print model, N0, rep, sum(basics) / len(basics), sum(nbcomps) / len(nbcomps), len(basics), len(nbcomps)
    return median, basics, nbcomps, necomps, harmcomps, vals, in_naive, in_corr


def do_table_ci(modelN0s, nsnps):
    #table = []
    cohort = "Newb"
    thres = 10
    nindivs = 50
    w = open("output/table-ci-%d.txt" % nsnps, "w")
    print >>w, "Corr Model N1 median mean stdDev medianTopCI meanTopCI stdDevTopCI aboveTop probTop medianTopErr medianBotCI meanBotCI stdDevBotCI belowBot probBot medianBotErr"
    for bname, model, N0 in modelN0s:
        nb = Nbs[(model, N0)]
        vals, ci, r2, poli = case[cohort][(model, N0)][(None, nindivs, nsnps, "SNP")]
        errs = []
        for v in vals:
            perc_err = abs(v - nb) / nb
            errs.append(perc_err)
        print bname, model, N0, nb, numpy.median(errs)
        for has_corr, corrections in get_corrs(bname, nindivs,
                                               vals, ci, r2, poli):
            cvals, cci = corrections
            topErr = [0, 0.0]
            botErr = [0, 0.0]
            probTop = 0
            probBot = 0
            if len(ci) > 0:
                bottoms, tops = zip(*cci)
                topProb = botProb = 0
                for bottom in bottoms:
                    if bottom > nb:
                        botProb += 1
                        botErr[0] += 1
                        botErr[1] += bottom - nb
                        if bottom - nb > thres:
                            probBot += 1
                for top in tops:
                    if top < nb or top > 10000:
                        topProb += 1
                        if top < nb:
                            topErr[0] += 1
                            topErr[1] += nb - top
                            if nb - top > thres:
                                probTop += 1
                topMean = numpy.mean(tops)
                botMean = numpy.mean(bottoms)
                topMedian = numpy.median(tops)
                botMedian = numpy.median(bottoms)
                topStd = numpy.std(tops)
                botStd = numpy.std(bottoms)
                topProb /= len(tops)
                botProb /= len(bottoms)
                topProb *= 100
                botProb *= 100
                probTop /= len(tops)
                probBot /= len(bottoms)
                probTop *= 100
                probBot *= 100
            else:
                topMedian = botMedian = topProb = botProb = topMean = botMean = topStd = botStd = probBot = probTop = "NA"
            print >>w, has_corr, bname, N0, numpy.median(cvals), numpy.mean(cvals), numpy.std(cvals),
            print >>w, topMedian, topMean, topStd, topProb, probTop, numpy.median(topErr[1]),
            print >>w, botMedian, botMean, botStd, botProb, probBot, numpy.median(botErr[1])
    w.close()


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
    print "naive", in_naive
    print "corr", in_corrs
    basicMeans = []
    basicStds = []
    compMeans = []
    compStds = []
    N0s = [180, 361, 722]
    for N0 in N0s:
        print "NbSTD", N0, numpy.std(allNbs[N0])
        meanBasic = numpy.mean(allBasics[N0])
        meanComp = numpy.mean(allNbComps[N0])
        stdBasic = numpy.std(allBasics[N0])
        stdComp = numpy.std(allNbComps[N0])
        basicMeans.append(meanBasic)
        basicStds.append(stdBasic)
        compMeans.append(meanComp)
        compStds.append(stdComp)
    pylab.clf()
    fig, ax = pylab.subplots()
    width = 0.35
    ind = numpy.arange(len(N0s))
    basicRects = pylab.bar(ind, basicMeans, width, color="r", yerr=basicStds)
    compRects = pylab.bar(ind + width, compMeans, width, color="y", yerr=compStds)
    pylab.legend((basicRects, compRects),
                 ("$\hat{N}_{e} - \\bar{N}_{b}$", "$\hat{N}_{e} - N_{b}$"))
    ax.set_xticks(ind + width)
    ax.set_xticklabels([str(x) for x in N0s])
    pylab.savefig("output/nb-err.png")


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
    #basicMeans = []
    #basicStds = []
    #compMeans = []
    #compStds = []
    N0s = [180, 361, 722]
    for N0 in N0s:
        meanBasics = [x for sl in allBasics[N0] for x in sl]
        meanComps = [x for sl in allNbComps[N0] for x in sl]
        neComps = [x for sl in allNeComps[N0] for x in sl]
        harmComps = [x for sl in allHarmComps[N0] for x in sl]

        pylab.clf()

        tops = {"X": 0, "Y": 0}
        stats = []

        def plotCase(case, tops):
            pylab.hist(case, 50)
            ymin, ymax = pylab.ylim()
            xmin, xmax = pylab.xlim()
            if tops["Y"] < ymax:
                tops["Y"] = ymax
            if tops["X"] < xmax:
                tops["X"] = xmax
            return numpy.mean(case), numpy.std(case)

        pylab.subplot(2, 2, 1)
        stats.append(plotCase(meanBasics, tops))

        pylab.subplot(2, 2, 2)
        stats.append(plotCase(meanComps, tops))

        pylab.subplot(2, 2, 3)
        stats.append(plotCase(neComps, tops))

        pylab.subplot(2, 2, 4)
        stats.append(plotCase(harmComps, tops))

        for i in range(4):
            pylab.subplot(2, 2, i + 1)
            pylab.xlim(0, tops["X"])
            pylab.ylim(0, tops["Y"])
            if i == 0:
                text = "$\hat{N}_{e} - \\bar{N}_{b}$"
            elif i == 1:
                text = "$\hat{N}_{e} - N_b$"
            elif i == 2:
                text = "$\hat{N}_{e} - N_e$"
            elif i == 3:
                text = "$\hat{N}_{e} - HM(\hat{N}_{e}, N_b)$"

            pylab.text(tops["X"], tops["Y"], text, va="bottom", ha="right")
            pylab.text(tops["X"], tops["Y"], "%.2f (%.2f)" % tuple(stats[i]), va="top", ha="right")
        pylab.figtext(0.02, 0.5, "bulltrout N1: %d Nb: %d Ne: %.2f" % (N0,
                                                                       Nbs[(model,
                                                                            N0)], Nes[N0]),
                      rotation="vertical", ha="left", va="center")
        pylab.savefig("output/nb-robin-%d.png" % N0)


def do_nb_linear(models, name, fun):
    pylab.clf()
    nbs = []
    pes = []
    tops = []
    bottoms = []
    nindivs = 50
    for model, n0 in models:
        nb = Nbs[(model, n0)]
        print model, n0, case["Newb"][(model, n0)].keys()
        vals, ci, r2, poli = case["Newb"][(model, n0)][(None, nindivs, 100, "SNP")]
        vals, ci = fun(model, nindivs, vals, ci)
        if len(vals) == 0:
            continue
        bottom, top = zip(*ci)
        nbs.append(nb)
        tops.append(numpy.median(top))
        pes.append(numpy.median(vals))
        bottoms.append(numpy.median(bottom))
    #pylab.yscale('log')
    pylab.plot(nbs, nbs, "k", label="Nb")
    pylab.plot(nbs, nbs, "ko")
    pylab.plot(nbs, tops, "b", label="Top")
    pylab.plot(nbs, tops, "bo-")
    pylab.plot(nbs, pes, "g", label="PE")
    pylab.plot(nbs, pes, "go")
    pylab.plot(nbs, bottoms, "b", label="Bottom")
    pylab.plot(nbs, bottoms, "bo-")
    pylab.legend()
    pylab.savefig("output/nb-linear-%s.png" % name)

load_nb(pref)

do_robin_nb_ne()

do_all_nb_ne()

do_ld_progress("bulltrout", [180, 361, 722])

try:
    os.remove("output/hz-cut.html")
except OSError:
    pass  # OK
for loc, ltype in [(0, "MSAT"), (100, "SNP")]:
    do_hz("bulltrout", ltype, loc, [180, 361, 722])
    do_hz("bullpred", ltype, loc, [193, 387, 775])
    do_hz("bullt2", ltype, loc, [3050, 6100])
    if ltype == "SNP":
        do_hz("bulltrout", ltype, loc, [90])
        do_hz("restricted", ltype, loc, [90, 180, 361, 722])
        do_hz("btrout", ltype, loc, [1619, 6476])
        do_hz("shepard", ltype, loc, [518, 1036])
        do_hz("fraley", ltype, loc, [641, 1282])
        do_hz("lake", ltype, loc, [18, 72])
        do_hz("bullt2", ltype, loc, [305, 610, 915, 1220, 1525,
                                     1830, 2440, 3050, 4575, 6100])
shutil.move("hz-cut.html", "output/hz-cut.html")
shutil.move("hz-cut.txt", "output/hz-cut.txt")

do_hz_comp("bulltrout", 361)

case = load_file(pref, 45)
do_lt_comp(200, "Newb")
do_lt_comp(200, "All")
case = load_file(pref, 40)
linear = [("bullt2", 305), ("bullt2", 610), ("bullt2", 915), ("bullt2", 1220),
          ("bullt2", 1525), ("bullt2", 1830), ("bullt2", 2440),
          ("bullt2", 3050), ("bullt2", 4575), ("bullt2", 6100)]
do_nb_linear(linear, "std",
             functools.partial(lambda model, nindivs, x, y: (x, y)))
do_nb_linear(linear, "Int0.9",
             functools.partial(correct_ci, fixed=-0.9))
#do_lt_comp(50, "All")
#do_lt_comp(100, "All")

do_lt_comp(25, "Newb")
do_lt_comp(50, "Newb")
do_lt_comp(100, "Newb")

do_pcrit("bulltrout", 180, True)
do_pcrit("bulltrout", 180, False)
do_pcrit("bulltrout", 361, True)
do_pcrit("bulltrout", 361, False)

cis = [("BuTrout", "bulltrout", 90), ("BuTrout", "bulltrout", 180),
       ("BuTrout", "bulltrout", 361), ("BuTrout", "bulltrout", 722),
       ("WCT-S", "shepard", 518), ("WCT-S", "shepard", 1036),
       ("WCT-F", "fraley", 641), ("WCT-F", "fraley", 1282),
       ("BuLong", "bullt2", 305), ("BuLong", "bullt2", 610),
       ("BuLong", "bullt2", 915), ("BuLong", "bullt2", 1220),
       ("BuLong", "bullt2", 1525), ("BuLong", "bullt2", 2440),
       ("BuLong", "bullt2", 3050), ("BuLong", "bullt2", 4575),
       ("BuLong", "bullt2", 6100),
       ("BuPred", "bullpred", 193), ("BuPred", "bullpred", 387)]

do_table_ci(cis, 100)
do_table_ci(cis, 200)

for cohort in cohorts:
    do_nb(cohort, [100], "")
    do_nb(cohort, nsnps, "allsnps-")

do_cohort("bulltrout", 361, 50)
do_cohort("bulltrout", 180, 50)

do_rel("bulltrout", 361, 50)
#do_rel("bulltrout", 180, 50)

do_nb_comp()
do_bias()
