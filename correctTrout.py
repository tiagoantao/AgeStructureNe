from __future__ import division

import sys

from scipy import stats

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

import pylab

from trout import *

pref = sys.argv[1]

case = load_file(pref, 45)


def plot_model(fig, model):
    bname = get_bname(model)
    vals = []
    ldnes = {}
    errs = {}
    labels = []
    cnb = case['Newb']
    nbks = sorted(Nbs.keys(), key=lambda x: x[1])
    for cname, cdata in get_corrs(bname, [], []):
        ldnes[cname] = []
        errs[cname] = []

    for name, N0 in nbks:
        if name != model:
            continue
        labels.append(str(N0))
        val = Nbs[(model, N0)]
        ldne, ci = cnb[(model, N0)][None, 50, 100, 'SNP']
        for cname, cdata in get_corrs(bname, ldne, ci):
            cldne, ccis = cdata
            hmean = stats.hmean([x if x > 0 else 10000 for x in cldne])
            ldnes[cname].append(hmean)
            err = hmean / val
            errs[cname].append(err)
        vals.append(val)

    ax = fig.add_subplot(2, 1, 1)
    ax.set_title("Nb and estimators %s" % bname)
    ax.plot(vals, '+')
    for name, lvals in ldnes.items():
        ax.plot(lvals, '-', label=name)
        print name
        print vals
        print lvals
    ax.set_xticklabels(labels)
    ax.legend()

    ax = fig.add_subplot(2, 1, 2)
    ax.set_title("Fraction of error %s" % bname)
    for name, cvals in errs.items():
        ax.plot(cvals, '-', label=name)
    ax.set_xticklabels(labels)
    ax.legend()

    fig.savefig("output/correct.png")

fig = pylab.figure(figsize=(10, 20))
plot_model(fig, 'bullt2')
