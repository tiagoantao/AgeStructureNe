from __future__ import division, print_function

import os
import sys

import matplotlib
matplotlib.use('AGG')
import pylab

if len(sys.argv) < 4:
    print(("python %s title baseDir models..." % (sys.argv[0],)))
    sys.exit(-1)

title = sys.argv[1]
baseDir = sys.argv[2]
models = sys.argv[3:]

h_dist = {}

all_sims = {}

cuts = [0.45, 0.4, 0.35, 0.3, 0.25]

fnames = os.listdir(baseDir)
hzf = open("hz-cut.html", "a")
hztf = open("hz-cut.txt", "a")
for model in models:
    allLs = []
    all_sims[model] = {}
    start = True
    for fname in fnames:
        if fname.startswith(model) and fname.endswith(".hz"):
            rep = fname.split("-")[-1][:-3]
            all_sims[model][rep] = []
            f = open(baseDir + "/" + fname)
            ls = f.readlines()
            f.close()
            cnt = 0
            for l in ls:
                res = [float(x) for x in
                       l.rstrip().replace("\t", " ").split(" ")]
                all_sims[model][rep].append(sum(res) / len(res))
                if cnt == 200:
                    bla = res
                if start:
                    allLs.append(res)
                else:
                    allLs[cnt].extend(res)
                cnt += 1
            # last line
            if start:
                h_dist[model] = []
            h_dist[model].extend(bla)
            start = False
    my = []
    monom = []
    for hzs in allLs:
        my.append(sum(hzs) / len(hzs))
        monom.append(len([x for x in hzs if x == 0.0]) / len(hzs))
    currMean = my[0]
    for g in range(len(my) - 1):
        gen = g + 1
        meanHz = my[gen]
        for cut in cuts:
            if cut < currMean and cut > meanHz:
                hzf.write("<tr><td>%s</td><td>%.2f</td><td>%d</td></tr>\n" %
                          (model, cut, gen))
                hztf.write("%s\t%.2f\t%d\n" % (model, cut, gen))
        currMean = meanHz
    pylab.plot(my, label=model.split("-")[0])
    pylab.plot(monom, label="Monom " + model.split("-")[0])
hzf.close()
hztf.close()
pylab.ylim(0, 0.5)
pylab.title(title)
pylab.legend(loc=3)
pylab.savefig("hz.eps")
pylab.savefig("hz.png")

pylab.clf()
pylab.title(title)
for model, hd in list(h_dist.items()):
    val_cnt = {}
    for val in hd:
        val = round(val, 2)
        val_cnt[val] = val_cnt.get(val, 0) + 1
    freqs = []
    cnts = []
    for freq, cnt in list(val_cnt.items()):
        freqs.append(freq)
        cnts.append(cnt)
    pylab.plot(freqs, cnts, '.', label=model)
pylab.legend()
pylab.xlim(-0.01, 0.51)
pylab.savefig("hhz.eps")
pylab.savefig("hhz.png")

for model, reps in list(all_sims.items()):
    pylab.clf()
    pylab.title("%s" % model)
    for rep, vals in list(reps.items()):
        pylab.plot(list(range(len(vals))), vals, label=rep)
    pylab.legend()
    pylab.savefig("ahz-%s.eps" % model)
    pylab.savefig("ahz-%s.png" % model)
    pylab.gca().legend().set_visible(False)
