from __future__ import division

import os
import shutil
import sys

project = sys.argv[1]


def do_hz(model, ltype, loc, N1s):
    title = "%s %s" % (model, ltype)
    print('python plotHz.py "%s" data/%s %s' %
          (title, " ".join([str(N1) + model + "-" + str(loc)
                            for N1 in N1s])))
    os.system('python plotHz.py "%s" data/%s %s' %
              (title, project, " ".join([str(N1) + model + "-" + str(loc)
                                         for N1 in N1s])))
    shutil.move("hz.png", "output/hz-%s-%s.png" % (model, ltype))
    shutil.move("hhz.png", "output/hhz-%s-%s.png" % (model, ltype))
    for N1 in N1s:
        shutil.move("ahz-%d%s-%d.png" % (N1, model, loc),
                    "output/ahz-%s-%d-%s.png" % (model, N1, ltype))
        os.remove("ahz-%d%s-%d.eps" % (N1, model, loc))


try:
    os.remove("output/hz-cut.html")
except OSError:
    pass  # OK
for loc, ltype in [(50, "MSAT"), (100, "SNP")]:
    do_hz("bulltrout", ltype, loc, [86, 174, 353, 713])
    do_hz("bullpred", ltype, loc, [189, 381, 765])
    do_hz("bullt2", ltype, loc, [3040, 6090])
    if ltype == "SNP":
        do_hz("bulltrout", ltype, loc, [86, 174, 353, 713])
        #do_hz("restricted", ltype, loc, [90, 180, 361, 722])
        #do_hz("btrout", ltype, loc, [1619, 6476])
        do_hz("shepard", ltype, loc, [513, 1030])
        do_hz("fraley", ltype, loc, [637, 1278])
        #do_hz("lake", ltype, loc, [18, 36])
        do_hz("wfrog", ltype, loc, [148, 297, 597])
        #do_hz("grizzly", ltype, loc, [23, 46])
        #do_hz("sagegrouse", ltype, loc, [28, 14])
        #do_hz("seaweed", ltype, loc, [70, 35])
        do_hz("synseaweed", ltype, loc, [81, 164, 207])
        do_hz("mosquito", ltype, loc, [102, 167, 210, 428])
        do_hz("bullt2", ltype, loc, [606, 1214, 1822, 2433, 3040, 4565, 6090])
shutil.move("hz-cut.html", "output/hz-cut.html")
shutil.move("hz-cut.txt", "output/hz-cut.txt")
