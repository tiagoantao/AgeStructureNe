from __future__ import division

import os
import shutil
import sys


pref = sys.argv[1]


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
        do_hz("lake", ltype, loc, [18, 36])
        do_hz("wfrog", ltype, loc, [600, 300, 150])
        do_hz("grizzly", ltype, loc, [23, 46])
        do_hz("sagegrouse", ltype, loc, [28, 14])
        do_hz("seaweed", ltype, loc, [70, 35])
        do_hz("mosquito", ltype, loc, [436, 218, 110])
        do_hz("bullt2", ltype, loc, [305, 610, 915, 1220, 1525,
                                     1830, 2440, 3050, 4575, 6100])
shutil.move("hz-cut.html", "output/hz-cut.html")
shutil.move("hz-cut.txt", "output/hz-cut.txt")
