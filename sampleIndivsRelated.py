from __future__ import division
from sys import stdin, argv, exit
import random

#Only from newborns

if len(argv) not in [5]:
    print "python %s frac indivs startGen endGen" % (argv[0],)
    print "indivs can be ALL"
    exit(-1)

frac = int(argv[1]) / 100
if argv[2] == "ALL":
    nindivs = None  # Actually we stop at 200
    maxInds = 200
else:
    nindivs = int(argv[2])
startGen = int(argv[3])
endGen = int(argv[4])

l = stdin.readline()

currGen = 0


def getBrothasAndSistas(sameParents, todo):
    rels = []
    done_parents = []
    while todo > 0:
        have_more = False
        for parents, vals in sameParents.items():
            if parents in done_parents:
                continue
            done_parents.append(parents)
            if len(vals) > 1:
                have_more = True
                rels.append(vals[0])
                rels.append(vals[1])
                todo -= 2
                break
        if not have_more:
            print "Not enough relateds"
            exit()
    return rels

sameParents = {}
indivs = []
while l != "":
    toks = l.rstrip().split(" ")
    rep = int(float(toks[1]))
    if rep > 0:
        break
    gen = int(float(toks[0]))
    id = int(float(toks[2]))
    sex = int(float(toks[3]))
    father = int(float(toks[4]))
    mother = int(float(toks[5]))
    age = int(float(toks[6]))
    if age > 1:
        l = stdin.readline()
        continue
    sameParents.setdefault((father, mother), []).append(id)
    #sameParents.setdefault(father, []).append(id)
    #sameParents.setdefault(mother, []).append(id)
    indivs.append(id)

    if gen > currGen:
        if currGen >= startGen:
            relateds = getBrothasAndSistas(sameParents, frac * nindivs)
            for related in relateds:
                indivs.remove(related)
            if nindivs:
                remaining = nindivs - len(relateds)
                print currGen, relateds + random.sample(indivs, remaining)
            else:
                if len(indivs) > maxInds:
                    remaining = maxInds - len(relateds)
                    print currGen, relateds + random.sample(indivs, remaining)
                else:
                    print currGen, relateds + indivs
        sameParents = {}
        indivs = []
        currGen = gen
        if currGen > endGen:
            break

    l = stdin.readline()

if currGen <= endGen:
    relateds = getBrothasAndSistas(sameParents, frac * nindivs)
    for related in relateds:
        indivs.remove(related)
    if nindivs:
        remaining = nindivs - len(relateds)
        print currGen, relateds + random.sample(indivs, remaining)
    else:
        if len(indivs) > maxInds:
            remaining = maxInds - len(relateds)
            print currGen, relateds + random.sample(indivs, remaining)
        else:
            print currGen, relateds + indivs
