from __future__ import print_function
from myUtils import getNonGenStats
from sys import stdin,argv

dout = argv[1]
bout = argv[2]
dMout = dout + "M"
dFout = dout + "F"

stats = getNonGenStats(stdin)

demoout=open(dout,"w")
breeout=open(bout,"w")
demoMout=open(dMout,"w")
demoFout=open(dFout,"w")

years={}
currGen = 0
prevMature = {}
matureAge = {}
parentAgeGen = {}
offspring = []
ageParents = {}
maxYear = 0
for stat in stats:
    if stat["gen"]>currGen:
        currGen = stat["gen"]
        parentAgeGen[currGen] = []
        ageParents[currGen]=[]
        for father, mother in offspring:
            try:
                parentAgeGen[currGen].append(prevMature[father])
                parentAgeGen[currGen].append(prevMature[mother])
                ageParents[currGen].extend(
                    [prevMature[father],prevMature[mother]])
            except KeyError:
                if currGen <= 1: print("bonkers", currGen)

            
        offspring = []
        prevMature = matureAge
        matureAge = {}
        currGen = stat["gen"]
    if stat["rep"] == 0 and stat["gen"]>0:
        if stat["gen"] not in years:
            years[stat["gen"]] = {}
            years[stat["gen"]][ "age"] = []
            years[stat["gen"]]["fage"] = []
            years[stat["gen"]]["mage"] = []
        di = years[stat["gen"]]
        maxYear = max((maxYear, stat["gen"]))
        age = stat["age"]
        di["age"].append(age)
        if stat["gender"] == 1:
            di["mage"].append(age)
        else:
            di["fage"].append(age)
        if age == 1:
            offspring.append((stat["father"], stat["mother"]))
        matureAge[stat["id"]] = age
                

ages = list(set([x for x in years[maxYear]["age"] if float(x)]))
ages.sort()
for year in years:
    print("%3d" %(year, ), end=' ', file=demoout)
    print("%3d" %(year, ), end=' ', file=demoMout)
    print("%3d" %(year, ), end=' ', file=demoFout)
    if len(ageParents[year]) > 0:
        print("%2.2f" % (1.0*sum(ageParents[year])/len(ageParents[year]),), end=' ', file=demoout)
    else:
        print("00.00", end=' ', file=demoout)

    for age in ages:
        print("%3d" % (len([x for x in years[year]["age"] if x==age]),), end=' ', file=demoout)
        print("%3d" % (len([x for x in years[year]["mage"] if x==age]),), end=' ', file=demoMout)
        print("%3d" % (len([x for x in years[year]["fage"] if x==age]),), end=' ', file=demoFout)
    print("", file=demoout)
    print("", file=demoMout)
    print("", file=demoFout)


for gen in parentAgeGen:
    if gen==1: continue
    print("%3d"% (gen,), end=' ', file=breeout)
    for age in ages:
        print("%3d" % (len([x for x in parentAgeGen[gen] if x==age]),), end=' ', file=breeout)
    print("", file=breeout)

demoout.close()
breeout.close()
demoMout.close()
demoFout.close()
