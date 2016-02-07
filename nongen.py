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
                if currGen <= 1: print "bonkers", currGen

            
        offspring = []
        prevMature = matureAge
        matureAge = {}
        currGen = stat["gen"]
    if stat["rep"] == 0 and stat["gen"]>0:
        if not years.has_key(stat["gen"]):
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
                

ages = list(set(filter(lambda x:float(x), years[maxYear]["age"])))
ages.sort()
for year in years:
    print >>demoout, "{0:3d}".format(year ),
    print >>demoMout, "{0:3d}".format(year ),
    print >>demoFout, "{0:3d}".format(year ),
    if len(ageParents[year]) > 0:
        print >>demoout, "{0:2.2f}".format(1.0*sum(ageParents[year])/len(ageParents[year])),
    else:
        print >>demoout, "00.00",

    for age in ages:
        print >>demoout, "{0:3d}".format(len(filter(lambda x: x==age, years[year]["age"]))),
        print >>demoMout, "{0:3d}".format(len(filter(lambda x: x==age, years[year]["mage"]))),
        print >>demoFout, "{0:3d}".format(len(filter(lambda x: x==age, years[year]["fage"]))),
    print >>demoout, ""
    print >>demoMout, ""
    print >>demoFout, ""


for gen in parentAgeGen:
    if gen==1: continue
    print >>breeout, "{0:3d}".format(gen),
    for age in ages:
        print >>breeout, "{0:3d}".format(len(filter(lambda x: x==age, parentAgeGen[gen]))),
    print >>breeout, ""

demoout.close()
breeout.close()
demoMout.close()
demoFout.close()
