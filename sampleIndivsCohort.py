from sys import stdin, stderr, argv, exit
import random

if len(argv) not in [5, 6]:
    print "python %s maxAge indivs startGen endGen [lp]" % (argv[0],)
    print "indivs can be ALL"
    exit(-1)

maxAge = int(argv[1])
if argv[2]=="ALL":
    nindivs=None # Actually we stop at 200
    maxInds = 200
else:
    nindivs = int(argv[2])
startGen = int(argv[3])
endGen = int(argv[4])
lp = len(argv) == 6

l = stdin.readline()

#currGen = cfg.startSave
currGen = 0

def getIndivs(ageInds):
    indivs = []
    maxIndsPerCohort = min(map(lambda x:len(x), ageInds.values()))
    for age, inds in ageInds.items():
        indivs.extend(random.sample(inds, maxIndsPerCohort))
    return indivs

ageCnt={}
indivsAge = {}
while l!= "":
    toks = l.rstrip().split(" ")
    rep = int(float(toks[1]))
    if rep>0: break
    gen = int(float(toks[0]))
    id = int(float(toks[2]))
    sex = int(float(toks[3]))
    father = int(float(toks[4]))
    mother = int(float(toks[5]))
    age = int(float(toks[6]))

    if gen>currGen:
        if currGen>=startGen:
            indivs = getIndivs(indivsAge)
            if nindivs:
                print currGen, random.sample(indivs, nindivs)
            else:
                if len(indivs)>maxInds:
                    print currGen, random.sample(indivs, maxInds)
                else:
                    print currGen, indivs
            if lp:
                if nindivs:
                    print random.sample(indivs, nindivs)
                else:
                    if len(indivs)>maxInds:
                        print random.sample(indivs, maxInds)
                    else:
                        print indivs
        indivsAge = {}
        currGen=gen
        if currGen>endGen:
            break
    
    if age<=maxAge:
        ageCnt[age] = ageCnt.get(age, 0) + 1
        indivsAge.setdefault(age,[]).append(id)
    l = stdin.readline()
if currGen<=endGen:
    indivs = getIndivs(indivsAge)
    if nindivs:
        print currGen, random.sample(indivs, nindivs)
    else:
        if len(indivs)>maxInds:
            print currGen, random.sample(indivs, maxInds)
        else:
            print currGen, indivs
    if lp:
        if nindivs:
            print random.sample(indivs, nindivs)
        else:
            if len(indivs)>maxInds:
                print random.sample(indivs, maxInds)
            else:
                print indivs
