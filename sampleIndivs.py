from sys import stdin, stderr, argv, exit
import random

if len(argv) not in [5, 6]:
    print "python %s condition indivs startGen endGen [lp]" % (argv[0],)
    print "indivs can be ALL"
    exit(-1)

cond = argv[1]
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

indivs = []
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
        indivs = []
        currGen=gen
        if currGen>endGen:
            break
    
    if eval(cond):
        indivs.append(id)
    l = stdin.readline()
if currGen<=endGen:
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
