from simuOpt import simuOptions
simuOptions["Quiet"] = True
from simuPOP import *
from myUtils import getConfig
import myUtils
from sys import argv, stderr
import math
import random
import numpy
import genepop
import os
#from mpimath import mp

if len(argv) not in [3, ]:
    print "Syntax:", argv[0], "<confFile> <prefout>"
    exit(-1)

cfg = getConfig(argv[1])
prefOut = argv[2]


def createGenome(size, numMSats, numSNPs):
    maxAlleleN = 100
    #print "Mutation model is most probably not correct", numMSats, numSNPs
    loci = (numMSats+numSNPs)*[1]
    initOps = []

    for msat in range(numMSats):
        diri = numpy.random.mtrand.dirichlet([1.0]*
                                             cfg.startAlleles)
        if type(diri[0]) == float:
            diriList = diri
        else:
            diriList = list(diri)

        initOps.append(
            InitGenotype(
                freq = [0.0] * ((maxAlleleN+1-8)//2) +
                    diriList +
                    [0.0] * ((maxAlleleN+1-8)//2),
                loci=msat
           )
        )

    for snp in range(numSNPs):
        freq = 0.5
        initOps.append(
            InitGenotype(
                #Position 0 is coded as 0, not good for genepop
                freq=[0.0, freq, 1 - freq],
                loci=numMSats + snp
           )
        )

    preOps = []
    if cfg.mutFreq > 0:
        preOps.append(StepwiseMutator(rates=cfg.mutFreq,
                      loci=range(numMSats)))
    return loci, initOps, preOps


def createSinglePop(popSize, nLoci, startLambda = 99999, lbd=1.0):
    initOps = [InitSex(maleFreq=cfg.maleProb)]
    if startLambda<99999:
        preOps = [ResizeSubPops(proportions=(float(lbd),), begin=startLambda),]
    else:
        preOps = []
    postOps = []
    pop = Population(popSize, ploidy=2,
                     loci = [1]*nLoci, chromTypes=[AUTOSOME]*nLoci,
                     infoFields=["ind_id", "father_id", "mother_id",
                                 "age", "breed", "rep_succ",
                                 "mate", "force_skip"])
    for ind in pop.individuals():
        ind.breed = -1000
    oExpr=('"%s/samp/%f/%%d/%%d/smp-%d-%%d-%%d.txt" %% ' +
            '(numIndivs, numLoci, gen, rep)')  % (
              cfg.dataDir, cfg.mutFreq, popSize)
    return pop, initOps, preOps, postOps, oExpr


def createSim(pop, reps):
    sim = Simulator (pop, rep = reps)
    return sim

def evolveSim(sim, gens, mateOp,
              genInitOps, genPreOps, popInitOps,
              ageInitOps, popPreOps, agePreOps, popPostOps, agePostOps,
              reportOps,
              oExpr):
    sim.evolve(
        initOps = genInitOps + popInitOps + ageInitOps,
        preOps = popPreOps + genPreOps + agePreOps,
        postOps = popPostOps + reportOps + agePostOps,
        matingScheme = mateOp,
        gen=gens)




(pop, popInitOps, popPreOps, popPostOps, oExpr) = createSinglePop(
    cfg.popSize, cfg.numMSats + cfg.numSNPs, cfg.startLambda, cfg.lbd)
(loci, genInitOps, genPreOps) = createGenome(cfg.popSize, cfg.numMSats, cfg.numSNPs)


def calcDemo(gen, pop):
    myAges = []
    for age in range(cfg.ages - 2):
        myAges.append(age+1)
    curr = 0
    for i in pop.individuals():
        if i.age in myAges:
            curr += 1
    if gen>=cfg.startLambda:
        cfg.N0 = cfg.N0*cfg.lbd
    return cfg.N0 + curr

def getRandomPos(arr):
    sumVal = sum(arr)
    rnd = random.random()
    acu = 0.0
    for i in range(len(arr)):
        acu += arr[i]
        if acu>=rnd*sumVal:
            return i

lSizes = [0,0,0,0,0,0]


def litterSkipGenerator(pop, subPop):
    fecms = cfg.fecundityMale
    fecfs = cfg.fecundityFemale
    nextFemales = []
    malesAge = {}
    femalesAge = {}
    availOfs = {}
    gen = pop.dvars().gen
    nLitter = None
    if cfg.litter and cfg.litter[0] < 0:
        nLitter = - cfg.litter[0]
    for ind in pop.individuals():
        if ind.sex() == 1: #male
            malesAge.setdefault(int(ind.age), []).append(ind)
        else:
            if nLitter is not None:
                availOfs[ind] = nLitter
            diff = int(gen - ind.breed)
            if diff > len(cfg.skip):
                available = True
                #print diff, len(cfg.skip), gen, ind.breed
            else:
                prob = random.random() * 100
                #print prob, cfg.skip, diff
                if prob > cfg.skip[diff - 1]:
                    available = True
                else:
                    available = False
            #print ind, available
            if available:
                femalesAge.setdefault(int(ind.age), []).append(ind)

    maleFec = []
    for i in range(len(fecms)):
        maleFec.append(fecms[i] * len(malesAge.get(i + 1, [])))
    femaleFec = []
    for i in range(len(fecfs)):
        if cfg.forceSkip > 0 and random.random() < cfg.forceSkip:
            femaleFec.append(0.0)
        else:
            femaleFec.append(fecfs[i] * len(femalesAge.get(i + 1, [])))

    while True:
        female = None
        if len(nextFemales) > 0:
            female = nextFemales.pop()
        while not female:
            age = getRandomPos(femaleFec) + 1
            if len(femalesAge.get(age, [])) > 0:
                female = random.choice(femalesAge[age])
                if nLitter is not None:
                    if availOfs[female] == 0:
                        female = None
                    else:
                        availOfs[female] -= 1
                elif cfg.litter:
                    lSize = getRandomPos(cfg.litter) + 1
                    lSizes[lSize] += 1
                    if lSize > 1:
                        nextFemales = [female] * (lSize - 1)
                    femalesAge[age].remove(female)

        male = None
        if cfg.isMonog:
            if female.mate > -1:
                male = pop.indByID(female.mate)
        while male is None:
            age = getRandomPos(maleFec) + 1
            if len(malesAge.get(age, [])) > 0:
                male = random.choice(malesAge[age])
            if cfg.isMonog:
                if male.mate > -1:
                    male = None
                else:
                    male.mate = female.ind_id

        female.breed = gen
        if cfg.isMonog:
            female.mate = male.ind_id
        yield male, female


def fitnessGenerator(pop, subPop):
    fecms = cfg.fecundityMale
    fecfs = cfg.fecundityFemale
    totFecMales = 0.0
    totFecFemales = 0.0
    availableFemales = []
    perAgeMaleNorm = {}
    perAgeFemaleNorm = {}
    gen = pop.dvars().gen
    ageCntMale = {}
    ageCntFemale = {}
    for ind in pop.individuals():
        if ind.sex() == 1: #male
            a = cfg.gammaAMale[int(ind.age)-1]
            b = cfg.gammaBMale[int(ind.age)-1]
            if a:
                gamma = numpy.random.gamma(a, b)
                ind.rep_succ = gamma
                #ind.rep_succ = numpy.random.poisson(gamma)
            else:
                ind.rep_succ = 1
            perAgeMaleNorm[int(ind.age)-1] = perAgeMaleNorm.get(int(ind.age)-1, 0.0) + ind.rep_succ
            ageCntMale[int(ind.age)-1] = ageCntMale.get(int(ind.age)-1, 0.0) + 1
        else:
            #if ind.age == 0: totFecFemales +=0
            a = cfg.gammaAFemale[int(ind.age)-1]
            b = cfg.gammaBFemale[int(ind.age)-1]
            if a:
                gamma = numpy.random.gamma(a, b)
                ind.rep_succ = gamma
                #ind.rep_succ = numpy.random.poisson(gamma)
            else:
                ind.rep_succ = 1
            perAgeFemaleNorm[int(ind.age)-1] = perAgeFemaleNorm.get(int(ind.age)-1, 0.0) + ind.rep_succ
            ageCntFemale[int(ind.age)-1] = ageCntFemale.get(int(ind.age)-1, 0.0) + 1
            availableFemales.append(ind.ind_id)

    for ind in pop.individuals():
        if ind.sex() == 1: #male
            if perAgeMaleNorm[int(ind.age)-1]==0.0:
                ind.rep_succ = 0.0
            else:
                ind.rep_succ = ageCntMale[int(ind.age)-1]*fecms[int(ind.age)-1] *ind.rep_succ/perAgeMaleNorm[int(ind.age)-1]
            totFecMales += ind.rep_succ
        else:
            if ind.ind_id not in availableFemales:
                continue
            if perAgeFemaleNorm[int(ind.age)-1]==0.0:
                ind.rep_succ = 0.0
            else:
                ind.rep_succ =  ageCntFemale[int(ind.age)-1]*fecfs[int(ind.age)-1]*ind.rep_succ/perAgeFemaleNorm[int(ind.age)-1]
            totFecFemales += ind.rep_succ

    nextFemales = []
    while True:
        mVal = random.random()*totFecMales
        fVal = random.random()*totFecFemales
        runMale = 0.0
        runFemale = 0.0
        male = False
        female = False
        if len(nextFemales) > 0:
            female = nextFemales.pop()
            female.breed = gen
        inds = list(pop.individuals())
        random.shuffle(inds)
        for ind in inds:
            if ind.age == 0: continue
            if ind.sex() == 1 and not male: #male
                runMale += ind.rep_succ
                if runMale>mVal:
                    male = ind
            elif ind.sex() == 2 and not female:
                if ind.ind_id not in availableFemales: continue
                runFemale += ind.rep_succ
                if runFemale>fVal:
                    female = ind
                    female.breed = gen
            if male and female:
                break
        yield male, female

def cull(pop):
    kills = []
    for i in pop.individuals():
        if i.age>0 and i.age<cfg.ages-1:
            if i.sex()==1:
                cut = cfg.survivalMale[int(i.age)-1]
            else:
                cut = cfg.survivalFemale[int(i.age)-1]
            if random.random()>cut:
                kills.append(i.ind_id)
    pop.removeIndividuals(IDs=kills)
    return True

def zeroC(v):
    a = str(v)
    while len(a)<3:
        a = "0" + a
    return a

def outputAge(pop):
    gen = pop.dvars().gen
    if gen<cfg.startSave: return True
    rep = pop.dvars().rep
    for i in pop.individuals():
        out.write("%d %d %d %d %d %d %d\n" % (gen, rep, i.ind_id, i.sex(),
            i.father_id, i.mother_id, i.age))
        if i.age == 1 or gen==0:
            err.write("%d %d "%  (i.ind_id, gen))
            for pos in range(len(i.genotype(0))):
                a1 = zeroC(i.allele(pos, 0))
                a2 = zeroC(i.allele(pos, 1))
                err.write(a1+a2+" ")
            err.write("\n")
    return True

def outputMega(pop):
    gen = pop.dvars().gen
    if gen<cfg.startSave: return True
    for i in pop.individuals():
        if i.age==0:
            megaDB.write( "%d %d %d %d %d\n" % (gen, i.ind_id, i.sex(),
                i.father_id, i.mother_id))
    return True

def setAge(pop):
    probMale = [1.0]
    for sv in cfg.survivalMale:
        probMale.append(probMale[-1]*sv)
    totMale = sum(probMale)
    probFemale = [1.0]
    for sv in cfg.survivalFemale:
        probFemale.append(probFemale[-1]*sv)
    totFemale = sum(probFemale)
    for ind in pop.individuals():
        if ind.sex()==1:
            prob = probMale
            tot = totMale
        else:
            prob = probFemale
            tot = totFemale
        cut = tot*random.random()
        acu = 0
        for i in range(len(prob)):
            acu += prob[i]
            if acu>cut:
                age = i
                break
        ind.age = age

    return True


def createAge(pop):
    ageInitOps = [
        #InitInfo(lambda: random.randint(0, cfg.ages-2), infoFields='age'),
        IdTagger(),
        #PyOperator(func=outputAge,at=[0]),
        PyOperator(func=setAge,at=[0]),
    ]
    agePreOps = [
        InfoExec("age += 1"),
        InfoExec("mate = -1"),
        InfoExec("force_skip = 0"),
        PyOperator(func=outputAge),
    ]
    mySubPops = []
    for age in range(cfg.ages - 2):
        mySubPops.append((0,age+1))
    mateOp = HeteroMating([
            HomoMating(
                    PyParentsChooser(fitnessGenerator if cfg.doNegBinom
                        else litterSkipGenerator),
                    OffspringGenerator(numOffspring=1,
                        ops = [
                            MendelianGenoTransmitter(),
                            IdTagger(),
                            PedigreeTagger(),
                        ], sexMode = (PROB_OF_MALES, cfg.maleProb)),
                    weight=1
                ),
            CloneMating(subPops=mySubPops,weight=-1),
            ],

            subPopSize=calcDemo
        )
    agePostOps = [
        PyOperator(func=outputMega),
        PyOperator(func=cull),
    ]
    pop.setVirtualSplitter(InfoSplitter(field='age', cutoff=
       range(1,cfg.ages)))
    return ageInitOps, agePreOps, mateOp, agePostOps

(ageInitOps, agePreOps, mateOp, agePostOps) = createAge(pop)


out=open(prefOut + ".sim" ,"w")
err=open(prefOut + ".gen" ,"w")
megaDB=open(prefOut + ".db" ,"w")
reportOps = [
    Stat(popSize=True),
    #PyEval(r'"gen %d\n" % gen', reps=0),
    #PyEval(r'"size %s\n" % subPopSize', reps=0),
]
sim = createSim(pop, cfg.reps)
evolveSim(sim, cfg.gens, mateOp,
              genInitOps,
              genPreOps, popInitOps, ageInitOps, popPreOps, agePreOps,
              popPostOps, agePostOps, reportOps,
              oExpr)
out.close()
err.close()
megaDB.close()
print lSizes
