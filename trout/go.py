import os
import sys

reps = int(sys.argv[1])
repe = int(sys.argv[2])
model = sys.argv[3]
agedesc = sys.argv[4]
agecond = sys.argv[5]
# ex 0 20 250ecology All True
# ex 0 20 250ecology Newb age==1
# ex 0 20 250ecology Mature age\>1
DDIR = "../data/trout"
# N1modelAgeIndivsLoci-rep

if model.endswith('mosquito') or model.endswith('weed'):
    gens = 5000
else:
    gens = 1000

# loci, indivs
sampleStratsMSAT = False, []
sampleStratsSNP = True, [(100, 50)]
#fig2
#sampleStratsMSAT = False, [(15, 15), (15, 25), (15, 50), (15, 100),
#                           (100, 15), (100, 25), (100, 50), (100, 100)]
#sampleStratsIndivsSNP = True, [(100, 15), (100, 25), (100, 50), (100, 100),
#                               (200, 15), (200, 25), (200, 50), (200, 100),
#                               (400, 15), (400, 25), (400, 50), (400, 100)]

myd = {"DDIR": DDIR, "MODEL": model,
       "AGECOND": agecond, "AGEDESC": agedesc, "GENS": gens}
print(model, agecond)

for rep in range(reps, repe):
    print("REP", rep)
    myd["rep"] = rep
    if agedesc in ["c2c", "c3c"]:
        myd["thres"] = 0.011
        myd["nindivs"] = 50
        myd["nloci"] = 15
        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivsCohort.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ne2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep}.r2'.format(**myd))
        myd["nloci"] = 100
        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivsCohort.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ne2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep}.r2'.format(**myd))
    elif agedesc in ["Mature", "All"]:
        myd["thres"] = 0.011
        myd["nindivs"] = 50
        myd["nloci"] = 15
        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ne2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep}.r2'.format(**myd))
        myd["nloci"] = 100
        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ne2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep}.r2'.format(**myd))
    else:
        for isSNP, strat in [sampleStratsSNP, sampleStratsMSAT]:
            isMSAT = not isSNP
            for nloci, nindivs in strat:
                print(nloci, nindivs, isSNP)
                if nindivs <= 15:
                    thres = 0.035
                elif nindivs <= 25:
                    thres = 0.021
                else:
                    thres = 0.011
                myd["thres"] = thres
                myd["nindivs"] = nindivs
                myd["nloci"] = nloci
                if isMSAT:
                    os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ne2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep}.r2'.format(**myd))
                else:
                    os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ne2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep}.r2'.format(**myd))
                # related individuals
                #if nindivs == 50:
                if False:
                    if isMSAT:
                        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivsRelated.py 20 {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ne2.py {thres} > ldout/{MODEL}All{nindivs}{nloci}-rel-{rep} 2> ldout/{MODEL}All{nindivs}{nloci}-rel-{rep}.r2'.format(**myd))
                    else:
                        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivsRelated.py 20 {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ne2.py {thres} > ldout/{MODEL}All{nindivs}{nloci}-snp-rel-{rep} 2> ldout/{MODEL}All{nindivs}{nloci}-snp-rel-{rep}.r2'.format(**myd))
                # pcrit Thresholds
                #if nindivs == 50 and isSNP:
                if False:
                    for thres in [0.021, 0.035, 0.05, 0.1]:
                        myd["thres"] = thres
                        if isSNP:
                            os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ne2.py {thres} > ldout/{thres}-{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep} 2> ldout/{thres}-{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep}.r2'.format(**myd))
                        else:
                            os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ne2.py {thres} > ldout/{thres}-{MODEL}{AGEDESC}{nindivs}{nloci}-{rep} 2> ldout/{thres}-{MODEL}{AGEDESC}{nindivs}{nloci}-{rep}.r2'.format(**myd))
