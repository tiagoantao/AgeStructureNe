import os
import sys

reps = int(sys.argv[1])
repe = int(sys.argv[2])
model = sys.argv[3]
agedesc = sys.argv[4]
agecond = sys.argv[5]
#ex 0 20 250ecology All True
#ex 0 20 250ecology Newb age==1
#ex 0 20 250ecology Mature age\>1
DDIR = "../data/trout"
#N1modelAgeIndivsLoci-rep

gens = 1000

sampleStratsLoci = False, [(10, 50), (15, 50), (25, 50), (50, 50), (100, 50)]
sampleStratsLoci = False, [(15, 50), (50, 50), (100, 50)]
sampleStratsLociSNP = True, [(100, 50), (200, 50), (400, 50)]
sampleStratsLociSNP = True, []
sampleStratsIndivs = False, [(15, 15), (15, 25), (15, 50), (15, 100)]
sampleStratsIndivs = False, []
sampleStratsIndivsSNP = True, [(100, 15), (100, 25), (100, 50), (100, 100)]
sampleStratsIndivsSNP = True, []

myd = {"DDIR": DDIR, "MODEL": model,
       "AGECOND": agecond, "AGEDESC": agedesc, "GENS": gens}
print model, agecond

for rep in range(reps, repe):
    print "REP", rep
    myd["rep"] = rep
    if agedesc in ["c2c", "c3c"]:
        if model not in ["180bulltrout", "361bulltrout"]:
            continue
        myd["thres"] = 0.011
        myd["nindivs"] = 50
        myd["nloci"] = 15
        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivsCohort.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ld2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep}.r2'.format(**myd))
        myd["nloci"] = 100
        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivsCohort.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ld2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep}.r2'.format(**myd))
    elif agedesc in ["Mature", "All"]:
        if model not in ["180bulltrout", "361bulltrout", "722bulltrout",
                         "518shepard", "193bullpred", "775bullpred", "1619btrout",
                         "6476btrout", "641fraley", "3050bullt2",
                         "18lake", "72lake"]:
            continue
        myd["thres"] = 0.011
        myd["nindivs"] = 50
        myd["nloci"] = 15
        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ld2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep}.r2'.format(**myd))
        myd["nloci"] = 100
        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ld2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep}.r2'.format(**myd))
    else:
        for isSNP, strat in [sampleStratsLociSNP, sampleStratsIndivsSNP, sampleStratsLoci, sampleStratsIndivs]:
            isMSAT = not isSNP
            for nloci, nindivs in strat:
                print nloci, nindivs, isSNP
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
                    if model not in ["180bulltrout", "361bulltrout"]:
                        continue
                    os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ld2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-{rep}.r2'.format(**myd))
                else:
                    #if model not in ["180bulltrout", "361bulltrout"]:
                    #    if nindivs != 50 or nloci != 100:
                    #        continue
                    os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ld2.py {thres} > ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep} 2> ldout/{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep}.r2'.format(**myd))
                # related individuals
                #if nindivs == 50:
                if False:
                    if model not in ["180bulltrout", "361bulltrout"]:
                        continue
                    if isMSAT and nloci != 15:
                        continue
                    if isMSAT:
                        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivsRelated.py 20 {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ld2.py {thres} > ldout/{MODEL}All{nindivs}{nloci}-rel-{rep} 2> ldout/{MODEL}All{nindivs}{nloci}-rel-{rep}.r2'.format(**myd))
                    else:
                        os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivsRelated.py 20 {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ld2.py {thres} > ldout/{MODEL}All{nindivs}{nloci}-snp-rel-{rep} 2> ldout/{MODEL}All{nindivs}{nloci}-snp-rel-{rep}.r2'.format(**myd))
                # Thresholds
                #if nindivs == 50 and isSNP:
                if True:
                    if model not in ["180bulltrout", "361bulltrout"]:
                        continue
                    for thres in [0.021, 0.035, 0.05, 0.1]:
                        myd["thres"] = thres
                        if isSNP:
                            os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 400 100|python ../ld2.py {thres} > ldout/{thres}-{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep} 2> ldout/{thres}-{MODEL}{AGEDESC}{nindivs}{nloci}-snp-{rep}.r2'.format(**myd))
                        else:
                            os.system('bzcat {DDIR}/{MODEL}{rep}.sim.bz2 |python ../sampleIndivs.py "{AGECOND}" {nindivs} 1 {GENS}|python ../sampleLoci.py {DDIR}/{MODEL}{rep}.gen.bz2 {nloci} 100|python ../ld2.py {thres} > ldout/{thres}-{MODEL}{AGEDESC}{nindivs}{nloci}-{rep} 2> ldout/{thres}-{MODEL}{AGEDESC}{nindivs}{nloci}-{rep}.r2'.format(**myd))
