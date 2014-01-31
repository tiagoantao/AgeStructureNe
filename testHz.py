from Bio.PopGen.GenePop.EasyController import EasyController
import sys
import os

if len(sys.argv) != 3:
    print "python %s startGen endGen" % (sys.argv[0],)
    sys.exit(-1)

startGen = int(sys.argv[1])
endGen = int(sys.argv[2])

fName = "/tmp/hz%s.txt" % (str(os.getpid()))

f = open(fName, "w")
l = sys.stdin.readline()
while l!="":
    f.write(l)
    l = sys.stdin.readline()
f.close()

ctrl = EasyController(fName)
pop_list, loci_list = ctrl.get_basic_info()
numPops = len(pop_list)
for pop in range(numPops):
    for locus in loci_list:
        #expho, obsho, exphe, obshe =  ctrl.get_heterozygosity_info(pop, locus)
        #print float(exphe)/(exphe+expho),
        obs, freqs = ctrl.get_allele_frequency(pop, locus)
        ho = 0.0
        for freq in freqs.values():
            ho += freq**2
        print 1-ho,
    print
os.remove(fName)
