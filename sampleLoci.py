from __future__ import print_function
from sys import stdin, argv, exit
import bz2
import random

if len(argv) not in [4, 5]:
    print("python %s genFile loci maxLoci [startLoci]" % (argv[0],))
    exit(-1)

genFile = argv[1]
nloci = int(argv[2])
maxLoci = int(argv[3])
startLoci = int(argv[4]) if len(argv) == 5 else 0

loci = [x + startLoci for x in random.sample(list(range(maxLoci)), nloci)]

l = stdin.readline()
gens = []

currGen = 0

# get individuals per generation
indivs = set()
while l != "":
    l = l.rstrip()
    point = l.find(" ")
    gen = l[:point]
    genIndivs = eval(l[point:])
    indivs = indivs.union(genIndivs)
    gens.append(genIndivs)
    l = stdin.readline()

# get genetic data
f = bz2.open(genFile, 'rt')
l = f.readline()
genetics = {}
while l != "":
    toks = l.rstrip().split(" ")
    id = int(float(toks[0]))
    gen = int(float(toks[1]))
    myAlleles = toks[2:]
    if id in indivs:
        myloci = []
        for locus in loci:
            myloci.append(myAlleles[locus])
        genetics[id] = myloci
    l = f.readline()
f.close()

# print >>stderr, genetics.keys()

# dump genepop file
print("lala land")
for locus in loci:
    print("l" + str(locus))
for indivs in gens:
    print("Pop")
    for indiv in indivs:
        print("i" + str(indiv) + ",", end=' ')
        print(" ".join(genetics[indiv]))
