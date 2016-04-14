from __future__ import print_function
from sys import stdin, argv, exit
import scipy

if len(argv) not in [4, 5, 6]:
    print("python %s demo startGen endGen matureAgeM matureAgeF" % (argv[0]))
    exit(-1)

# stdin is ofs
demo = argv[1]  # from nogen.py
startGen = int(argv[2])
endGen = int(argv[3])

if len(argv) == 6:
    matureAgeM = int(argv[4])
    matureAgeF = int(argv[5])
else:
    matureAgeM = None
    matureAgeF = None

L = {}
N1s = {}
popSizes = {}

f = open(demo)
l = f.readline()
while l != "":
    toks = [float(x) for x in [y for y in l.rstrip().split(" ") if y != ""]]
    toks[0] = int(toks[0])
    popSizes[toks[0]] = sum(toks[2:])
    N1s[toks[0]] = toks[2]
    L[toks[0]] = toks[1]
    l = f.readline()
f.close()

popSizesM = {}
popSizesF = {}


def getGenderDemo(name, sizes):
    f = open(name)
    l = f.readline()
    while l != "":
        toks = [float(x) for x in [y for y in l.rstrip().split(" ")
                                   if y != ""]]
        toks[0] = int(toks[0])
        sizes[toks[0]] = sum(toks[1:])
        l = f.readline()
    f.close()
getGenderDemo(demo + "M", popSizesM)
getGenderDemo(demo + "F", popSizesF)

l = stdin.readline()
ofGens = {}
myOfs = {}
currGen = startGen
while l != "":
    toks = l.rstrip().split(" ")
    myGen = int(toks[0])
    if myGen < startGen:
        l = stdin.readline()
        continue
    if myGen > endGen:
        break
    if myGen > currGen:
        ofGens[currGen] = myOfs
        myOfs = {}
        currGen = myGen
    bornGen = int(toks[1])
    id = int(toks[2])
    gender = int(toks[3])
    nofs = int(toks[4])
    age = myGen - bornGen - 1
    myOfs.setdefault(age, []).append(nofs)
    l = stdin.readline()

NeHills = []
for gen in ofGens:
    tot = []
    ofs = ofGens[gen]
    for age in ofs:
        if len(ofs[age]) == 0:
            continue
        a = scipy.array(ofs[age])
        tot.extend(ofs[age])
        # print age, a.mean(), a.var()
    tota = scipy.array(tot)
    k, Vk = tota.mean(), tota.var()
    NeHill = 4 * N1s[gen] * L[gen] / (Vk + 2)
    NeHills.append(NeHill)
    print(gen, k, Vk, (k*popSizes[gen] - 2) / (k-1 + (Vk/k)), NeHill)
