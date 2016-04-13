from __future__ import print_function
from myUtils import getNonGenStats
from sys import stdin, stderr

stats = getNonGenStats(stdin)

gens={}
currGen = 0
alive = []
aliveNow = []
offs = {}
offsGen = {}
istats = {}
for stat in stats:
    if stat["gen"]>currGen:
        currGen = int(float(stat["gen"]))
        rems = []
        for dead in alive:
            if dead not in aliveNow:
                print(currGen, istats[dead][0], dead, istats[dead][1], offs[dead], offsGen.get(dead,0))
                del offs[dead]
                del istats[dead]
                rems.append(dead)
        for rem in rems:
            alive.remove(rem)
        aliveNow = []
        offsGen = {}
    id = stat["id"]
    gender = stat["gender"]
    age = stat["age"]
    father = stat["father"]
    mother = stat["mother"]
    aliveNow.append(id)
    if age == 1:
        alive.append(id)
        istats[id] = currGen, gender
        offs[id] = 0
        if father>0:
            offs[father] = offs.get(father,0)+1
            offsGen[father] = offsGen.get(father,0)+1
            offs[mother] = offs.get(mother,0)+1
            offsGen[mother] = offsGen.get(mother,0)+1
