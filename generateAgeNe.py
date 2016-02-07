from __future__ import division
import sys

from myUtils import ages, survivalFemale, survivalMale, fecFemale, fecMale

model = sys.argv[1]
N1 = int(sys.argv[2])
with_bsf = len(sys.argv) > 3

print(model)
print('{0:d} {1:d} 0.5'.format(ages[model] - 1, N1))
survivalFemale[model].append(0)
survivalMale[model].append(0)
curr_survival = N1 / 2, N1 / 2
for age in range(ages[model] - 1):
    if with_bsf:
        if curr_survival[0] < 1:
            bf = 0
        else:
            bf = (curr_survival[0] - 1) / curr_survival[0]
        if curr_survival[1] < 1:
            bm = 0
        else:
            bm = (curr_survival[1] - 1) / curr_survival[1]
        curr_survival = curr_survival[0] * survivalFemale[model][age], curr_survival[1] * survivalMale[model][age]
    else:
        bf = bm = 1
    print('{0:d} {1:f} {2:f} {3:f} {4:f} {5:f} {6:f}'.format(age + 1, survivalFemale[model][age], fecFemale[model][age], bf,
           survivalMale[model][age], fecMale[model][age], bm))
