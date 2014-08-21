import sys

from myUtils import ages, survivalFemale, survivalMale, fecFemale, fecMale

model = sys.argv[1]
N1 = int(sys.argv[2])

print(model)
print('%d %d 0.5' % (ages[model] - 1, N1))
survivalFemale[model].append(0)
survivalMale[model].append(0)
for age in range(ages[model] - 1):
    print('%d %f %f 1 %f %f 1' %
          (age + 1, survivalFemale[model][age], fecFemale[model][age],
           survivalMale[model][age], fecMale[model][age]))
