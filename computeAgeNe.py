from __future__ import division

import os
import subprocess
import sys

agene = '/home/tra/AgeNeV2.exe'
model = sys.argv[1]
N1 = int(sys.argv[2])
with_bsf = len(sys.argv) > 3

os.system('python generateAgeNe.py %s %d %s > tmp.agene' % (
    model, N1, 'yes' if with_bsf else ''))

proc = subprocess.Popen(['wine', agene],
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE)
proc.stdin.write('tmp.agene\r\nout.tmp\r\n\r\n')
ls = proc.stdout.readlines()
for l in ls[-2:]:
    l = l.rstrip()
    print(l)
ne = float(ls[-2].rstrip().split(' ')[-1])
nb = float(ls[-1].rstrip().split(' ')[-1])
print (nb / ne)
