from __future__ import division, print_function

import os
import subprocess
import sys
import time

agene = '/home/tra/AgeNeV2.exe'
model = sys.argv[1]
N1 = int(sys.argv[2])
with_bsf = len(sys.argv) > 3

os.system('python generateAgeNe.py %s %d %s > tmp.agene' % (
    model, N1, 'yes' if with_bsf else ''))

os.remove('out.tmp')
proc = subprocess.Popen(['wine', agene],
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE)
proc.stdin.write('tmp.agene\r\nout.tmp\r\n\r\n')
proc.stdin.flush()
time.sleep(3)
f = open('out.tmp')
for l in f:
    if l.find('Ne (Eq 2) ') > -1:
        ne = float(l.rstrip().split(' ')[-1])
    if l.find('Nb ') > -1:
        # We want the last
        nb = float(l.rstrip().split(' ')[-1])
f.close()
print((nb, ne, nb / ne))
