#!/usr/bin/env python
import os
import sys

line = sys.argv[1]
for arg in sys.argv[2:]:
    if arg.find(' ') > -1:
        arg = "'" + arg + "'"
    line += ' ' + arg
print(line)
os.system(line)
