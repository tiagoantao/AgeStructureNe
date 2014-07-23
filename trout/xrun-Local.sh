#!/usr/bin/env python
import os
import sys

line = sys.argv[0]
for arg in sys.argv[1:]:
    if arg.find(' ') > -1:
        arg = "'" + arg + "'"
    line += ' ' + arg
os.system(line)
