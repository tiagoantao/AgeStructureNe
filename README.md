AgeStructureNe
==============

Software to simulate and analyse Age Structured PopGen Data


Requirements
------------

Linux (64-bit)

Python 2.7

simuPOP

Genepop (1)

LDNe (1,2)


1 - Only if you do genetic analysis

2 - Binary is supplied here


Instructions
------------

There are 3 parts to the analysis

1. Data simulation
2. Non-genetic analysis
3. Genetic analysis

You might want to skip 2 or 3 (e.g to repeat the Evolution paper, you do not
need genetic analysis)


Data Simulation
---------------


You can start by running the simulator on one scenario, for example

python sim.py ldne/100grizzly.conf example

Will run an age structured population of Grizzlies with a N1 of 100 and write
the output on example.* files

Non-genetic Analysis
--------------------

Genetic Analysis
----------------
