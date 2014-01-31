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


You can start by running the simulator on one scenario, for example:

python sim.py ldne/100grizzly.conf example

Will run an age structured population of Grizzlies with a N1 of 100 and write
the output on example.\* files

This will generate 2 files:

example.sim - non-genetic output (generation, replicate, individual id sex and
parent id)

example.gen - the genome of each individual

Typically you will want to run several replicates. For that you can run

bash doSim.sh 100grizzly 0 20 ldne

This will do 20 replicates (0 to 19) and store them on data/ldne/100grizzly\* .
Ignore the rm errors.

If you desire you can run replicates in parallel (i.e. run several instances
of the above, see e.g doSimAll.sh which runs 20 replicates - 10 at a time).
You will need to have enough CPU power and memory for this!

Non-genetic Analysis
--------------------

Generate non-genetic raw information (vk, offspring count, ...)

python totalAll.py ldne/var-100grizzly.conf

(see ldne/var\*conf for examples)

Genetic Analysis
----------------
