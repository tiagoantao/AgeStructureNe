AgeStructureNe
==============

Software to simulate and analyse Age Structured PopGen Data


Requirements
------------

Linux (64-bit)

Python 2.7

simuPOP

Biopython

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

Some of these calculations can be very slow or very computationally
demanding...

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

Generate summary statistics

python demo.py ldne/var-100grizzly.conf > demo.txt

python nb3.py ldne/var-100grizzly.conf > nb3.txt

These files can be loaded into Excel (space delimited)


Genetic Analysis
----------------

Heterozygosity

Do:

bash doH.sh 100grizzly 0 20 ldne 100 0

0 20 are the start and end replicates (0 to 19)

ldne is the directory name

100 is the number of generations

0 is the loci to start the analysis (0 - 100 are normally MSats, 100 onwards
are SNPs). 100 loci will be used to compute H.

This will generate hz computations on
data/ldne/1000grizzly-\*hz

The models on the directory hz simulate 600 years, probably more appropriate
to test Hz than the standard models (100 years)


Trout studies
-------------

[This assumes LDNe was run]

python nb3.py trout/var.conf > output/trt-nb3.txt 2> output/trt-nb.txt

Heterozygosity levels on the simulations need to be computed BEFORE the main
study is done. Do:

python sumTrout.py trout/var.conf trt

Errors starting with err can be ignored. There will be an exception at the
end.

Run:

python plotTrout.py trt

This will fail. Re-run both commands in order again:

python sumTrout.py trout/var.conf trt

python plotTrout.py trt


And now all should be fine (again, ignore the err messages). This is
convoluted but works!
