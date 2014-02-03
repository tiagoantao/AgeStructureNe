#!/bin/bash

QUEUE=put_queue_here
echo \#!/bin/bash > b.$$
echo PYTHONPATH=path_here b.$$
echo export PYTHONPATH >> b.$$
echo $@ >> b.$$
qsub -S /bin/bash -V -cwd -P $QUEUE -l h_vmem=2g b.$$
