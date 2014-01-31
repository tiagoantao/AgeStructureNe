MODEL=$1
REPST=$2
REPED=$3
DIR=$4
MAXGEN=$5
STARTLOCUS=$6

DDIR=../data/$DIR
mkdir $$
cd $$
for ((rep=${REPST} ; rep<${REPED}; rep++))
do
    bzcat ${DDIR}/${MODEL}${rep}.sim.bz2 |python ../sampleIndivs.py True ALL 0 $MAXGEN|python ../sampleLoci.py ${DDIR}/${MODEL}${rep}.gen.bz2 100 100 $STARTLOCUS|python ../testHz.py 0 $MAXGEN > ${DDIR}/${MODEL}-${STARTLOCUS}-${rep}.hz
    #bzcat ${DDIR}/${MODEL}${rep}.sim.bz2 |python ../sampleIndivs.py age==1 ALL 0 $MAXGEN|python ../sampleLoci.py ${DDIR}/${MODEL}${rep}.gen.bz2 100 100 $STARTLOCUS|python ../testHz.py 0 $MAXGEN > ${DDIR}/${MODEL}-${STARTLOCUS}-${rep}.hz-newb
done
cd ..
rm -rf $$
