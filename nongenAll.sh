REPS=$1
REPE=$2
MODEL=$3
DIR=$4
STARTGEN=$5
ENDGEN=$6
MATUREAGE=$7
DDIR=data/$DIR

for ((rep = $REPS; rep<$REPE; rep++))
do
  bzcat ${DDIR}/${MODEL}${rep}.sim.bz2 |python nongen.py ${DDIR}/${MODEL}${rep}.demo ${DDIR}/${MODEL}${rep}.bree
  bzcat ${DDIR}/${MODEL}${rep}.sim.bz2 |python ofs.py > ${DDIR}/${MODEL}${rep}.ofs
  cat ${DDIR}/${MODEL}${rep}.ofs |python vk.py ${DDIR}/${MODEL}${rep}.demo  $STARTGEN $ENDGEN > ${DDIR}/${MODEL}${rep}.vk
  cat ${DDIR}/${MODEL}${rep}.ofs |python vk.py ${DDIR}/${MODEL}${rep}.demo  $STARTGEN $ENDGEN $MATUREAGE > ${DDIR}/${MODEL}${rep}.vk.mature
done

