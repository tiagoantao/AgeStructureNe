MODEL=$1
DIR=$2
MAXGEN=$3
STARTLOCUS=$4
nohup bash doH.sh $MODEL 0  2 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 2  4 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 4  6 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 6  8 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 8  10 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 10 12 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 12 14 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 14 16 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 16 18 $DIR $MAXGEN $STARTLOCUS &
nohup bash doH.sh $MODEL 18 20 $DIR $MAXGEN $STARTLOCUS &
