REPS=$1
REPE=$2
MODEL=$3
MATURE=$4

#./xrun.sh python go.py $REPS $REPE $MODEL All True 
./xrun.sh python go.py $REPS $REPE $MODEL Newb age==1 
#./xrun.sh python go.py $REPS $REPE $MODEL Mature \"$MATURE\" 

#./xrun.sh python go.py $REPS $REPE $MODEL c2c 2
#./xrun.sh python go.py $REPS $REPE $MODEL c3c 3
