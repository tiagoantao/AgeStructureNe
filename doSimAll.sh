M=$1
project=$2
echo $M
bash doSim.sh $M  0  2 $project &
bash doSim.sh $M  2  4 $project &
bash doSim.sh $M  4  6 $project &
bash doSim.sh $M  6  8 $project &
bash doSim.sh $M  8 10 $project &
bash doSim.sh $M 10 12 $project &
bash doSim.sh $M 12 14 $project &
bash doSim.sh $M 14 16 $project &
bash doSim.sh $M 16 18 $project &
bash doSim.sh $M 18 20 $project &
