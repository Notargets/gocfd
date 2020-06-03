#!/bin/bash

rm -f sysout_* VAL_*

for ((i=2; i<=6; i++))
do
	echo "Running N = "$i
	for ((k=50;k<=200;k*=2))
	do
		#CASE="VAL_GK_LAX_N="$i
		#MODEL=2
		#gocfd -model $MODEL -CFL 1 -K $k -N $i -FinalTime 0.2 >> $CASE &
		CASE="VAL_DFR_ROE="$i
		MODEL=5
		gocfd -model $MODEL -CFL 0.5 -K $k -N $i -FinalTime 0.1 -Case 1 >> $CASE &
		CASE="VAL_DFR_LAX_N="$i
		MODEL=6
		gocfd -model $MODEL -CFL 0.5 -K $k -N $i -FinalTime 0.1 -Case 1 >> $CASE &
	done
	wait
done
