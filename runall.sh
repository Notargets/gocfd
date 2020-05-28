#!/bin/bash

rm -f sysout_*

for ((i=2; i<3; i++))
do
	echo "Running N = "$i
	gocfd -model 6 -CFL 1 -K 100 -N $i -FinalTime 0.2 >> sysout_$i
	gocfd -model 6 -CFL 1 -K 200 -N $i -FinalTime 0.2 >> sysout_$i
	gocfd -model 6 -CFL 1 -K 400 -N $i -FinalTime 0.2 >> sysout_$i
	gocfd -model 6 -CFL 1 -K 800 -N $i -FinalTime 0.2 >> sysout_$i
	gocfd -model 6 -CFL 1 -K 1600 -N $i -FinalTime 0.2 >> sysout_$i
	gocfd -model 6 -CFL 1 -K 3200 -N $i -FinalTime 0.2 >> sysout_$i
done
