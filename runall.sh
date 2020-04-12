#!/bin/bash

for ((i=1; i<5; i++))
do
	echo $i
	gocfd -K 8 -N $i -delay 150 -graph >> sysout &
done
