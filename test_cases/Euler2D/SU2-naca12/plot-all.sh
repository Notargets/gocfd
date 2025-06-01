#!/bin/bash

for order in 0 1 2 3 4
do
   	./plot-cp.sh plotfile-$order.dat
	mv plot_output.png naca0012-comparetoSU2-$order.png
done
