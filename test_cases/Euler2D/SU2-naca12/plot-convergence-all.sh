#!/bin/bash

for order in 0 1 2 3 4
do
   	cat sysout-$order | ~/bin/plot_sysout.sh
	mv plot_output.png naca0012-conv-comparetoSU2-$order.png
done
