#!/bin/bash

for pts in 100pts 500pts
do
	filepre=shocktube-$pts-order
	for order in 0 1 2 3 4
	do
    	./plot.sh $filepre$order.dat
		mv plot_output.png $filepre$order.png
	done
done
