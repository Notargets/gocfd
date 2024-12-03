#!/bin/bash

for p in 1 2 4 8 16 32
do
	gocfd 2D -I input-bench.yaml -F mesh/naca12_2d-medium-lal.su2 -s 100 -p $p
done

for p in 1 2 4 8 16 32
do
	gocfd 2D -I input-bench.yaml -F mesh/nacaAirfoil-fine.su2 -s 100 -p $p
done
