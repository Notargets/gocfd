#!/bin/bash

#MESH=nacaAirfoil-fine.su2
MESH=naca12_2d-medium-lal.su2

for ((ii=0; ii<2; ii++))
do
	taskset -c 1 gocfd 2D -I input-bench.yaml -F mesh/$MESH -s 100
	taskset -c 1-2 gocfd 2D -I input-bench.yaml -F mesh/$MESH -s 100
	taskset -c 1-4 gocfd 2D -I input-bench.yaml -F mesh/$MESH -s 100
	taskset -c 1-8 gocfd 2D -I input-bench.yaml -F mesh/$MESH -s 100
	taskset -c 1-16 gocfd 2D -I input-bench.yaml -F mesh/$MESH -s 100
	taskset -c 0-31 gocfd 2D -I input-bench.yaml -F mesh/$MESH -s 100
	gocfd 2D -I input-bench.yaml -F mesh/$MESH -s 100
done
