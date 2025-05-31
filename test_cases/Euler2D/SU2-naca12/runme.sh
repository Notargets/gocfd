#!/bin/bash

#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2
for order in 0 1 2 3 4
#for order in 4
do
	cat input-base.yaml > tmp.yaml
	echo "PolynomialOrder: $order" >> tmp.yaml
	if [ $order -eq 0 ]; then
		echo "CFL: 2.0" >> tmp.yaml
		echo "MaxIterations: 6000" >> tmp.yaml
	elif [ $order -eq 1 ]; then
		echo "CFL: 2.0" >> tmp.yaml
		echo "MaxIterations: 15000" >> tmp.yaml
	elif [ $order -eq 2 ]; then
		echo "MaxIterations: 20000" >> tmp.yaml
	elif [ $order -eq 4 ]; then
		echo "CFL: 2.5" >> tmp.yaml
		echo "MaxIterations: 60000" >> tmp.yaml
	fi
	gocfd 2D -I tmp.yaml -F mesh_NACA0012_inv.su2 >& sysout-$order
    if [ -f meshfile.gobcfd ]; then
    	mv meshfile.gobcfd meshfile-$order.gobcfd
    	mv solutionfile.gobcfd solutionfile-$order.gobcfd
    	mv plotfile.dat plotfile-$order.dat
	fi
done
