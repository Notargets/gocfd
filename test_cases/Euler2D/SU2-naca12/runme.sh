#!/bin/bash

#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2
#for order in 0 1 2 3 4
for order in 0
do
	cat input-base.yaml > tmp.yaml
	echo "PolynomialOrder: $order" >> tmp.yaml
	if [ $order -eq 0 ]; then
		echo "CFL: 2.0" >> tmp.yaml
		echo "MaxIterations: 2000" >> tmp.yaml
	elif [ $order -eq 1 ]; then
		echo "CFL: 2.0" >> tmp.yaml
		echo "MaxIterations: 15000" >> tmp.yaml
	elif [ $order -eq 2 ]; then
		echo "MaxIterations: 20000" >> tmp.yaml
	elif [ $order -eq 4 ]; then
		echo "CFL: 2.5" >> tmp.yaml
		echo "MaxIterations: 60000" >> tmp.yaml
	fi
	#gocfd 2D -I tmp.yaml -F mesh_NACA0012_inv.su2 >& sysout-$order
	gocfd 2D -I tmp.yaml -F mesh_NACA0012_inv.su2
    if [ -f meshfile.gobcfd ]; then
    	mv meshfile.gobcfd meshfile-$order.gobcfd
    	mv solutionfile.gobcfd solutionfile-$order.gobcfd
    	mv plotfile.dat plotfile-$order.dat
	fi
done
if [ ! `which gnuplot` ]; then
	echo "You can install gnuplot and then use ./plot-cp.sh to plot a Mach and CP distribution"
	echo "To install gnuplot, use 'sudo apt install gnuplot'"
else
	echo "You can plot a CP and Mach chart using: ./plot-cp.sh plotfile-0.dat"
fi
