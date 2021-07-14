#!/bin/bash

printf "NACA 0012 Airfoil, Cp (Pressure Coefficient)\n"
((NUM_ITERATIONS_PER_FRAME=20))
gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -g -z 0.04 -l -1 -k 1.5 -q 7
