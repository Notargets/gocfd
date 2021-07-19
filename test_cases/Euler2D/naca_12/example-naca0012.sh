#!/bin/bash

printf "NACA 0012 Airfoil, Cp (Pressure Coefficient)\n"
((NUM_ITERATIONS_PER_FRAME=300))
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -g -z 0.08 -l -0.7 -k 0.3 -q 7
gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -g -l -0.7 -k 0.3 -q 7
