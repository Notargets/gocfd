#!/bin/bash

printf "NACA 0012 Airfoil, Cp (Pressure Coefficient)\n"
((NUM_ITERATIONS_PER_FRAME=50))
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -g -z 0.08 -l -0.7 -k 0.3 -q 7
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -z 0.03 -g -l -0.9 -k 0.5 -q 7
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -z 0.06 -g -q 7
gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -z 0.06 -k 0.1 -l 0.4 -g -q 4
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -z 0.03 -g -l 0 -k 0.7 -q 4
