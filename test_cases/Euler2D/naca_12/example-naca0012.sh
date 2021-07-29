#!/bin/bash

printf "NACA 0012 Airfoil, Cp (Pressure Coefficient)\n"
((NUM_ITERATIONS_PER_FRAME=5))
gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -g -z 0.08 -k -0.7 -l 0.3 -q 7
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -z 0.03 -g -k -0.9 -l 0.5 -q 7
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -z 0.06 -g -q 7
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -z 0.06 -k 0.1 -l 0.4 -g -q 4
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -g -q 7 -k -20 -l 20 -z 0.06
#gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -s "$NUM_ITERATIONS_PER_FRAME" -z 0.03 -g -k 0 -l 0.7 -q 4
