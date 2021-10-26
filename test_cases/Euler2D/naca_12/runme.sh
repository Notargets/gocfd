#!/bin/bash

gocfd 2D -I input-wall.yaml -F mesh/nacaAirfoil-base.su2 -g -s20 -q 4 -z 0.06 -x .25 -y .25 -k 0.3 -l 1.5
