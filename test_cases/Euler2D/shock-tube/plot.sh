#!/bin/bash

rm /tmp/tmp1 /tmp/tmp2

if [ -z "$1" ]; then
  file="shocktube.dat"
else
  file="$1"
fi

ln -s `pwd`/$file /tmp/tmp1
ln -s `pwd`/shocktube_analytic.dat /tmp/tmp2

gnuplot -persist plot.gp 
