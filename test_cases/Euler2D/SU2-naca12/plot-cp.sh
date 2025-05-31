#!/bin/bash

filename=$1
if [ ! -f "$filename" ]; then
  echo "Must provide a plotfile.dat filename"
  exit 1
fi

#(awk -F, '$2 < 0' $filename | sort -t, -k1,1g; awk -F, '$2 >= 0' $filename | sort -t, -k1,1g)
cat $filename | (awk -F, '$2 < 0' $filename | sort -t, -k1,1g; awk -F, '$2 >= 0' $filename | sort -t, -k1,1gr) > plot.dat

gnuplot plot2.gp
