#!/bin/bash
MeshSize=100pts
#for order in 0 1 2 3 4 5 6 7
for order in 6 7
do
  echo "Running order: $order"
  cat input-base.yaml > tmp.yaml
  if [ $order -eq 0 ]; then
	echo "CFL: 2.5" >> tmp.yaml
  elif [ $order -eq 1 ]; then
	echo "CFL: 1.0" >> tmp.yaml
  elif [ $order -eq 5 ]; then
	echo "CFL: 0.25" >> tmp.yaml
	echo "Kappa: 4.5" >> tmp.yaml
  elif [ $order -gt 5 ]; then
	echo "CFL: 0.15" >> tmp.yaml
	echo "Kappa: 4.5" >> tmp.yaml
  else
	echo "CFL: 5" >> tmp.yaml
  fi
  echo "PolynomialOrder: "$order >> tmp.yaml
  gocfd 2D -I tmp.yaml -F sod-aligned-$MeshSize.su2
  if [ -f shocktube.dat ]; then
  	./plot.sh
  	mv shocktube.dat shocktube-$MeshSize-order$order.dat
	echo "Enter to continue..."
	read arg
  else
    echo "Run failed"
  fi
done

#gocfd 2D -I input-wall.yaml -F sod-aligned-100pts.su2
#gocfd 2D -I input-wall.yaml -F sod-aligned-500pts.su2
