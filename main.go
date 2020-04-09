package main

import (
	"math"

	"github.com/notargets/gocfd/model_problems"

	"github.com/notargets/gocfd/DG1D"
)

const (
	K      = 10 // Number of elements
	N      = 8  // Polynomial degree
	NFaces = 2  // Number of faces per element
	Nfp    = 1  // Number of points per face
)

func main() {
	VX, EToV := DG1D.SimpleMesh1D(0, 2, K)
	e1D := DG1D.NewElements1D(N, VX, EToV)
	c := model_problems.NewConvection(2*math.Pi, 0.75, 10., e1D)
	c.Run()
}
