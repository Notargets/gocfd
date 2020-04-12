package main

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"

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
	VX, EToV := DG1D.SimpleMesh1D(0, 20, K)
	e1D := DG1D.NewElements1D(N, VX, EToV)
	c := model_problems.NewConvection(2*math.Pi, 0.75, 10., e1D)
	c.Run()
	fmt.Printf("X = \n%v\n", mat.Formatted(e1D.X, mat.Squeeze()))
}
