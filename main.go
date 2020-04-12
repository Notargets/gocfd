package main

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"

	"github.com/notargets/gocfd/model_problems"

	"github.com/notargets/gocfd/DG1D"
)

const (
	K = 10 // Number of elements
	N = 4  // Polynomial degree
)

func main() {
	VX, EToV := DG1D.SimpleMesh1D(0, 2*math.Pi, K)
	e1D := DG1D.NewElements1D(N, VX, EToV)
	c := model_problems.NewConvection(2*math.Pi, 0.75, 100000., e1D)
	c.Run()
	fmt.Printf("X = \n%v\n", mat.Formatted(e1D.X, mat.Squeeze()))
}
