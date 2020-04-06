package main

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
)

const (
	K      = 10
	N      = 8
	NFaces = 2
	Nfp    = 1
)

func main() {
	X := DG1D.Startup1D(K, N, NFaces, Nfp)
	U := utils.MatApply(X, math.Sin)
	fmt.Printf("U = \n%v\n", mat.Formatted(U, mat.Squeeze()))
}
