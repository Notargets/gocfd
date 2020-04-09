package main

import (
	"math"

	"github.com/notargets/gocfd/model_problems"

	"github.com/notargets/gocfd/DG1D"
)

const (
	K      = 10
	N      = 8
	NFaces = 2
	Nfp    = 1
)

func main() {
	VX, EToV := DG1D.SimpleMesh1D(0, 2, K)
	e1D := DG1D.NewElements1D(K, N, NFaces, Nfp, VX, EToV)
	c := model_problems.NewConvection(2*math.Pi, 0.75, 10., e1D)
	c.Run()
}
