package DG1D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

func SimpleMesh1D(xmin, xmax float64, K int) (VX utils.Vector, EToV utils.Matrix) {
	// K is the number of elements, there are K+1 vertices
	var (
		x             = make([]float64, K+1)
		elementVertex = make([]float64, K*2)
	)
	for i := 0; i < K+1; i++ {
		x[i] = (xmax-xmin)*float64(i)/float64(K) + xmin
	}
	/*
		Example: K=4, 5 vertices
			EToV =
			⎡0.0000  1.0000⎤
			⎢1.0000  2.0000⎥
			⎢2.0000  3.0000⎥
			⎣3.0000  4.0000⎦
	*/
	var iter int
	for i := 0; i < K; i++ {
		elementVertex[iter] = float64(i)
		elementVertex[iter+1] = float64(i + 1)
		iter += 2
	}
	return utils.NewVector(K+1, x), utils.NewMatrix(K, 2, elementVertex)
}

func Gamma0(alpha, beta float64) float64 {
	ab1 := alpha + beta + 1.
	a1 := alpha + 1.
	b1 := beta + 1.
	return math.Gamma(a1) * math.Gamma(b1) * math.Pow(2, ab1) / ab1 / math.Gamma(ab1)
}
func Gamma1(alpha, beta float64) float64 {
	ab := alpha + beta
	a1 := alpha + 1.
	b1 := beta + 1.
	return a1 * b1 * Gamma0(alpha, beta) / (ab + 3.0)
}
