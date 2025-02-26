package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

func TestNewRTBasisSimplex(t *testing.T) {
	var (
		N = 1
	)
	tol := 0.0000001
	rt1 := NewRTElement(N, SimplexRTBasis)
	rt2 := NewRTElement(N, ErvinBasisRT)
	Np := rt1.Np
	Diff := utils.NewMatrix(Np, Np)
	for i := 0; i < rt1.Np; i++ {
		for j := 0; j < rt1.Np; j++ {
			d1 := rt1.V.At(i, j)
			d2 := rt2.V.At(i, j)
			if math.Abs(d1) < tol && math.Abs(d2) < tol {
				continue
			} else if (math.Abs(d1) > tol && math.Abs(d2) < tol) ||
				(math.Abs(d1) < tol && math.Abs(d2) > tol) {
				fmt.Printf("one is zero: %5.2f, %5.2f\n", d1, d2)
			}
			// fmt.Printf("d1,d2[%d,%d]: %5.2f, %5.2f - Ratio: %5.2f\n",
			// 	i, j, d1, d2, d1/d2)
			Diff.Set(i, j, d1/d2)
		}
	}
	Diff.Print("Diff")
	for i := 0; i < Np; i++ {
		for j := 0; j < Np; j++ {
			if math.Abs(Diff.At(i, j)) > tol {
				assert.InDeltaf(t, 1., Diff.At(i, j), tol, "")
			}
		}
	}
}
