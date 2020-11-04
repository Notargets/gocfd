package Euler2D

import (
	"fmt"
	"math"
	"strconv"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

func TestEuler(t *testing.T) {
	c := NewEuler(1, 1, 1, "../../DG2D/test_tris_5.neu", FLUX_Average, FREESTREAM)
	K := c.dfr.K
	Np := c.dfr.SolutionElement.Np
	for n := 0; n < 4; n++ {
		for i := 0; i < Np; i++ {
			for k := 0; k < K; k++ {
				c.Q[n].Data()[k+i*K] = float64(k + 1)
			}
		}
	}
	// Interpolate from solution points to edges using precomputed interpolation matrix
	for n := 0; n < 4; n++ {
		c.Q_Face[n] = c.dfr.FluxInterpMatrix.Mul(c.Q[n])
	}
	Nedge := c.dfr.FluxElement.Nedge
	for n := 0; n < 4; n++ {
		for i := 0; i < 3*Nedge; i++ {
			for k := 0; k < K; k++ {
				assert.True(t, near(float64(k+1), c.Q_Face[n].Data()[k+i*K], 0.000001))
			}
		}
	}
}

func PrintQ(Q [4]utils.Matrix) {
	var (
		label string
	)
	for ii := 0; ii < 4; ii++ {
		switch ii {
		case 0:
			label = "Rho"
		case 1:
			label = "RhoU"
		case 2:
			label = "RhoV"
		case 3:
			label = "RhoE"
		}
		fmt.Println(Q[ii].Print(label))
	}
}
func PrintFlux(F []utils.Matrix) {
	for ii := 0; ii < len(F); ii++ {
		label := strconv.Itoa(ii)
		fmt.Println(F[ii].Print("F" + "[" + label + "]"))
	}
}

func nearVec(a, b []float64, tol float64) (l bool) {
	for i, val := range a {
		if !near(b[i], val, tol) {
			fmt.Printf("Diff = %v, Left[%d] = %v, Right[%d] = %v\n", math.Abs(val-b[i]), i, val, i, b[i])
			return false
		}
	}
	return true
}

func near(a, b float64, tolI ...float64) (l bool) {
	var (
		tol float64
	)
	if len(tolI) == 0 {
		tol = 1.e-08
	} else {
		tol = tolI[0]
	}
	bound := math.Max(tol, tol*math.Abs(a))
	if math.Abs(a-b) <= bound {
		l = true
	}
	return
}
