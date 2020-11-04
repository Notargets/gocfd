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
	{ // Test interpolation of solution to edges for all supported orders
		/*
			Solver approach:
			0) Solution is stored on sol points as Q
			0a) Flux is computed and stored in X, Y component projections in the 2*Nint front of F_RT_DOF
			1) Solution is extrapolated to edge points in Q_Face from Q
			2) Edges are traversed, flux is calculated and projected onto edge face normals, scaled and placed into F_RT_DOF
		*/
		Nmax := 1
		for N := 1; N <= Nmax; N++ {
			c := NewEuler(1, 1, N, "../../DG2D/test_tris_5.neu", FLUX_Average, FREESTREAM, false)
			K := c.dfr.K
			Nint := c.dfr.SolutionElement.Np
			for n := 0; n < 4; n++ {
				for i := 0; i < Nint; i++ {
					for k := 0; k < K; k++ {
						c.Q[n].Data()[k+i*K] = float64(k + 1)
					}
				}
			}
			// Calculate flux and project into R and S (transformed) directions
			var ind int
			for k := 0; k < c.dfr.K; k++ {
				for i := 0; i < Nint; i++ {
					Fr, Fs := c.CalculateFluxTransformed(k, i)
					for n := 0; n < 4; n++ {
						rtD := c.F_RT_DOF[n].Data()
						rtD[ind], rtD[ind+Nint*K] = Fr[n], Fs[n]
					}
					ind++
				}
			}
			PrintQ(c.F_RT_DOF)
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
