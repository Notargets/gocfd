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
		Nmax := 7
		for N := 1; N <= Nmax; N++ {
			c := NewEuler(1, 1, N, "../../DG2D/test_tris_5.neu", FLUX_Average, FREESTREAM, false)
			K := c.dfr.K
			Nint := c.dfr.FluxElement.Nint
			Nedge := c.dfr.FluxElement.Nedge
			for n := 0; n < 4; n++ {
				for i := 0; i < Nint; i++ {
					for k := 0; k < K; k++ {
						c.Q[n].Data()[k+i*K] = float64(k + 1)
					}
				}
			}
			// Interpolate from solution points to edges using precomputed interpolation matrix
			for n := 0; n < 4; n++ {
				c.Q_Face[n] = c.dfr.FluxInterpMatrix.Mul(c.Q[n])
			}
			for n := 0; n < 4; n++ {
				for i := 0; i < 3*Nedge; i++ {
					for k := 0; k < K; k++ {
						assert.True(t, near(float64(k+1), c.Q_Face[n].Data()[k+i*K], 0.000001))
					}
				}
			}
		}
	}
	{ // Test solution process
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
			Nint := c.dfr.FluxElement.Nint
			Nedge := c.dfr.FluxElement.Nedge
			_ = Nedge
			// Calculate flux and project into R and S (transformed) directions
			for k := 0; k < c.dfr.K; k++ {
				for i := 0; i < Nint; i++ {
					Fr, Fs := c.CalculateFluxTransformed(k, i, c.Q)
					for n := 0; n < 4; n++ {
						rtD := c.F_RT_DOF[n].Data()
						ind := k + i*K
						rtD[ind], rtD[ind+Nint*K] = Fr[n], Fs[n]
					}
				}
			}
			// Interpolate from solution points to edges using precomputed interpolation matrix
			for n := 0; n < 4; n++ {
				c.Q_Face[n] = c.dfr.FluxInterpMatrix.Mul(c.Q[n])
			}
			/*
				for k := 0; k < c.dfr.K; k++ {
					for i := 0; i < 3*Nedge; i++ {
						FrE, FsE := c.CalculateFluxTransformed(k, i, c.Q_Face)
						for n := 0; n < 4; n++ {
							rtD := c.F_RT_DOF[n].Data()
							ind := k + (i+2*Nint)*K
							rtD[ind] = FrE[n] + FsE[n] // TODO: replace garbage here with actual edge calculations
						}
					}
				}
			*/
			PrintQ(c.Q, "Q_Int")
			PrintQ(c.Q_Face, "Q_Face")
			PrintQ(c.F_RT_DOF, "F_RT")
		}
	}
}

func PrintQ(Q [4]utils.Matrix, l string) {
	var (
		label string
	)
	for ii := 0; ii < 4; ii++ {
		switch ii {
		case 0:
			label = l + "_0"
		case 1:
			label = l + "_1"
		case 2:
			label = l + "_2"
		case 3:
			label = l + "_3"
		}
		fmt.Println(Q[ii].Transpose().Print(label))
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
