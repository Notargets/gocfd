package Euler2D

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"sync"
	"testing"

	"github.com/notargets/gocfd/InputParameters"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/types"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

var ipDefault = &InputParameters.InputParameters2D{
	CFL:      2.5,
	FluxType: "Roe",
	InitType: "Freestream",
	Minf:     2,
	Gamma:    1.4,
	Limiter:  "PerssonC0",
	Kappa:    3,
}

func TestFluidFunctions(t *testing.T) {
	N := 1
	ip := *ipDefault
	ip.Minf = 2.
	ip.PolynomialOrder = N
	c := NewEuler(&ip, "../../DG2D/test_data/test_tris_6.neu", 1, false, false)
	funcs := []FlowFunction{Density, XMomentum, YMomentum, Energy, Mach, StaticPressure}
	values := make([]float64, len(funcs))
	for i, plotField := range []FlowFunction{Density, XMomentum, YMomentum, Energy, Mach, StaticPressure} {
		values[i] = c.FSFar.GetFlowFunction(c.Q[0], 0, plotField)
		// fmt.Printf("%s[%d] = %8.5f\n", plotField.String(), ik, values[i])
	}
	assert.InDeltaSlicef(t, []float64{1, 2, 0, 3.78571, 2, 0.71429}, values, 0.00001, "err msg %s")
}

func TestEuler(t *testing.T) {
	var (
		msg = "err msg %s"
		tol = 0.000001
	)
	ip := *ipDefault
	ip.FluxType = "average"
	if true {
		{ // Test interpolation of solution to edges for all supported orders
			Nmax := 7
			for N := 0; N <= Nmax; N++ {
				ip.PolynomialOrder = N
				// c := NewEuler(1, N, "../../DG2D/test_tris_5.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, None, false, false, false)
				c := NewEuler(&ip, "../../DG2D/test_data/test_tris_5.neu", 1, false, false)
				Kmax := c.DFR.K
				Nint := c.DFR.FluxElement.NpInt
				Nedge := c.DFR.FluxElement.NpEdge
				var Q, Q_Face [4]utils.Matrix
				for n := 0; n < 4; n++ {
					Q[n] = utils.NewMatrix(Nint, Kmax)
					Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
				}
				for n := 0; n < 4; n++ {
					for i := 0; i < Nint; i++ {
						for k := 0; k < Kmax; k++ {
							ind := k + i*Kmax
							Q[n].DataP[ind] = float64(k + 1)
						}
					}
				}
				// Interpolate from solution points to edges using precomputed interpolation matrix
				for n := 0; n < 4; n++ {
					Q_Face[n] = c.DFR.FluxEdgeInterp.Mul(Q[n])
				}
				for n := 0; n < 4; n++ {
					for i := 0; i < 3*Nedge; i++ {
						for k := 0; k < Kmax; k++ {
							ind := k + i*Kmax
							assert.InDeltaf(t, float64(k+1), Q_Face[n].DataP[ind], tol, msg)
						}
					}
				}
			}
		}
	}
	if true {
		{ // Test solution process
			/*
				Solver approach:
				0) Solution is stored on sol points as Q
				0a) Flux is computed and stored in R, S component projections in the 2*NpInt front of F_RT_DOF
				1) Solution is extrapolated to edge points in Q_Face from Q
				2) Edges are traversed, flux is calculated and projected onto edge face normals, scaled and placed into F_RT_DOF
			*/
			Nmax := 7
			for N := 0; N <= Nmax; N++ {
				// c := NewEuler(1, N, "../../DG2D/test_tris_5.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, None, false, false, false)
				ip.PolynomialOrder = N
				c := NewEuler(&ip, "../../DG2D/test_data/test_tris_5.neu", 1, false, false)
				Kmax := c.DFR.K
				Nint := c.DFR.FluxElement.NpInt
				Nedge := c.DFR.FluxElement.NpEdge
				NpFlux := c.DFR.FluxElement.Np // Np = 2*NpInt+3*NpEdge
				var Q_Face, F_RT_DOF [4]utils.Matrix
				for n := 0; n < 4; n++ {
					Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
					F_RT_DOF[n] = utils.NewMatrix(NpFlux, Kmax)
				}
				Q := c.Q[0]
				// Mark the initial state with the element number
				qD := [4][]float64{Q[0].DataP, Q[1].DataP, Q[2].DataP, Q[3].DataP}
				for i := 0; i < Nint; i++ {
					for k := 0; k < Kmax; k++ {
						ind := k + i*Kmax
						qD[0][ind] = float64(k + 1)
						qD[1][ind] = 0.1 * float64(k+1)
						qD[2][ind] = 0.05 * float64(k+1)
						qD[3][ind] = 2.00 * float64(k+1)
					}
				}
				// Flux values for later checks are invariant with i (i=0)
				Fr_check, Fs_check := make([][4]float64, Kmax), make([][4]float64, Kmax)
				for k := 0; k < Kmax; k++ {
					Fr_check[k], Fs_check[k] = c.CalculateFluxTransformed(k, Kmax, 0, c.DFR.Jdet, c.DFR.Jinv, Q)
				}
				// Interpolate from solution points to edges using precomputed interpolation matrix
				for n := 0; n < 4; n++ {
					Q_Face[n] = c.DFR.FluxEdgeInterp.Mul(Q[n])
				}
				// Calculate flux and project into R and S (transformed) directions
				rtD := [4][]float64{F_RT_DOF[0].DataP, F_RT_DOF[1].DataP, F_RT_DOF[2].DataP, F_RT_DOF[3].DataP}
				for n := 0; n < 4; n++ {
					for i := 0; i < Nint; i++ {
						for k := 0; k < Kmax; k++ {
							ind := k + i*Kmax
							Fr, Fs := c.CalculateFluxTransformed(k, Kmax, i, c.DFR.Jdet, c.DFR.Jinv, Q)
							rtD[n][ind], rtD[n][ind+Nint*Kmax] = Fr[n], Fs[n]
						}
					}
					// Check to see that the expected values are in the right place (the internal locations)
					rtTD := F_RT_DOF[n].Transpose().DataP
					for k := 0; k < Kmax; k++ {
						val0, val1 := Fr_check[k][n], Fs_check[k][n]
						is := k * NpFlux
						assert.True(t, nearVecScalar(rtTD[is:is+Nint], val0, 0.000001))
						is += Nint
						assert.True(t, nearVecScalar(rtTD[is:is+Nint], val1, 0.000001))
					}
					// Set normal flux to a simple addition of the two sides to use as a check in assert()
					for k := 0; k < Kmax; k++ {
						for i := 0; i < 3*Nedge; i++ {
							ind := k + (2*Nint+i)*Kmax
							Fr, Fs := c.CalculateFluxTransformed(k, Kmax, i, c.DFR.Jdet, c.DFR.Jinv, Q_Face)
							rtD[n][ind] = Fr[n] + Fs[n]
						}
					}
					// Check to see that the expected values are in the right place (the edge locations)
					rtTD = F_RT_DOF[n].Transpose().DataP
					for k := 0; k < Kmax; k++ {
						val := Fr_check[k][n] + Fs_check[k][n]
						is := k * NpFlux
						ie := (k + 1) * NpFlux
						assert.True(t, nearVecScalar(rtTD[is+2*Nint:ie], val, 0.000001))
					}
				}
			}
		}
	}
	if true {
		{ // Test solution process part 2 - Freestream divergence should be zero
			Nmax := 7
			for N := 0; N <= Nmax; N++ {
				// c := NewEuler(1, N, "../../DG2D/test_tris_5.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, None, false, false, false)
				ip.PolynomialOrder = N
				// ip.Minf = 0.301
				// ip.Alpha = 2
				// c := NewEuler(&ip, "../../DG2D/test_tris_5.neu", 1, plotMesh, false, false)
				c := NewEuler(&ip, "../../DG2D/test_data/test_tris_6_nowall.neu", 1, false, false)
				c.FSIn = c.FSFar
				Kmax := c.DFR.K
				Nint := c.DFR.FluxElement.NpInt
				Nedge := c.DFR.FluxElement.NpEdge
				NpFlux := c.DFR.FluxElement.Np // Np = 2*NpInt+3*NpEdge
				// Mark the initial state with the element number
				var Q_Face, Q_Face_P0, F_RT_DOF [4]utils.Matrix
				var QMean [4]utils.Vector
				var Flux, Flux_Face [2][4]utils.Matrix
				for n := 0; n < 4; n++ {
					F_RT_DOF[n] = utils.NewMatrix(NpFlux, Kmax)
					Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
					Q_Face_P0[n] = utils.NewMatrix(3*Nedge, Kmax)
					Flux_Face[0][n] = utils.NewMatrix(3*Nedge, Kmax)
					Flux_Face[1][n] = utils.NewMatrix(3*Nedge, Kmax)
					Flux[0][n] = utils.NewMatrix(Nint, Kmax)
					Flux[1][n] = utils.NewMatrix(Nint, Kmax)
					QMean[n] = utils.NewVector(Kmax)
				}
				Q := c.Q[0]
				c.SetRTFluxInternal(Kmax, c.DFR.Jdet, c.DFR.Jinv, F_RT_DOF, Q)
				c.UpdateElementMean(Q, QMean)
				c.InterpolateSolutionToEdges(Q, Q_Face, Q_Face_P0)
				EdgeQ1 := make([][4]float64, Nedge)
				EdgeQ2 := make([][4]float64, Nedge)
				c.CalculateEdgeEulerFlux(0,
					[][4]utils.Matrix{Q_Face},
					[][4]utils.Vector{QMean},
					[][2][4]utils.Matrix{Flux_Face},
					EdgeQ1, EdgeQ2, c.SortedEdgeKeys[0])
				c.SetRTFluxOnEdges(0, Kmax, F_RT_DOF)
				// Check that freestream divergence on this mesh is zero
				for n := 0; n < 4; n++ {
					var div utils.Matrix
					div = c.DFR.FluxElement.DivInt.Mul(F_RT_DOF[n])
					for k := 0; k < Kmax; k++ {
						for i := 0; i < Nint; i++ {
							ind := k + i*Kmax
							div.DataP[ind] /= c.DFR.Jdet.DataP[k]
						}
					}
					fmt.Printf("checking variable[%d]...", n)
					assert.True(t, nearVecScalar(div.DataP, 0., 0.000001))
					fmt.Printf("done.\n")
				}
			}
		}
		if true { // Test divergence of polynomial initial condition against
			// analytic values
			/*
				Note: the Polynomial flux is asymmetric around the R and S axes - it uses abs(x) and abs(y)
				Elements should not straddle the axes if a perfect polynomial flux capture is needed
			*/
			Nmax := 7
			for N := 2; N <= Nmax; N++ {
				// Single triangle test case
				var c *Euler
				// c = NewEuler(1, N, "../../DG2D/test_tris_1tri.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, None, plotMesh, false, false)
				ip.PolynomialOrder = N
				c = NewEuler(&ip, "../../DG2D/test_data/test_tris_1tri.neu", 1, false, false)
				CheckFlux0(c, t)
				// Two widely separated triangles - no shared faces
				// c = NewEuler(1, N, "../../DG2D/test_tris_two.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, None, plotMesh, false, false)
				c = NewEuler(&ip, "../../DG2D/test_data/test_tris_two.neu", 1, false, false)
				CheckFlux0(c, t)
				// Two widely separated triangles - no shared faces - one tri listed in reverse order
				// c = NewEuler(1, N, "../../DG2D/test_tris_twoR.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, None, plotMesh, false, false)
				c = NewEuler(&ip, "../../DG2D/test_data/test_tris_twoR.neu", 1, false, false)
				CheckFlux0(c, t)
				// Connected tris, sharing one edge
				// plotMesh = true
				// c = NewEuler(1, N, "../../DG2D/test_tris_6_nowall.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, None, plotMesh, false, false)
				c = NewEuler(&ip, "../../DG2D/test_data/test_tris_6_nowall.neu", 1, false, false)
				CheckFlux0(c, t)
			}
		}
	}
	if true {
		{ // Test divergence of Isentropic Vortex initial condition against analytic values - density equation only
			N := 1
			{
				// Nmax := 7
				// for N := 0; N <= Nmax; N++ {
				ip.PolynomialOrder = N
				ip.InitType = "ivortex"
				// c := NewEuler(1, N, "../../DG2D/test_tris_6.neu", 1, FLUX_Average, IVORTEX, 1, 0, 1.4, 0, false, 5000, None, plotMesh, false, false)
				c := NewEuler(&ip, "../../DG2D/test_data/test_tris_6.neu", 1, false, false)
				for _, e := range c.DFR.Tris.Edges {
					if e.BCType == types.BC_IVortex {
						e.BCType = types.BC_None
					}
				}
				Kmax := c.DFR.K
				Nint := c.DFR.FluxElement.NpInt
				Nedge := c.DFR.FluxElement.NpEdge
				NpFlux := c.DFR.FluxElement.Np // Np = 2*NpInt+3*NpEdge
				// Mark the initial state with the element number
				var Q_Face, Q_Face_P0, F_RT_DOF [4]utils.Matrix
				var QMean [4]utils.Vector
				var Flux, Flux_Face [2][4]utils.Matrix
				for n := 0; n < 4; n++ {
					F_RT_DOF[n] = utils.NewMatrix(NpFlux, Kmax)
					Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
					Q_Face_P0[n] = utils.NewMatrix(3*Nedge, Kmax)
					Flux_Face[0][n] = utils.NewMatrix(3*Nedge, Kmax)
					Flux_Face[1][n] = utils.NewMatrix(3*Nedge, Kmax)
					Flux[0][n] = utils.NewMatrix(Nint, Kmax)
					Flux[1][n] = utils.NewMatrix(Nint, Kmax)
					QMean[n] = utils.NewVector(Kmax)
				}
				Q := c.Q[0]
				X, Y := c.DFR.FluxX, c.DFR.FluxY
				c.SetRTFluxInternal(Kmax, c.DFR.Jdet, c.DFR.Jinv, F_RT_DOF, Q)
				c.InterpolateSolutionToEdges(Q, Q_Face, Q_Face_P0)
				c.UpdateElementMean(Q, QMean)
				EdgeQ1 := make([][4]float64, Nedge)
				EdgeQ2 := make([][4]float64, Nedge)
				c.CalculateEdgeEulerFlux(0, [][4]utils.Matrix{Q_Face},
					[][4]utils.Vector{QMean},
					[][2][4]utils.Matrix{Flux_Face},
					EdgeQ1, EdgeQ2, c.SortedEdgeKeys[0])
				c.SetRTFluxOnEdges(0, Kmax, F_RT_DOF)
				var div utils.Matrix
				// Density is the easiest equation to match with a polynomial
				n := 0
				// fmt.Printf("component[%d]\n", n)
				div = c.DFR.FluxElement.DivInt.Mul(F_RT_DOF[n])
				// c.DivideByJacobian(Kmax, NpInt, c.DFR.Jdet, div.DataP, 1)
				for k := 0; k < Kmax; k++ {
					for i := 0; i < Nint; i++ {
						ind := k + i*Kmax
						div.DataP[ind] /= -c.DFR.Jdet.DataP[k]
					}
				}
				// Get the analytic values of divergence for comparison
				for k := 0; k < Kmax; k++ {
					for i := 0; i < Nint; i++ {
						ind := k + i*Kmax
						x, y := X.DataP[ind], Y.DataP[ind]
						qc1, qc2, qc3, qc4 := c.AnalyticSolution.GetStateC(0, x, y)
						q1, q2, q3, q4 := Q[0].DataP[ind], Q[1].DataP[ind], Q[2].DataP[ind], Q[3].DataP[ind]
						assert.InDeltaSlicef(t, []float64{q1, q2, q3, q4}, []float64{qc1, qc2, qc3, qc4}, tol, msg)
						divC := c.AnalyticSolution.GetDivergence(0, x, y)
						divCalc := div.DataP[ind]
						assert.InDeltaf(t, divCalc/qc1, divC[n]/qc1, 0.001, msg) // 0.1 percent match
					}
				}
			}
		}
	}
}

func TestFluxInterpolation(t *testing.T) {
	N := 2
	ip := *ipDefault
	ip.Minf = 1.
	ip.PolynomialOrder = N
	ip.FluxType = "Roe"
	c := NewEuler(&ip, "../../DG2D/test_data/test_tris_6.neu", 1, false, false)
	rk := c.NewRungeKuttaSSP()
	c.InterpolateSolutionToEdges(c.Q[0], rk.Q_Face[0], rk.Q_Face_P0[0])
	el := c.DFR.SolutionElement
	/*
		el.JB2D.V.String("V")
		el.JB2D.Vinv.String("Vinv")
		el.MassMatrix.String("M")
	*/
	RR := [][2]float64{{-1 / 3, -1 / 3}, {-1, -1}, {-1, 1}, {1, -1}}
	NpInt := el.Np
	locations := utils.NewMatrix(len(RR), NpInt)
	for ii, rr := range RR {
		var sm int
		for i := 0; i <= N; i++ {
			for j := 0; j <= N-i; j++ {
				// fmt.Printf("i,j = %d,%d\n", i, j)
				r, s := rr[0], rr[1]
				// fmt.Printf("Basis value at [%3.1f,%3.1f][%d,%d] = %5.3f\n", r, s, i, j, el.JB2D.PolynomialTerm(r, s, i, j))
				locations.Set(ii, sm, el.JB2D.PolynomialTerm(r, s, i, j))
				sm++
			}
		}
	}
	filter := make([]float64, 10)
	filter[0] = 1
	filter[1] = 1.0
	filter[2] = 1.0
	filter[3] = 0.00
	filter[4] = 0.00
	filter[5] = 0.00
	filter[6] = 0.00
	filter[7] = 0.00
	filter[8] = 0.00
	filter[9] = 0.00
	truncate := utils.NewDiagMatrix(NpInt, filter[:NpInt])

	// solution := utils.NewMatrix(NpInt, 1, []float64{1, 1, 1, 1, 1, 1})
	// modes := el.JB2D.Vinv.Mul(solution)
	// solution.Transpose().String("solution-const")
	// modes.Transpose().String("modes1")
	// locations.Mul(modes).Transpose().String("Locations1")
	// solution = utils.NewMatrix(NpInt, 1, []float64{1, 1, 1, 0.5, 0.5, 1.5})
	solution := utils.NewMatrix(NpInt, 1, []float64{1.0, 0.5, 1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1})
	modes := el.JB2D.Vinv.Mul(solution)
	solution.Transpose().Print("solution-variable")
	modes.Transpose().Print("modes2")
	locations.Mul(modes).Transpose().Print("Locations2")

	VinvFiltered := truncate.Mul(el.JB2D.Vinv)
	modes = VinvFiltered.Mul(solution)
	// modes.Transpose().String("modes3")
	locations.Mul(modes).Transpose().Print("Locations3")

	R, S := utils.NewVector(4, []float64{-1 / 3, -1, -1, 1}), utils.NewVector(4, []float64{-1 / 3, -1, 1, -1})
	linterp := el.JB2D.GetInterpMatrix(R, S)
	linterp.Mul(solution).Transpose().Print("linterp-locations")
}

func TestFluxJacobian(t *testing.T) {
	var (
		tol = 0.000001
		msg = "err msg %s"
	)
	{ // Flux Jacobian calculation
		c := Euler{}
		c.FSFar = NewFreeStream(0.1, 1.4, 0)
		Qinf := c.FSFar.Qinf
		Fx, Gy := c.FluxJacobianCalc(Qinf[0], Qinf[1], Qinf[2], Qinf[3])
		// Matlab: using the FSFar Q = [1,.1,0,1.79071]
		assert.InDeltaSlicef(t, []float64{
			0, 1.0000, 0, 0,
			-0.0080, 0.1600, 0, 0.4000,
			0, 0, 0.1000, 0,
			-0.2503, 2.5010, 0, 0.1400,
		}, Fx[:], tol, msg)
		assert.InDeltaSlicef(t, []float64{
			0, 0, 1.0000, 0,
			0, 0, 0.1000, 0,
			0.0020, -0.0400, 0, 0.4000,
			0, 0, 2.5050, 0,
		}, Gy[:], tol, msg)
	}
}

func TestDissipation(t *testing.T) {
	{
		dfr := DG2D.NewDFR2D(1, false, "../../DG2D/test_data/test_tris_9.neu")
		VtoE := NewVertexToElement(dfr.Tris.EToV)
		vToE_test := VertexToElement{{0, 0, 0}, {0, 1, 0}, {1, 3, 0}, {1,
			1, 0}, {1, 2, 0}, {2, 4, 0}, {2, 3, 0}, {3, 5, 0}, {3, 0, 0}, {4,
			2, 0}, {4, 5, 0}, {4, 6, 0}, {4, 1, 0}, {4, 0, 0}, {4, 7, 0},
			{5, 4, 0}, {5, 3, 0}, {5, 9, 0}, {5, 8, 0}, {5, 2, 0}, {5, 7, 0},
			{6, 4, 0}, {6, 9, 0}, {7, 5, 0}, {7, 6, 0}, {8, 8, 0}, {8, 6, 0},
			{8, 7, 0}, {9, 8, 0}, {9, 9, 0}}
		sort.Sort(vToE_test)
		assert.Equal(t, vToE_test, VtoE)
		vepFinal := [3]int32{9, 9, 0}
		for NPar := 1; NPar < 10; NPar += 2 {
			pm := utils.NewPartitionMap(NPar, dfr.K)
			sd := NewScalarDissipation(0, dfr, pm)
			var vep [3]int32
			for np := 0; np < NPar; np++ {
				for _, val := range sd.VtoE[np] {
					vep = val
				}
			}
			assert.Equal(t, vepFinal[0], vep[0])
		}
		KMax := dfr.K
		assert.Equal(t, len(dfr.Jdet.DataP), KMax)
		for k := 0; k < KMax; k++ {
			area := 2. * dfr.Jdet.DataP[k] // Area of element is 2x Determinant, because the unit triangle is area=2
			assert.InDeltaf(t, area, 0.25, 0.000001, "err msg %s")
		}
	}
	if false { // Turn off value check tests while working on the constants in the artificial dissipation
		dfr := DG2D.NewDFR2D(2, false, "../../DG2D/test_data/test_tris_9.neu")
		Np, KMax := dfr.SolutionElement.Np, dfr.K
		pMap := utils.NewPartitionMap(1, KMax)
		Q := make([][4]utils.Matrix, 1)
		n := 0
		Q[0][n] = utils.NewMatrix(dfr.SolutionElement.Np, KMax)
		var val float64
		for k := 0; k < KMax; k++ {
			for i := 0; i < Np; i++ {
				if i < Np/2 {
					val = 2
				} else {
					val = 1
				}
				ind := k + KMax*i
				Q[0][n].DataP[ind] = val
			}
		}
		sd := NewScalarDissipation(0, dfr, pMap)
		sd.Kappa = 4.
		sd.CalculateElementViscosity(0, Q)
		// assert.InDeltaSlicef(t, []float64{0.09903, 0.09903, 0.09903, 0.09903, 0.09903, 0.09903, 0.09903, 0.09903, 0.09903, 0.09903},
		assert.InDeltaSlicef(t, []float64{0.09903, 0.09903, 0.09903, 0.07003, 0.09903, 0.07003, 0.09903, 0.07003, 0.09903, 0.07003},
			sd.EpsilonScalar[0], 0.00001, "err msg %s")
		sd.Kappa = 0.75
		sd.CalculateElementViscosity(0, Q)
		// assert.InDeltaSlicef(t, []float64{0.01270, 0.01270, 0.01270, 0.01270, 0.01270, 0.01270, 0.01270, 0.01270, 0.01270, 0.01270},
		assert.InDeltaSlicef(t, []float64{0.01270, 0.01270, 0.01270, 0.00898, 0.01270, 0.00898, 0.01270, 0.00898, 0.01270, 0.00898},
			sd.EpsilonScalar[0], 0.00001, "err msg %s")
	}
	{
		dfr := DG2D.NewDFR2D(1, false, "../../DG2D/test_data/test_tris_9.neu")
		_, KMax := dfr.SolutionElement.Np, dfr.K
		for NP := 1; NP < 5; NP++ {
			pm := utils.NewPartitionMap(NP, KMax)
			sd := NewScalarDissipation(0, dfr, pm)
			/*
				assert.InDeltaSlicef(t, []float64{0.666666, 0.166666, 0.166666, 0.166666, 0.666666, 0.166666, 0.166666, 0.166666,
					0.666666, 0.666666, 0.166666, 0.166666, 0.166666, 0.666666, 0.166666, 0.166666, 0.166666, 0.666666, 0.827326,
					0.172673, 0, 0.5, 0.5, 0, 0.172673, 0.827326, 0, 0, 0.827326, 0.172673, 0, 0.5, 0.5, 0, 0.172673, 0.827326,
					0.172673, 0, 0.827326, 0.5, 0, 0.5, 0.827326, 0, 0.172673},
					sd.BaryCentricCoords.DataP, 0.00001, "err msg %s")
			*/
			// Set the epsilon scalar value to the element ID and check the vertex aggregation
			for np := 0; np < NP; np++ {
				KMaxLocal := pm.GetBucketDimension(np)
				for k := 0; k < KMaxLocal; k++ {
					kGlobal := pm.GetGlobalK(k, np)
					sd.EpsilonScalar[np][k] = float64(kGlobal)
				}
			}
			wg := &sync.WaitGroup{}
			for np := 0; np < NP; np++ {
				wg.Add(1)
				sd.propagateEpsilonMaxToVertices(np)
			}
			assert.Equal(t, []float64{1, 3, 4, 5, 7, 9, 9, 6, 8, 9}, sd.EpsVertex)
		}
	}
}

func TestDissipation2(t *testing.T) {
	// Test C0 continuity of Epsilon using element vertex aggregation
	{
		var (
			dfr = DG2D.NewDFR2D(1, false, "../../DG2D/test_data/test_tris_9.neu")
		)
		NP := 1
		_, KMax := dfr.SolutionElement.Np, dfr.K
		pm := utils.NewPartitionMap(NP, KMax)
		sd := NewScalarDissipation(0, dfr, pm)
		for np := 0; np < NP; np++ {
			KMax = sd.PMap.GetBucketDimension(np)
			for k := 0; k < KMax; k++ {
				kGlobal := pm.GetGlobalK(k, np)
				sd.EpsilonScalar[np][k] = float64(kGlobal)
			}
		}
		wg := &sync.WaitGroup{}
		for np := 0; np < NP; np++ {
			wg.Add(1)
			sd.propagateEpsilonMaxToVertices(np)
		}
		assert.Equal(t, []float64{1, 3, 4, 5, 7, 9, 9, 6, 8, 9}, sd.EpsVertex)
		for np := 0; np < NP; np++ {
			sd.linearInterpolateEpsilon(np)
			// sd.baryCentricInterpolateEpsilon(np)
		}
		/*
			assert.InDeltaSlicef(t, []float64{
				5.66667, 3.33333, 7.66667, 7.16667, 8.16667, 6.00000, 7.50000, 8.00000, 8.83333, 9.00000, 4.66667, 5.33333,
				6.66667, 4.16667, 8.16667, 5.50000, 6.50000, 7.50000, 8.33333, 9.00000, 2.66667, 2.33333, 4.66667, 4.66667,
				5.66667, 6.50000, 7.00000, 8.50000, 8.83333, 9.00000, 5.66667, 3.33333, 7.66667, 7.16667, 8.16667, 6.00000,
				7.50000, 8.00000, 8.83333, 9.00000, 4.66667, 5.33333, 6.66667, 4.16667, 8.16667, 5.50000, 6.50000, 7.50000,
				8.33333, 9.00000, 2.66667, 2.33333, 4.66667, 4.66667, 5.66667, 6.50000, 7.00000, 8.50000, 8.83333, 9.00000,
				6.65465, 3.69069, 8.65465, 7.96396, 9.00000, 5.82733, 7.65465, 7.82733, 8.82733, 9.00000, 6.00000, 5.00000,
				8.00000, 6.00000, 9.00000, 5.50000, 7.00000, 7.50000, 8.50000, 9.00000, 5.34535, 6.30931, 7.34535, 4.03604,
				9.00000, 5.17267, 6.34535, 7.17267, 8.17267, 9.00000, 4.30931, 5.96396, 6.30931, 3.17267, 8.13663, 5.34535,
				6.17267, 7.34535, 8.17267, 9.00000, 3.00000, 4.00000, 5.00000, 3.50000, 6.50000, 6.00000, 6.50000, 8.00000,
				8.50000, 9.00000, 1.69069, 2.03604, 3.69069, 3.82733, 4.86337, 6.65465, 6.82733, 8.65465, 8.82733, 9.00000,
				2.03604, 1.34535, 4.03604, 4.86337, 4.86337, 6.82733, 7.17267, 8.82733, 9.00000, 9.00000, 4.00000, 2.00000,
				6.00000, 6.50000, 6.50000, 6.50000, 7.50000, 8.50000, 9.00000, 9.00000, 5.96396, 2.65465, 7.96396, 8.13663,
				8.13663, 6.17267, 7.82733, 8.17267, 9.00000, 9.00000},
				sd.Epsilon[0].DataP, 0.00001, "err msg %s")
		*/
		// fmt.Printf(sd.Epsilon[0].String("Epsilon"))
	}
	// Gradient test using GetSolutionGradient()
	{
		ip := *ipDefault
		ip.FluxType = "average"
		// Testing to fourth order in R and S
		ip.PolynomialOrder = 4
		c := NewEuler(&ip, "../../DG2D/test_data/test_tris_5.neu", 1, false, false)
		var (
			dfr                      = c.DFR
			Kmax                     = dfr.K
			fel                      = dfr.FluxElement
			NpInt                    = fel.NpInt
			Nedge                    = fel.NpEdge
			NpFlux                   = fel.Np
			Q, Q_Face                [4]utils.Matrix
			QGradXCheck, QGradYCheck [4]utils.Matrix
		)
		for n := 0; n < 4; n++ {
			Q[n] = utils.NewMatrix(NpInt, Kmax)
			Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
			QGradXCheck[n] = utils.NewMatrix(NpFlux, Kmax)
			QGradYCheck[n] = utils.NewMatrix(NpFlux, Kmax)
		}
		// Each variable [n] is an nth order polynomial in R and S
		for i := 0; i < NpFlux; i++ {
			for k := 0; k < Kmax; k++ {
				ind := k + i*Kmax
				X, Y := dfr.FluxX.DataP[ind], dfr.FluxY.DataP[ind]
				if i < NpInt {
					Q[0].DataP[ind] = X + Y
					Q[1].DataP[ind] = X*X + Y*Y
					Q[2].DataP[ind] = X*X*X + Y*Y*Y
					Q[3].DataP[ind] = X*X*X*X + Y*Y*Y*Y
				}
				// First order in R and S
				QGradXCheck[0].DataP[ind] = 1.
				QGradYCheck[0].DataP[ind] = 1.
				// Second order in R and S
				QGradXCheck[1].DataP[ind] = 2 * X
				QGradYCheck[1].DataP[ind] = 2 * Y
				// Third order in R and S
				QGradXCheck[2].DataP[ind] = 3 * X * X
				QGradYCheck[2].DataP[ind] = 3 * Y * Y
				// Fourth order in R and S
				QGradXCheck[3].DataP[ind] = 4 * X * X * X
				QGradYCheck[3].DataP[ind] = 4 * Y * Y * Y
			}
		}
		// Interpolate from solution points to edges using precomputed interpolation matrix
		fI := dfr.FluxEdgeInterp.Mul
		for n := 0; n < 4; n++ {
			Q_Face[n] = fI(Q[n])
		}
		// Test hard coded loop GradX and GradY
		{
			for n := 0; n < 4; n++ {
				var (
					DXmd, DYmd   = dfr.DXMetric.DataP, dfr.DYMetric.DataP
					DOFX, DOFY   = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
					DOFXd, DOFYd = DOFX.DataP, DOFY.DataP
				)
				var Un float64
				for k := 0; k < Kmax; k++ {
					for i := 0; i < NpFlux; i++ {
						ind := k + i*Kmax
						switch {
						case i < NpInt: // The first NpInt points are the solution element nodes
							Un = Q[n].DataP[ind]
						case i >= NpInt && i < 2*NpInt: // The second NpInt points are duplicates of the first NpInt values
							Un = Q[n].DataP[ind-NpInt*Kmax]
						case i >= 2*NpInt:
							Un = Q_Face[n].DataP[ind-2*NpInt*Kmax] // The last 3*NpEdge points are the edges in [0-1,1-2,2-0] order
						}
						DOFXd[ind] = DXmd[ind] * Un
						DOFYd[ind] = DYmd[ind] * Un
					}
				}
				DX := dfr.FluxElement.Div.Mul(DOFX) // R Derivative, divergence x RT_DOF is R derivative for this DOF
				DY := dfr.FluxElement.Div.Mul(DOFY) // S Derivative, divergence x RT_DOF is S derivative for this DOF
				fmt.Printf("Order[%d] check ...", n+1)
				assert.Equal(t, len(DX.DataP), len(QGradXCheck[n].DataP))
				assert.Equal(t, len(DY.DataP), len(QGradYCheck[n].DataP))
				assert.InDeltaSlicef(t, DX.DataP, QGradXCheck[n].DataP, 0.000001, "err msg %s")
				assert.InDeltaSlicef(t, DY.DataP, QGradYCheck[n].DataP, 0.000001, "err msg %s")
				fmt.Printf("... validates\n")
			}
		}
		// Test function for GradX and GradY
		{
			var (
				DX, DY     = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
				DOFX, DOFY = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
			)
			// Before this call, we need to load edge data into the edge store
			edgeValues := make([][4]float64, Nedge)
			myThread := -1
			for k := 0; k < Kmax; k++ {
				kGlobal := c.Partitions.GetGlobalK(k, myThread)
				for edgeNum := 0; edgeNum < 3; edgeNum++ {
					en := dfr.EdgeNumber[kGlobal+Kmax*edgeNum]
					primeElement := c.DFR.Tris.Edges[en].ConnectedTris[0]
					if int(primeElement) == kGlobal {
						for i := 0; i < Nedge; i++ {
							ind := k + (i+edgeNum*Nedge)*Kmax
							for n := 0; n < 4; n++ {
								edgeValues[i][n] = Q_Face[n].DataP[ind]
							}
						}
						c.EdgeStore.PutEdgeValues(en, EdgeQValues, edgeValues)
					}
				}
			}
			for n := 0; n < 4; n++ {
				fmt.Printf("Order[%d] check ...", n+1)
				c.GetSolutionGradient(-1, n, Q, DX, DY, DOFX, DOFY)
				assert.Equal(t, len(DX.DataP), len(QGradXCheck[n].DataP))
				assert.Equal(t, len(DY.DataP), len(QGradYCheck[n].DataP))
				assert.InDeltaSlicef(t, DX.DataP, QGradXCheck[n].DataP, 0.000001, "err msg %s")
				assert.InDeltaSlicef(t, DY.DataP, QGradYCheck[n].DataP, 0.000001, "err msg %s")
				fmt.Printf("... validates\n")
			}
		}
	}
}

func TestEuler_GetSolutionGradientUsingRTElement(t *testing.T) {
	{
		ip := *ipDefault
		ip.FluxType = "average"
		// Testing to fourth order in R and S
		ip.PolynomialOrder = 4
		ip.Minf = 0.8
		ip.Alpha = 2.
		c := NewEuler(&ip, "../../DG2D/test_data/test_tris_9.neu", 1, false, false)
		rk := c.NewRungeKuttaSSP()
		myThread := 0
		var (
			dfr                      = c.DFR
			Kmax                     = dfr.K
			fel                      = dfr.FluxElement
			Nedge                    = fel.NpEdge
			NpFlux                   = fel.Np
			QGradXCheck, QGradYCheck [4]utils.Matrix
			Q0                       = c.Q[myThread]
			SortedEdgeKeys           = c.SortedEdgeKeys[myThread]
			EdgeQ1                   = make([][4]float64, Nedge) // Local working memory
			EdgeQ2                   = make([][4]float64, Nedge) // Local working memory
		)
		for n := 0; n < 4; n++ {
			QGradXCheck[n] = utils.NewMatrix(NpFlux, Kmax)
			QGradYCheck[n] = utils.NewMatrix(NpFlux, Kmax)
		}
		// Flow is freestream, graident should be 0 in all variables
		for n := 0; n < 4; n++ {
			for i := 0; i < NpFlux; i++ {
				for k := 0; k < Kmax; k++ {
					ind := k + i*Kmax
					QGradXCheck[n].DataP[ind] = 0.
					QGradYCheck[n].DataP[ind] = 0.
				}
			}
		}
		// Interpolate from solution points to edges using precomputed interpolation matrix
		// Test function for GradX and GradY
		{
			var (
				DX, DY     = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
				DOFX, DOFY = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
			)
			c.CalculateEdgeEulerFlux(rk.Time, rk.Q_Face, rk.QMean,
				rk.Flux_Face,
				EdgeQ1, EdgeQ2, SortedEdgeKeys) // Global
			for n := 0; n < 4; n++ {
				fmt.Printf("Variable[%d] check ...", n+1)
				c.GetSolutionGradient(-1, n, Q0, DX, DY, DOFX, DOFY)
				assert.Equal(t, len(DX.DataP), len(QGradXCheck[n].DataP))
				assert.Equal(t, len(DY.DataP), len(QGradYCheck[n].DataP))
				assert.InDeltaSlicef(t, DX.DataP, QGradXCheck[n].DataP, 0.000001, "err msg %s")
				assert.InDeltaSlicef(t, DY.DataP, QGradYCheck[n].DataP, 0.000001, "err msg %s")
				fmt.Printf("... validates\n")
			}
		}
	}
}

func TestInputParameters_Parse(t *testing.T) {
	var (
		err error
	)
	fileInput := []byte(`
Title: Test Case
CFL: 1.
InitType: Freestream # Can be IVortex or Freestream
FluxType: Roe
PolynomialOrder: 2
FinalTime: 4.
BCs: 
  Inflow:
      37:
         NPR: 4.0
  Outflow:
      22:
         P: 1.5
`)
	var input InputParameters.InputParameters2D
	if err = input.Parse(fileInput); err != nil {
		panic(err)
	}
	// Check Inflow BC number 37
	assert.Equal(t, input.BCs["Inflow"][37]["NPR"], 4.)
	// Check Outflow BC number 22
	assert.Equal(t, input.BCs["Outflow"][22]["P"], 1.5)
	input.Print()
	assert.Equal(t, input.FinalTime, 4.)
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

func nearVecScalar(a []float64, b float64, tol float64) (l bool) {
	near := func(a, b float64, tolI ...float64) (l bool) {
		var (
			tol float64
		)
		if len(tolI) == 0 {
			tol = 1.e-08
		} else {
			tol = tolI[0]
		}
		bound := math.Max(tol, tol*math.Abs(a))
		val := math.Abs(a - b)
		if val <= bound {
			l = true
		} else {
			fmt.Printf("Diff = %v, Left = %v, Right = %v\n", val, a, b)
		}
		return
	}
	for i, val := range a {
		if !near(b, val, tol) {
			fmt.Printf("Diff = %v, Left[%d] = %v, Right[%d] = %v\n", math.Abs(val-b), i, val, i, b)
			return false
		}
	}
	return true
}

func InitializePolynomial(X, Y utils.Matrix) (Q [4]utils.Matrix) {
	var (
		Np, Kmax = X.Dims()
	)
	for n := 0; n < 4; n++ {
		Q[n] = utils.NewMatrix(Np, Kmax)
	}
	for ii := 0; ii < Np*Kmax; ii++ {
		x, y := X.DataP[ii], Y.DataP[ii]
		rho, rhoU, rhoV, E := GetStatePoly(x, y)
		Q[0].DataP[ii] = rho
		Q[1].DataP[ii] = rhoU
		Q[2].DataP[ii] = rhoV
		Q[3].DataP[ii] = E
	}
	return
}

func GetStatePoly(x, y float64) (rho, rhoU, rhoV, E float64) {
	/*
		Matlab script:
				syms a b c d x y gamma
				%2D Polynomial field
				rho=a*abs(x)+b*abs(y);
				u = c*x; v = d*y;
				rhou=rho*u; rhov=rho*v;
				p=rho^gamma;
				q=0.5*rho*(u^2+v^2);
				E=p/(gamma-1)+q;
				U = [ rho, rhou, rhov, E];
				F = [ rhou, rho*u^2+p, rho*u*v, u*(E+p) ];
				G = [ rhov, rho*u*v, rho*v^2+p, v*(E+p) ];
				div = diff(F,x)+diff(G,y);
				fprintf('Code for divergence of F and G FluxIndex\n%s\n',ccode(div));
				fprintf('Code for U \n%s\n%s\n%s\n%s\n',ccode(U));
	*/
	var (
		a, b, c, d = 1., 1., 1., 1.
		pow        = math.Pow
		fabs       = math.Abs
		gamma      = 1.4
	)
	rho = a*fabs(x) + b*fabs(y)
	rhoU = c * x * (a*fabs(x) + b*fabs(y))
	rhoV = d * y * (a*fabs(x) + b*fabs(y))
	E = ((c*c)*(x*x)+(d*d)*(y*y))*((a*fabs(x))/2.0+(b*fabs(y))/2.0) + pow(a*fabs(x)+b*fabs(y), gamma)/(gamma-1.0)
	return
}

func GetDivergencePoly(t, x, y float64) (div [4]float64) {
	var (
		gamma      = 1.4
		pow        = math.Pow
		fabs       = math.Abs
		a, b, c, d = 1., 1., 1., 1.
	)
	div[0] = c*(a*fabs(x)+b*fabs(y)) + d*(a*fabs(x)+b*fabs(y)) + a*c*x*(x/fabs(x)) + b*d*y*(y/fabs(y))
	div[1] = (c*c)*x*(a*fabs(x)+b*fabs(y))*2.0 + c*d*x*(a*fabs(x)+b*fabs(y)) + a*(c*c)*(x*x)*(x/fabs(x)) + a*gamma*(x/fabs(x))*pow(a*fabs(x)+b*fabs(y), gamma-1.0) + b*c*d*x*y*(y/fabs(y))
	div[2] = (d*d)*y*(a*fabs(x)+b*fabs(y))*2.0 + c*d*y*(a*fabs(x)+b*fabs(y)) + b*(d*d)*(y*y)*(y/fabs(y)) + b*gamma*(y/fabs(y))*pow(a*fabs(x)+b*fabs(y), gamma-1.0) + a*c*d*x*y*(x/fabs(x))
	div[3] = c*(((c*c)*(x*x)+(d*d)*(y*y))*((a*fabs(x))/2.0+(b*fabs(y))/2.0)+pow(a*fabs(x)+b*fabs(y), gamma)+pow(a*fabs(x)+b*fabs(y), gamma)/(gamma-1.0)) + d*(((c*c)*(x*x)+(d*d)*(y*y))*((a*fabs(x))/2.0+(b*fabs(y))/2.0)+pow(a*fabs(x)+b*fabs(y), gamma)+pow(a*fabs(x)+b*fabs(y), gamma)/(gamma-1.0)) + c*x*((c*c)*x*((a*fabs(x))/2.0+(b*fabs(y))/2.0)*2.0+(a*(x/fabs(x))*((c*c)*(x*x)+(d*d)*(y*y)))/2.0+a*gamma*(x/fabs(x))*pow(a*fabs(x)+b*fabs(y), gamma-1.0)+(a*gamma*(x/fabs(x))*pow(a*fabs(x)+b*fabs(y), gamma-1.0))/(gamma-1.0)) + d*y*((b*(y/fabs(y))*((c*c)*(x*x)+(d*d)*(y*y)))/2.0+(d*d)*y*((a*fabs(x))/2.0+(b*fabs(y))/2.0)*2.0+b*gamma*(y/fabs(y))*pow(a*fabs(x)+b*fabs(y), gamma-1.0)+(b*gamma*(y/fabs(y))*pow(a*fabs(x)+b*fabs(y), gamma-1.0))/(gamma-1.0))
	return
}

func FluxCalcMomentumOnly(rho, rhoU, rhoV, E float64) (Fx, Fy [4]float64) {
	Fx, Fy =
		[4]float64{rhoU, rhoU, rhoU, rhoU},
		[4]float64{rhoV, rhoV, rhoV, rhoV}
	return
}

func CheckFlux0(c *Euler, t *testing.T) {
	/*
	   		Conditions of this test:
	            - Two duplicated triangles, removes the question of transformation Jacobian making the results differ
	            - Flux is calculated identically for each equation (only density components), removes the question of flux
	              accuracy being different for the more complex equations
	            - Flowfield is initialized to a freestream for a polynomial field, interpolation to edges is not done,
	              instead, analytic initialization values are put into the edges
	             Result:
	            - No test of different triangle shapes and orientations
	            - No test of accuracy of interpolation to edges
	            - No accuracy test of the complex polynomial fluxes in Q[1-3]
	*/
	if c.Partitions.ParallelDegree != 1 {
		panic("parallel degree should be 1 for this test")
	}
	c.FluxCalcMock = FluxCalcMomentumOnly // For testing, only consider the first component of flux for all [4]
	// Initialize
	X, Y := c.DFR.FluxX, c.DFR.FluxY
	QFlux := InitializePolynomial(X, Y)
	Kmax := c.DFR.K
	Nint := c.DFR.FluxElement.NpInt
	Nedge := c.DFR.FluxElement.NpEdge
	NpFlux := c.DFR.FluxElement.Np
	var Q, Q_Face, F_RT_DOF [4]utils.Matrix
	var QMean [4]utils.Vector
	var Flux_Face [2][4]utils.Matrix
	for n := 0; n < 4; n++ {
		Q[n] = utils.NewMatrix(Nint, Kmax)
		Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
		Flux_Face[0][n] = utils.NewMatrix(3*Nedge, Kmax)
		Flux_Face[1][n] = utils.NewMatrix(3*Nedge, Kmax)
		F_RT_DOF[n] = utils.NewMatrix(NpFlux, Kmax)
		QMean[n] = utils.NewVector(Kmax)
		for k := 0; k < Kmax; k++ {
			for i := 0; i < Nint; i++ {
				ind := k + i*Kmax
				Q[n].DataP[ind] = QFlux[n].DataP[ind]
			}
			for i := 0; i < 3*Nedge; i++ {
				ind := k + i*Kmax
				ind2 := k + (i+2*Nint)*Kmax
				Q_Face[n].DataP[ind] = QFlux[n].DataP[ind2]
			}
		}
	}
	c.SetRTFluxInternal(Kmax, c.DFR.Jdet, c.DFR.Jinv, F_RT_DOF, Q)
	EdgeQ1 := make([][4]float64, Nedge)
	EdgeQ2 := make([][4]float64, Nedge)
	c.UpdateElementMean(Q, QMean)
	// No need to interpolate to the edges, they are left at initialized state in Q_Face
	c.CalculateEdgeEulerFlux(0, [][4]utils.Matrix{Q_Face},
		[][4]utils.Vector{QMean},
		[][2][4]utils.Matrix{Flux_Face}, EdgeQ1, EdgeQ2, c.SortedEdgeKeys[0])
	c.SetRTFluxOnEdges(0, Kmax, F_RT_DOF)

	var div utils.Matrix
	for n := 0; n < 4; n++ {
		div = c.DFR.FluxElement.DivInt.Mul(F_RT_DOF[n])
		d1, d2 := div.Dims()
		assert.Equal(t, d1, Nint)
		assert.Equal(t, d2, Kmax)
		for k := 0; k < Kmax; k++ {
			Jdet := c.DFR.Jdet.At(k, 0)
			for i := 0; i < Nint; i++ {
				ind := k + i*Kmax
				div.DataP[ind] /= Jdet
			}
		}
		// Get the analytic values of divergence for comparison
		nn := 0 // Use only the density component of divergence to check
		for k := 0; k < Kmax; k++ {
			for i := 0; i < Nint; i++ {
				ind := k + i*Kmax
				x, y := X.DataP[ind], Y.DataP[ind]
				divC := GetDivergencePoly(0, x, y)
				divCalc := div.DataP[ind]
				normalizer := Q[nn].DataP[ind]
				// test := near(divCalc/normalizer, divC[nn]/normalizer, 0.0001) // 1% of field value
				assert.InDeltaf(t, divCalc/normalizer, divC[nn]/normalizer, 0.0001, "err msg %s") // 1% of field value
			}
		}
	}
}
