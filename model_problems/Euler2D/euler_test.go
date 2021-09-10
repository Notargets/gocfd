package Euler2D

import (
	"fmt"
	"math"
	"strconv"
	"testing"

	"github.com/notargets/gocfd/types"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

func TestEuler(t *testing.T) {
	var (
		msg = "err msg %s"
		tol = 0.000001
	)
	if true {
		{ // Test interpolation of solution to edges for all supported orders
			Nmax := 7
			for N := 1; N <= Nmax; N++ {
				c := NewEuler(1, N, "../../DG2D/test_tris_5.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, false, false, false)
				Kmax := c.dfr.K
				Nint := c.dfr.FluxElement.Nint
				Nedge := c.dfr.FluxElement.Nedge
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
					Q_Face[n] = c.dfr.FluxEdgeInterpMatrix.Mul(Q[n])
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
				0a) Flux is computed and stored in X, Y component projections in the 2*Nint front of F_RT_DOF
				1) Solution is extrapolated to edge points in Q_Face from Q
				2) Edges are traversed, flux is calculated and projected onto edge face normals, scaled and placed into F_RT_DOF
			*/
			Nmax := 7
			for N := 1; N <= Nmax; N++ {
				c := NewEuler(1, N, "../../DG2D/test_tris_5.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, false, false, false)
				Kmax := c.dfr.K
				Nint := c.dfr.FluxElement.Nint
				Nedge := c.dfr.FluxElement.Nedge
				NpFlux := c.dfr.FluxElement.Np // Np = 2*Nint+3*Nedge
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
					Fr_check[k], Fs_check[k] = c.CalculateFluxTransformed(k, Kmax, 0, c.dfr.Jdet, c.dfr.Jinv, Q)
				}
				// Interpolate from solution points to edges using precomputed interpolation matrix
				for n := 0; n < 4; n++ {
					Q_Face[n] = c.dfr.FluxEdgeInterpMatrix.Mul(Q[n])
				}
				// Calculate flux and project into R and S (transformed) directions
				rtD := [4][]float64{F_RT_DOF[0].DataP, F_RT_DOF[1].DataP, F_RT_DOF[2].DataP, F_RT_DOF[3].DataP}
				for n := 0; n < 4; n++ {
					for i := 0; i < Nint; i++ {
						for k := 0; k < Kmax; k++ {
							ind := k + i*Kmax
							Fr, Fs := c.CalculateFluxTransformed(k, Kmax, i, c.dfr.Jdet, c.dfr.Jinv, Q)
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
							Fr, Fs := c.CalculateFluxTransformed(k, Kmax, i, c.dfr.Jdet, c.dfr.Jinv, Q_Face)
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
			for N := 1; N <= Nmax; N++ {
				c := NewEuler(1, N, "../../DG2D/test_tris_5.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, false, false, false)
				Kmax := c.dfr.K
				Nint := c.dfr.FluxElement.Nint
				Nedge := c.dfr.FluxElement.Nedge
				NpFlux := c.dfr.FluxElement.Np // Np = 2*Nint+3*Nedge
				// Mark the initial state with the element number
				var Q_Face, F_RT_DOF [4]utils.Matrix
				for n := 0; n < 4; n++ {
					F_RT_DOF[n] = utils.NewMatrix(NpFlux, Kmax)
					Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
				}
				Q := c.Q[0]
				c.SetNormalFluxInternal(Kmax, c.dfr.Jdet, c.dfr.Jinv, F_RT_DOF, Q)
				c.InterpolateSolutionToEdges(Q, Q_Face)
				EdgeQ1, EdgeQ2 := make([][4]float64, Nedge), make([][4]float64, Nedge)
				c.SetNormalFluxOnEdges(0, false, nil, nil, [][4]utils.Matrix{F_RT_DOF}, [][4]utils.Matrix{Q_Face}, c.SortedEdgeKeys[0], EdgeQ1, EdgeQ2)
				// Check that freestream divergence on this mesh is zero
				for n := 0; n < 4; n++ {
					var div utils.Matrix
					div = c.dfr.FluxElement.DivInt.Mul(F_RT_DOF[n])
					for k := 0; k < Kmax; k++ {
						for i := 0; i < Nint; i++ {
							ind := k + i*Kmax
							div.DataP[ind] /= c.dfr.Jdet.DataP[k]
						}
					}
					assert.True(t, nearVecScalar(div.DataP, 0., 0.000001))
				}
			}
		}
		{ // Test divergence of polynomial initial condition against analytic values
			/*
				Note: the Polynomial flux is asymmetric around the X and Y axes - it uses abs(x) and abs(y)
				Elements should not straddle the axes if a perfect polynomial flux capture is needed
			*/
			Nmax := 7
			for N := 1; N <= Nmax; N++ {
				plotMesh := false
				// Single triangle test case
				var c *Euler
				c = NewEuler(1, N, "../../DG2D/test_tris_1tri.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, plotMesh, false, false)
				CheckFlux0(c, t)
				// Two widely separated triangles - no shared faces
				c = NewEuler(1, N, "../../DG2D/test_tris_two.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, plotMesh, false, false)
				CheckFlux0(c, t)
				// Two widely separated triangles - no shared faces - one tri listed in reverse order
				c = NewEuler(1, N, "../../DG2D/test_tris_twoR.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, plotMesh, false, false)
				CheckFlux0(c, t)
				// Connected tris, sharing one edge
				// plotMesh = true
				c = NewEuler(1, N, "../../DG2D/test_tris_6_nowall.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, plotMesh, false, false)
				CheckFlux0(c, t)
			}
		}
	}
	if true {
		{ // Test divergence of Isentropic Vortex initial condition against analytic values - density equation only
			N := 1
			plotMesh := false
			c := NewEuler(1, N, "../../DG2D/test_tris_6.neu", 1, FLUX_Average, IVORTEX, 1, 0, 1.4, 0, false, 5000, plotMesh, false, false)
			for _, e := range c.dfr.Tris.Edges {
				if e.BCType == types.BC_IVortex {
					e.BCType = types.BC_None
				}
			}
			Kmax := c.dfr.K
			Nint := c.dfr.FluxElement.Nint
			Nedge := c.dfr.FluxElement.Nedge
			NpFlux := c.dfr.FluxElement.Np // Np = 2*Nint+3*Nedge
			// Mark the initial state with the element number
			var Q_Face, F_RT_DOF [4]utils.Matrix
			for n := 0; n < 4; n++ {
				F_RT_DOF[n] = utils.NewMatrix(NpFlux, Kmax)
				Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
			}
			Q := c.Q[0]
			X, Y := c.dfr.FluxX, c.dfr.FluxY
			c.SetNormalFluxInternal(Kmax, c.dfr.Jdet, c.dfr.Jinv, F_RT_DOF, Q)
			c.InterpolateSolutionToEdges(Q, Q_Face)
			EdgeQ1, EdgeQ2 := make([][4]float64, Nedge), make([][4]float64, Nedge)
			c.SetNormalFluxOnEdges(0, false, nil, nil, [][4]utils.Matrix{F_RT_DOF}, [][4]utils.Matrix{Q_Face}, c.SortedEdgeKeys[0], EdgeQ1, EdgeQ2)
			var div utils.Matrix
			// Density is the easiest equation to match with a polynomial
			n := 0
			fmt.Printf("component[%d]\n", n)
			div = c.dfr.FluxElement.DivInt.Mul(F_RT_DOF[n])
			//c.DivideByJacobian(Kmax, Nint, c.dfr.Jdet, div.DataP, 1)
			for k := 0; k < Kmax; k++ {
				for i := 0; i < Nint; i++ {
					ind := k + i*Kmax
					div.DataP[ind] /= -c.dfr.Jdet.DataP[k]
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

func TestFluxJacobian(t *testing.T) {
	var (
		tol = 0.000001
		msg = "err msg %s"
	)
	{ // Flux Jacobian calculation
		c := Euler{}
		c.FS = NewFreeStream(0.1, 1.4, 0)
		Qinf := c.FS.Qinf
		Fx, Gy := c.FluxJacobianCalc(Qinf[0], Qinf[1], Qinf[2], Qinf[3])
		// Matlab: using the FS Q = [1,.1,0,1.79071]
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
	// Test build system matrix (for implicit solver)
	{
		N := 1
		c := NewEuler(1, N, "../../DG2D/test_tris_1tri.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, false, false, false)
		ei := c.NewElementImplicit()
		var (
			myThread         = 0
			Q0               = c.Q[myThread]
			Kmax, Jdet, Jinv = ei.Kmax[myThread], ei.Jdet[myThread], ei.Jinv[myThread]
			Q_Face, F_RT_DOF = ei.Q_Face[myThread], ei.F_RT_DOF[myThread]
			FluxJac          = ei.FluxJac[myThread]
		)
		c.PrepareEdgeFlux(Kmax, Jdet, Jinv, F_RT_DOF, Q0, Q_Face)
		c.SetFluxJacobian(Kmax, Jdet, Jinv, Q0, Q_Face, FluxJac)
		// Compose system matrix, one for each element
	}
	// Test implicit solver
	{
		N := 1
		c := NewEuler(1, N, "../../DG2D/test_tris_1tri.neu", 1, FLUX_Average, FREESTREAM, 1, 0, 1.4, 0, false, 5000, false, false, false)
		ei := c.NewElementImplicit()
		var (
			myThread         = 0
			Q0               = c.Q[myThread]
			Kmax, Jdet, Jinv = ei.Kmax[myThread], ei.Jdet[myThread], ei.Jinv[myThread]
			Q_Face, F_RT_DOF = ei.Q_Face[myThread], ei.F_RT_DOF[myThread]
			FluxJac          = ei.FluxJac[myThread]
			RHSQ, Residual   = ei.RHSQ[myThread], ei.Residual[myThread]
			Nedge            = ei.Nedge
			DT               = ei.DT[myThread]
			SM               = ei.SM[myThread]
			EdgeQ1, EdgeQ2   = make([][4]float64, Nedge), make([][4]float64, Nedge) // Local working memory
		)
		c.PrepareEdgeFlux(Kmax, Jdet, Jinv, F_RT_DOF, Q0, Q_Face)
		c.SetNormalFluxOnEdges(ei.Time, true, ei.Jdet, ei.DT, ei.F_RT_DOF, ei.Q_Face, c.SortedEdgeKeys[myThread], EdgeQ1, EdgeQ2) // Global
		// Create the RHS for all flux points in every element
		c.RHSInternalPoints(Kmax, Jdet, F_RT_DOF, RHSQ)
		c.SetFluxJacobian(Kmax, Jdet, Jinv, Q0, Q_Face, FluxJac)
		ei.SolveForResidual(Kmax, c.LocalTimeStepping, DT, Jdet, c.dfr.FluxElement.DivInt, SM, RHSQ, Residual, FluxJac)
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
				fprintf('Code for Divergence of F and G Fluxes\n%s\n',ccode(div));
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
	X, Y := c.dfr.FluxX, c.dfr.FluxY
	QFlux := InitializePolynomial(X, Y)
	Kmax := c.dfr.K
	Nint := c.dfr.FluxElement.Nint
	Nedge := c.dfr.FluxElement.Nedge
	NpFlux := c.dfr.FluxElement.Np
	var Q, Q_Face, F_RT_DOF [4]utils.Matrix
	for n := 0; n < 4; n++ {
		Q[n] = utils.NewMatrix(Nint, Kmax)
		Q_Face[n] = utils.NewMatrix(3*Nedge, Kmax)
		F_RT_DOF[n] = utils.NewMatrix(NpFlux, Kmax)
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
	c.SetNormalFluxInternal(Kmax, c.dfr.Jdet, c.dfr.Jinv, F_RT_DOF, Q)
	EdgeQ1, EdgeQ2 := make([][4]float64, Nedge), make([][4]float64, Nedge)
	// No need to interpolate to the edges, they are left at initialized state in Q_Face
	c.SetNormalFluxOnEdges(0, false, nil, nil, [][4]utils.Matrix{F_RT_DOF}, [][4]utils.Matrix{Q_Face}, c.SortedEdgeKeys[0], EdgeQ1, EdgeQ2)

	var div utils.Matrix
	for n := 0; n < 4; n++ {
		div = c.dfr.FluxElement.DivInt.Mul(F_RT_DOF[n])
		d1, d2 := div.Dims()
		assert.Equal(t, d1, Nint)
		assert.Equal(t, d2, Kmax)
		for k := 0; k < Kmax; k++ {
			Jdet := c.dfr.Jdet.At(k, 0)
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
				//test := near(divCalc/normalizer, divC[nn]/normalizer, 0.0001) // 1% of field value
				assert.InDeltaf(t, divCalc/normalizer, divC[nn]/normalizer, 0.0001, "err msg %s") // 1% of field value
			}
		}
	}
}
