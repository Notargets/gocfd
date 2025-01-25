package DG2D

import (
	"fmt"
	"math"
	"strconv"
	"testing"
	"time"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"

	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"

	"github.com/stretchr/testify/assert"
)

func TestRTElementConstruction(t *testing.T) {
	// Define an RT element at order P
	P := 2
	rt := NewRTElement(P)

	Print := func(label string, iip *int) {
		fmt.Printf("%s", label)
		for i := 0; i < rt.NpEdge; i++ {
			fmt.Printf("[%f,%f] ", rt.R.DataP[*iip], rt.S.DataP[*iip])
			*iip++
		}
		fmt.Printf("\n")
	}

	ii := 2 * rt.NpInt
	Print("Edge 1:", &ii)
	Print("Edge 2:", &ii)
	Print("Edge 3:", &ii)

	// Test polynomial bases

	// TEST1: Demonstrate each of the ways to calculate Polynomial Terms
	// 2D Interior Polynomials
	PInt := rt.P - 1
	RInt, SInt := NodesEpsilon(PInt)
	BasisTest := func(basis *JacobiBasis2D) (PSI utils.Vector, P_Alt []float64) {
		PSI = basis.GetAllPolynomials()
		CalcTerm := func(r, s float64, P int) (psi float64) {
			for i := 0; i <= P; i++ {
				for j := 0; j <= (P - i); j++ {
					pTerm := basis.PolynomialTerm(r, s, i, j)
					psi += pTerm
				}
			}
			return
		}
		r1, s1 := RInt.AtVec(0), SInt.AtVec(0)
		r2, s2 := RInt.AtVec(1), SInt.AtVec(1)
		r3, s3 := RInt.AtVec(2), SInt.AtVec(2)
		P_Alt = []float64{
			CalcTerm(r1, s1, PInt),
			CalcTerm(r2, s2, PInt),
			CalcTerm(r3, s3, PInt),
		}
		for i := 0; i < len(P_Alt); i++ {
			fmt.Printf("P_Alt[%d] = %f\n", i, P_Alt[i])
		}
		return
	}
	// R direction basis
	PSI, P_Alt := BasisTest(rt.RTPolyBasis2D_A)
	P1 := PSI.AtVec(0)
	assert.True(t, near(-.264298, P1, 0.00001))
	for i := 0; i < rt.NpInt; i++ {
		assert.True(t, nearVec(P_Alt, PSI.DataP, 0.00001))
	}

	// S direction basis
	PSI, P_Alt = BasisTest(rt.RTPolyBasis2D_B)
	P1 = PSI.AtVec(0)
	assert.True(t, near(-3.033715, P1, 0.00001))
	for i := 0; i < rt.NpInt; i++ {
		assert.True(t, nearVec(P_Alt, PSI.DataP, 0.00001))
	}

	// 1D Edge Polynomials
	// Get the edge values for edge1,2,3
	assert.Panics(t, func() { rt.getEdgeCoordinates(0) })
	// The edge distribution of edge2 should be the reverse direction of
	// edge1, and symmetric, shorthand is to take the neg of edge2 to get edge1
	edge1 := rt.getEdgeCoordinates(2)
	for i := range edge1.DataP {
		edge1.DataP[i] *= -1
	}
	assert.True(t, nearVec(rt.getEdgeCoordinates(1).DataP,
		edge1.DataP, 0.000001))
	assert.True(t, nearVec(DG1D.JacobiP(edge1, 0, 0, rt.P),
		rt.RTPolyBasis1D_Edge1, 0.00001))
	edge2 := rt.getEdgeCoordinates(2)
	assert.True(t, nearVec(DG1D.JacobiP(edge2, 0, 0, rt.P),
		rt.RTPolyBasis1D_Edge2, 0.00001))
	edge3 := rt.getEdgeCoordinates(3)
	assert.True(t, nearVec(DG1D.JacobiP(edge3, 0, 0, rt.P),
		rt.RTPolyBasis1D_Edge3, 0.00001))
}

/*
This follows the development in the V.J. Ervin paper: Computational Bases for
RTk and BDMk on Triangles
All are defined in terms of Eta and Xi which go from 0 to 1 compared to here
where R and S are defined from -1 to 1

	xi = (s + 1)/2.
	eta = (r + 1)/2.

There are three "Normal" or "edge" functions and on each of the 3 triangle edges
there are some number of points with a single vector basis function.

For the RT2 element, there are 3 points on each of the 3 edges. The vectors are
defined below such that each column is one edge, having three rows, each
defining the basis function for that point on the edge:
Φ1 (ξ, η) = q1 (η) ê1 (ξ, η) , Φ2 (ξ, η) = q2 (η) ê1 (ξ, η) , Φ3 (ξ, η) = q3 (η) ê1 (ξ, η) ,
Φ1 (ξ, η) = q3 (η) ê2 (ξ, η) , Φ2 (ξ, η) = q2 (η) ê2 (ξ, η) , Φ3 (ξ, η) = q1 (η) ê2 (ξ, η) ,
Φ1 (ξ, η) = q1 (ξ) ê3 (ξ, η) , Φ2 (ξ, η) = q2 (ξ) ê3 (ξ, η) , Φ3 (ξ, η) = q3 (ξ) ê3 (ξ, η) .

For the RT2 element, there are 3 interior points, and in the below example, we
see there are two basis functions (vectors) for each of the 3 points. Each
point is a column, and each column has two rows defining the basis vectors for
that point:
Φ1 (ξ, η) = (1 − ξ − η) ê4 (ξ, η) , Φ2 (ξ, η) = ξ ê4 (ξ, η) , Φ3 (ξ, η) = η ê4 (ξ, η) ,
Φ1 (ξ, η) = (1 − ξ − η) ê5 (ξ, η) , Φ2 (ξ, η) = ξ ê5 (ξ, η) , Φ3 (ξ, η) = η ê5 (ξ, η) .

Given that the edge functions are only evaluated on the edges, and that the
lagrange functions multiplying them cause them to vanish on all but their
defining point, it seems unnecessary to evaluate anything but the e1-e3
functions in practice.

For the general case of the interior basis functions, they take the form:
Φj (ξ, η) = bj (ξ, η) ê4 (ξ, η) , j = 1, 2, . . . , k(k + 1)/2 ,
Φj (ξ, η) = bj (ξ, η) ê5 (ξ, η) , j = 1, 2, . . . , k(k + 1)/2 .

where bj is a 2D polynomial of order P-1
*/
func baseBasisFunctions(r, s float64, fNum int) (ef [2]float64) {
	// Mapping:
	// Ervin Edge 1 => Edge 2 (Hypotenuse)
	// Ervin Edge 2 => Edge 1 (Left edge)
	// Ervin Edge 3 => Edge 3 (Bottom edge)
	// Note: Here we fix an error in Ervin - the midpoint of RT0 edges are now
	// correctly unit normal
	var (
		sr2 = math.Sqrt(2.)
		// Transform our unit triangle coords to Ervin:
		xi  = 0.5 * (r + 1)
		eta = 0.5 * (s + 1)
	)
	switch fNum {
	case 2:
		ef = [2]float64{sr2 * xi, sr2 * eta}
	case 1:
		ef = [2]float64{xi - 1, eta - 0.5}
	case 3:
		ef = [2]float64{xi - 0.5, eta - 1}
	case 4:
		ef = [2]float64{eta * xi, eta * (eta - 1)}
	case 5:
		ef = [2]float64{xi * (xi - 1), xi * eta}
	default:
		panic("wrong basis function number (1-5)")
	}
	return
}

func TestErvinBasisFunctions2(t *testing.T) {
	R := []float64{1. / 3., 0.5, 2. / 3.}
	assert.Equal(t, 1., DG1D.Lagrange1DPoly(1./3., R, 0))
	assert.Equal(t, 0., DG1D.Lagrange1DPoly(1./3., R, 1))
	assert.Equal(t, 0., DG1D.Lagrange1DPoly(1./3., R, 2))
	assert.Panics(t, func() { DG1D.Lagrange1DPoly(1./3., R, 3) })
	assert.InDeltaf(t, -9., DG1D.Lagrange1DPoly(1./3., R, 0, 1), 0.000001, "")
	assert.InDeltaf(t, 12., DG1D.Lagrange1DPoly(1./3., R, 1, 1), 0.000001, "")
	assert.InDeltaf(t, -3., DG1D.Lagrange1DPoly(1./3., R, 2, 1), 0.000001, "")

	// Generate Gauss Lobato points for P=5 to compare with the online article:
	// https://math.stackexchange.com/questions/1105160/evaluate-derivative-of-lagrange-polynomials-at-construction-points
	R = DG1D.JacobiGL(0, 0, 6).DataP
	// One row (i) is the evaluation of the j-th derivative at each i-th point
	validation_deriv := make([][]float64, len(R))
	for i := range validation_deriv {
		validation_deriv[i] = make([]float64, len(R))
	}
	validation_deriv[0] = []float64{-10.5, -2.4429, 0.6253, -0.3125, 0.2261,
		-0.2266, 0.5}
	validation_deriv[1] = []float64{14.2016, 0, -2.2158, 0.9075, -0.6164,
		0.6022, -1.3174}
	validation_deriv[2] = []float64{-5.669, 3.4558, 0, -2.007, 1.0664, -0.9613,
		2.05}
	validation_deriv[3] = []float64{3.2, -1.5986, 2.2667, 0, -2.2667, 1.5986,
		-3.2}
	validation_deriv[4] = []float64{-2.05, 0.9613, -1.0664, 2.007, 0, -3.4558,
		5.669}
	validation_deriv[5] = []float64{1.3174, -0.6022, 0.6164, -0.9075, 2.2158, 0,
		-14.2016}
	validation_deriv[6] = []float64{-0.5, 0.2266, -0.2261, 0.3125, -0.6253,
		2.4429, 10.5}
	testVec := make([]float64, len(R))
	for j, validate := range validation_deriv {
		for i, r := range R {
			testVec[i] = DG1D.Lagrange1DPoly(r, R, j, 1)
		}
		assert.True(t, nearVec(testVec, validate, 0.0001))
	}

}
func TestErvinBasisFunctions1(t *testing.T) {
	// This tests the basic basis functions e1,e2,e3 for edges and e4,
	// e5 interior
	var (
		sr2   = math.Sqrt(2.)
		oosr2 = 1. / sr2
	)
	conv := func(r, s float64) (xiEta [2]float64) {
		xiEta = [2]float64{0.5 * (r + 1), 0.5 * (s + 1)}
		return
	}
	assert.Equal(t, conv(-1, -1), [2]float64{0., 0.})
	assert.Equal(t, conv(-1, 1), [2]float64{0., 1.})
	assert.Equal(t, conv(1, -1), [2]float64{1., 0.})
	// Test the midpoint of each edge for values of the base edge functions
	// Edge midpoints
	ef1 := baseBasisFunctions(-1, 0, 1)
	ef2 := baseBasisFunctions(0, 0, 2)
	ef3 := baseBasisFunctions(0, -1, 3)
	assert.Equal(t, [2]float64{-1., 0}, ef1)
	assert.InDeltaf(t, oosr2, ef2[0], 0.0000001, "")
	assert.InDeltaf(t, oosr2, ef2[1], 0.0000001, "")
	assert.Equal(t, [2]float64{0, -1}, ef3)

	// Test the dot product of the e4 and e5 interior basis vectors against all
	// edge normals at the edge normal locations. All interior basis functions
	// should have zero dot product on all edges
	dot := func(v1, v2 [2]float64) (dp float64) {
		dp = v1[0]*v2[0] + v1[1]*v2[1]
		return
	}
	// Edge 1 midpoint
	ef4 := baseBasisFunctions(-1, 0, 4)
	ef5 := baseBasisFunctions(-1, 0, 5)
	assert.Equal(t, 0., dot(ef4, ef1))
	// Edge 1 endpoint 1
	ef4 = baseBasisFunctions(-1, -1, 4)
	ef5 = baseBasisFunctions(-1, -1, 5)
	assert.Equal(t, 0., dot(ef4, ef1))
	assert.Equal(t, 0., dot(ef5, ef1))
	// Edge 1 endpoint 2
	ef4 = baseBasisFunctions(-1, 1, 4)
	ef5 = baseBasisFunctions(-1, 1, 5)
	assert.Equal(t, 0., dot(ef4, ef1))
	assert.Equal(t, 0., dot(ef5, ef1))

	// Edge 2 midpoint
	ef4 = baseBasisFunctions(0, 0, 4)
	ef5 = baseBasisFunctions(0, 0, 5)
	assert.Equal(t, 0., dot(ef4, ef2))
	assert.Equal(t, 0., dot(ef5, ef2))
	// Edge 2 endpoint 1
	ef4 = baseBasisFunctions(1, -1, 4)
	ef5 = baseBasisFunctions(1, -1, 5)
	assert.Equal(t, 0., dot(ef4, ef2))
	assert.Equal(t, 0., dot(ef5, ef2))
	// Edge 2 endpoint 2
	ef4 = baseBasisFunctions(-1, 1, 4)
	ef5 = baseBasisFunctions(-1, 1, 5)
	assert.Equal(t, 0., dot(ef4, ef2))
	assert.Equal(t, 0., dot(ef5, ef2))

	// Edge 3 midpoint
	ef4 = baseBasisFunctions(0, -1, 4)
	ef5 = baseBasisFunctions(0, -1, 5)
	assert.Equal(t, 0., dot(ef4, ef3))
	assert.Equal(t, 0., dot(ef5, ef3))
	// Edge 3 endpoint 1
	ef4 = baseBasisFunctions(-1, -1, 4)
	ef5 = baseBasisFunctions(-1, -1, 5)
	assert.Equal(t, 0., dot(ef4, ef3))
	assert.Equal(t, 0., dot(ef5, ef3))
	// Edge 3 endpoint 2
	ef4 = baseBasisFunctions(1, -1, 4)
	ef5 = baseBasisFunctions(1, -1, 5)
	assert.Equal(t, 0., dot(ef4, ef3))
	assert.Equal(t, 0., dot(ef5, ef3))
}

func TestLagrangePolynomial(t *testing.T) {
	numSamples := 100
	rd := make([]float64, numSamples)
	xmin, xmax := -1., 1.
	fmin, fmax := -0.5, 1.25
	inc := (xmax - xmin) / float64(numSamples-1.)
	for i := 0; i < numSamples; i++ {
		rd[i] = xmin + float64(i)*inc
	}
	SamplesR := utils.NewVector(numSamples, rd)
	f := make([]float64, SamplesR.Len())
	var plot bool
	plot = false
	if plot {
		chart := utils.NewLineChart(1920, 1080, xmin, xmax, fmin, fmax)
		// TODO: Make a pluggable basis underneath the RT (and Lagrange) elements - Lagrange, Hesthaven, Spectral?
		var delay time.Duration
		Nmax := 4
		lineInc := 2. / float64(Nmax-2)
		lineColor := -1. // colormap goes from -1,1
		for n := 1; n < Nmax; n++ {
			Np := n + 1
			R := utils.NewVector(Np)
			inc = (xmax - xmin) / float64(Np-1)
			for i := 0; i < Np; i++ {
				R.DataP[i] = xmin + float64(i)*inc
			}
			lp := NewLagrangeBasis1D(R.DataP)
			_ = lp
			for j := 0; j < R.Len(); j++ {
				for i, r := range SamplesR.DataP {
					// f[i] = DG1D.JacobiP(utils.NewVector(1, []float64{r}), 0, 0, j)[0]
					f[i] = lp.BasisPolynomial([]float64{r}, j)[0]
				}
				if n == Nmax-1 && j == R.Len()-1 {
					delay = 120 * time.Second
				}
				name := "JacobiP[" + strconv.Itoa(n) + "," + strconv.Itoa(j) + "]"
				fmt.Printf("Chart Name: [%s], lineColor = %5.3f\n", name, lineColor)
				chart.Plot(delay, SamplesR.DataP, f, lineColor, name)
			}
			lineColor += lineInc
		}
	}

}

func TestRTElementLagrange(t *testing.T) {
	// Test 2D Lagrange polynomial basis
	if true {
		Nmax := 2
		Np := (Nmax + 1) * (Nmax + 2) / 2
		R, S := NodesEpsilon(Nmax)
		lb2d := NewLagrangeBasis2D(Nmax, utils.NewVector(Np, R.DataP), utils.NewVector(Np, S.DataP))
		RR := utils.NewVector(4, []float64{-1, -0.5, 0.5, 1.})
		SS := utils.NewVector(4, []float64{-1, -0.5, 0.5, 1.})
		assert.InDeltaSlicef(t, lb2d.BasisPolynomial(RR, SS, 0, 0), []float64{
			1.8736592735117479, 0.08527444844638511, 2.157995170376074, 0.1385595874119674,
		}, 0.0000001, "blah")
		Interp := lb2d.GetInterpMatrix(RR, SS)
		assert.InDeltaSlicef(t, Interp.SumRows().DataP, []float64{1, 1, 1, 1}, 0.000001, "blah")
		InterpDR, InterpDS := lb2d.GetGradInterpMatrices(RR, SS)
		assert.InDeltaSlicef(t, InterpDR.SumRows().DataP, []float64{0, 0, 0, 0}, 0.0000001, "blah")
		assert.InDeltaSlicef(t, InterpDS.SumRows().DataP, []float64{0, 0, 0, 0}, 0.0000001, "blah")
	}
	errorCheck := func(N int, div, divCheck []float64) (minInt, maxInt, minEdge, maxEdge float64) {
		var (
			Npm    = len(div)
			errors = make([]float64, Npm)
		)
		for i := 0; i < Npm; i++ {
			// var ddr, dds float64
			errors[i] = div[i] - divCheck[i]
		}
		minInt, maxInt = errors[0], errors[0]
		Nint := N * (N + 1) / 2
		minEdge, maxEdge = errors[Nint], errors[Nint]
		for i := 0; i < Nint; i++ {
			errAbs := math.Abs(errors[i])
			if minInt > errAbs {
				minInt = errAbs
			}
			if maxInt < errAbs {
				maxInt = errAbs
			}
		}
		for i := Nint; i < Npm; i++ {
			errAbs := math.Abs(errors[i])
			if minEdge > errAbs {
				minEdge = errAbs
			}
			if maxEdge < errAbs {
				maxEdge = errAbs
			}
		}
		fmt.Printf("Order = %d, ", N)
		fmt.Printf("Min, Max Int Err = %8.5f, %8.5f, Min, Max Edge Err = %8.5f, %8.5f\n", minInt, maxInt, minEdge, maxEdge)
		return
	}
	checkSolution := func(rt *RTElement, Order int) (s1, s2, divCheck []float64) {
		var (
			Np = rt.Np
		)
		s1, s2 = make([]float64, Np), make([]float64, Np)
		divCheck = make([]float64, Np)
		var ss1, ss2 float64
		for i := 0; i < Np; i++ {
			r := rt.R.DataP[i]
			s := rt.S.DataP[i]
			ccf := float64(Order)
			s1[i] = utils.POW(r, Order)
			s2[i] = utils.POW(s, Order)
			ss1, ss2 = ccf*utils.POW(r, Order-1), ccf*utils.POW(s, Order-1)
			divCheck[i] = ss1 + ss2
		}
		return
	}
	// Test the polynomial term function
	{
		lb2d := &LagrangeBasis2D{P: 0}
		assert.Equal(t, 0, lb2d.getTermNumber(0, 0))
		assert.Equal(t, -1, lb2d.getTermNumber(1, 0))
		assert.Equal(t, -1, lb2d.getTermNumber(0, 1))
	}
	if true { // Check Divergence for polynomial vector fields of order < N against analytical solution
		Nend := 8
		for N := 1; N < Nend; N++ {
			rt := NewRTElement(N)
			for cOrder := 0; cOrder <= N; cOrder++ {
				fmt.Printf("Check Order = %d, ", cOrder)
				// [s1,s2] values for each location in {R,S}
				s1, s2, divCheck := checkSolution(rt, cOrder)
				sp := rt.ProjectFunctionOntoDOF(s1, s2)
				sm := utils.NewMatrix(rt.Np, 1, sp)
				divM := rt.Div.Mul(sm)
				// fmt.Println(divM.Print("divM"))
				minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, divM.DataP, divCheck)
				assert.True(t, near(minerrInt, 0.0, 0.00001))
				assert.True(t, near(maxerrInt, 0.0, 0.00001))
				assert.True(t, near(minerrEdge, 0.0, 0.00001))
				assert.True(t, near(maxerrEdge, 0.0, 0.00001))
			}
		}
	}
}

func TestRTElement(t *testing.T) {
	{
		// Check term-wise orthogonal 2D polynomial basis
		N := 2
		R, S := NodesEpsilon(N - 1)
		JB2D := NewJacobiBasis2D(N-1, R, S, 0, 0)
		ii, jj := 1, 1
		p := JB2D.Simplex2DP(R, S, ii, jj)
		ddr, dds := JB2D.GradSimplex2DP(R, S, ii, jj)
		Np := R.Len()
		pCheck, ddrCheck, ddsCheck := make([]float64, Np), make([]float64, Np), make([]float64, Np)
		for i, rVal := range R.DataP {
			sVal := S.DataP[i]
			ddrCheck[i] = JB2D.PolynomialTermDr(rVal, sVal, ii, jj)
			ddsCheck[i] = JB2D.PolynomialTermDs(rVal, sVal, ii, jj)
			pCheck[i] = JB2D.PolynomialTerm(rVal, sVal, ii, jj)
		}
		assert.True(t, nearVec(pCheck, p, 0.000001))
		assert.True(t, nearVec(ddrCheck, ddr, 0.000001))
		assert.True(t, nearVec(ddsCheck, dds, 0.000001))
	}
	errorCheck := func(N int, div, divCheck []float64) (minInt, maxInt, minEdge, maxEdge float64) {
		var (
			Npm    = len(div)
			errors = make([]float64, Npm)
		)
		for i := 0; i < Npm; i++ {
			// var ddr, dds float64
			errors[i] = div[i] - divCheck[i]
		}
		minInt, maxInt = errors[0], errors[0]
		Nint := N * (N + 1) / 2
		minEdge, maxEdge = errors[Nint], errors[Nint]
		for i := 0; i < Nint; i++ {
			errAbs := math.Abs(errors[i])
			if minInt > errAbs {
				minInt = errAbs
			}
			if maxInt < errAbs {
				maxInt = errAbs
			}
		}
		for i := Nint; i < Npm; i++ {
			errAbs := math.Abs(errors[i])
			if minEdge > errAbs {
				minEdge = errAbs
			}
			if maxEdge < errAbs {
				maxEdge = errAbs
			}
		}
		fmt.Printf("Order = %d, ", N)
		fmt.Printf("Min, Max Int Err = %8.5f, %8.5f, Min, Max Edge Err = %8.5f, %8.5f\n", minInt, maxInt, minEdge, maxEdge)
		return
	}
	checkSolution := func(rt *RTElement, Order int) (s1, s2, divCheck []float64) {
		var (
			Np = rt.Np
		)
		s1, s2 = make([]float64, Np), make([]float64, Np)
		divCheck = make([]float64, Np)
		var ss1, ss2 float64
		for i := 0; i < Np; i++ {
			r := rt.R.DataP[i]
			s := rt.S.DataP[i]
			ccf := float64(Order)
			s1[i] = utils.POW(r, Order)
			s2[i] = utils.POW(s, Order)
			ss1, ss2 = ccf*utils.POW(r, Order-1), ccf*utils.POW(s, Order-1)
			divCheck[i] = ss1 + ss2
		}
		return
	}
	if true { // Check Divergence for polynomial vector fields of order < N against analytical solution
		Nend := 8
		for N := 1; N < Nend; N++ {
			rt := NewRTElement(N)
			for cOrder := 0; cOrder <= N; cOrder++ {
				fmt.Printf("Check Order = %d, ", cOrder)
				// [s1,s2] values for each location in {R,S}
				s1, s2, divCheck := checkSolution(rt, cOrder)
				sp := rt.ProjectFunctionOntoDOF(s1, s2)
				sm := utils.NewMatrix(rt.Np, 1, sp)
				divM := rt.Div.Mul(sm)
				// fmt.Println(divM.Print("divM"))
				minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, divM.DataP, divCheck)
				assert.True(t, near(minerrInt, 0.0, 0.00001))
				assert.True(t, near(maxerrInt, 0.0, 0.00001))
				assert.True(t, near(minerrEdge, 0.0, 0.00001))
				assert.True(t, near(maxerrEdge, 0.0, 0.00001))
			}
		}
	}
	plot := false
	if plot {
		N := 2
		rt := NewRTElement(N)
		s1, s2 := make([]float64, rt.R.Len()), make([]float64, rt.R.Len())
		for i := range rt.R.DataP {
			s1[i] = 1
			s2[i] = 1
		}
		if plot {
			chart := PlotTestTri(true)
			points := utils.ArraysToPoints(rt.R.DataP, rt.S.DataP)
			f := utils.ArraysTo2Vector(s1, s2, 0.1)
			_ = chart.AddVectors("test function", points, f, chart2d.Solid, utils.GetColor(utils.Green))
			utils.SleepFor(500000)
		}
	}
}

func PlotTestTri(plotGeom bool) (chart *chart2d.Chart2D) {
	var (
		points  []graphics2D.Point
		trimesh graphics2D.TriMesh
		K       = 1
	)

	points = make([]graphics2D.Point, 3)
	points[0].X[0], points[0].X[1] = -1, -1
	points[1].X[0], points[1].X[1] = 1, -1
	points[2].X[0], points[2].X[1] = -1, 1

	trimesh.Triangles = make([]graphics2D.Triangle, K)
	colorMap := utils2.NewColorMap(0, 1, 1)
	trimesh.Triangles[0].Nodes[0] = 0
	trimesh.Triangles[0].Nodes[1] = 1
	trimesh.Triangles[0].Nodes[2] = 2
	trimesh.Geometry = points
	box := graphics2D.NewBoundingBox(trimesh.GetGeometry())
	chart = chart2d.NewChart2D(1024, 1024, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
	chart.AddColorMap(colorMap)
	go chart.Plot()

	if plotGeom {
		if err := chart.AddTriMesh("TriMesh", trimesh,
			chart2d.CrossGlyph, 0.1, chart2d.Solid,
			utils.GetColor(utils.Black)); err != nil {
			panic("unable to add graph series")
		}
	}
	return
}

/*
	checkSolutionM := func(rt *RTElement, Order int) (s1, s2, divCheck []float64) {
		var (
			Npm = rt.Npm
		)
		s1, s2 = make([]float64, Npm), make([]float64, Npm)
		divCheck = make([]float64, Npm)
		var ss1, ss2 float64
		for i := 0; i < Npm; i++ {
			r := rt.Rm.DataP[i]
			s := rt.Sm.DataP[i]
			ccf := float64(Order)
			s1[i] = utils.POW(r, Order)
			s2[i] = utils.POW(s, Order)
			ss1, ss2 = ccf*utils.POW(r, Order-1), ccf*utils.POW(s, Order-1)
			divCheck[i] = ss1 + ss2
		}
		return
	}
	if false { // Check Divergence for polynomial vector fields of order < N against analytical solution
		Nend := 2
		for N := 1; N < Nend; N++ {
			R, S := NodesEpsilon(N - 1)
			rt := NewRTElement(N, R, S)
			fmt.Printf("calculating Dr and Ds...\n")
			//Dr := rt.Vrm[0].Mul(rt.VmInv[0])
			//Ds := rt.Vsm[1].Mul(rt.VmInv[1])
			Dr := rt.Vrm[0].Mul(rt.VmInv[0])
			Ds := rt.Vsm[1].Mul(rt.VmInv[1])
			//Dr := rt.Vrm[0]
			//Ds := rt.Vsm[1]
			for cOrder := 0; cOrder <= N; cOrder++ {
				fmt.Printf("Check Order = %d, ", cOrder)
				s1, s2, divCheck := checkSolutionM(rt, cOrder)
				s1p, s2p := rt.ProjectFunctionOntoBasis2(s1, s2)
				s1m, s2m := utils.NewMatrix(rt.Npm, 1, s1p), utils.NewMatrix(rt.Npm, 1, s2p)
				div := Dr.Mul(s1m).Add(Ds.Mul(s2m)).DataP
				minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, div, divCheck)
				assert.True(t, near(minerrInt, 0.0, 0.00001))
				assert.True(t, near(maxerrInt, 0.0, 0.00001))
				assert.True(t, near(minerrEdge, 0.0, 0.00001))
				assert.True(t, near(maxerrEdge, 0.0, 0.00001))
			}
		}
	}
*/
