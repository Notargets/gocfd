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

	// Test construction
	// For each basis polynomial (column), evaluate at appropriate positions
	// Interior block
	// R side of interior first
	for j := 0; j < rt.NpInt; j++ {
		basis := rt.RTPolyBasis2D_A
		PSI = rt.GetAllPolynomials()
		for i := 0; i < rt.NpInt; i++ {
			r, s := rt.R.AtVec(i), rt.S.AtVec(i)
			rt.A.Set(i, j)
		}
	}
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
