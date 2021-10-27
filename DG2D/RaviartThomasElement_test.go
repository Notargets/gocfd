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

func TestLagrangePolynomial(t *testing.T) {
	{
		//N := 2
		//R, _ := NodesEpsilon(N - 1)
		numSamples := 5000
		rd := make([]float64, numSamples)
		rd[0] = -1
		inc := 2. / float64(numSamples-1.)
		for i := 1; i < numSamples; i++ {
			rd[i] = rd[0] + float64(i)*inc
		}
		// TODO: Make a pluggable basis underneath the RT (and Lagrange) elements - Lagrange, Hesthaven, Spectral?
		R := utils.NewVector(numSamples, rd)
		//fmt.Printf("R = %5.3f\n", R.DataP)
		var plot bool
		plot = false
		if plot { // This is false by default
			chart := utils.NewLineChart(1920, 1080, -1, 1, -4, 4)
			var delay time.Duration
			lineColor := -1.
			nMax := 10
			inc = 2. / float64(nMax-1.)
			var alpha, beta float64
			for n := 0; n < nMax; n++ {
				switch n % 4 {
				case 0:
					//alpha, beta = 0, 0
				case 1:
					alpha, beta = 2*float64(n)+1, 0
				case 2:
					alpha, beta = 0, 2*float64(n)+1
				case 3:
					//alpha, beta = -.5, -.49
				}
				f := DG1D.JacobiP(R, alpha, beta, n)
				if n == nMax-1 {
					delay = 120 * time.Second
				}
				chart.Plot(delay, R.DataP, f, lineColor, "JacobiP["+strconv.Itoa(n)+"]")
				lineColor += inc
			}
		}
	}
}

func LagrangeInterpolationMatrix(R utils.Vector) (I utils.Matrix) {
	var (
		Np = R.Len()
	)
	I = utils.NewMatrix(Np, Np)
	for j := 0; j < Np; j++ {
		Col := LagrangeJthPolynomial(R.DataP, j)
		fmt.Printf("Col = %v\n", Col)
		I.SetCol(j, Col)
	}
	return
}

func LagrangeJthPolynomial(R []float64, J int) (L []float64) {
	/*
		Given a set of k+1 points (R), evaluate the jth lagrange polynomial at point t
	*/
	var (
		Np = len(R)
	)
	L = make([]float64, Np)
	for i, r := range R {
		L[i] = LagrangePolyAtJ(r, R, J)
	}
	return
}

func LagrangePolyAtJ(r float64, R []float64, j int) (f float64) {
	/*
		From https://en.wikipedia.org/wiki/Lagrange_polynomial
		This evaluates the Lagrange polynomial at term J for location R[j]

		The equivalent Newton polynomial is more efficient for repetitive usage

		Given a set of k points in (R), evaluate the jth lagrange polynomial at point t
		Note: j starts at 0
	*/
	var (
		km1 = len(R)
	)
	if j > km1-1 || j < 0 {
		panic("value of j larger than array or less than zero")
	}
	xj := R[j]
	f = 1
	for i, xi := range R {
		if i == j {
			continue
		}
		f *= (r - xi) / (xj - xi)
	}
	return
}

func TestRTElement(t *testing.T) {
	{
		// Check term-wise orthogonal 2D polynomial basis
		N := 2
		R, S := NodesEpsilon(N - 1)
		ii, jj := 1, 1
		p := Simplex2DP(R, S, ii, jj)
		ddr, dds := GradSimplex2DP(R, S, ii, jj)
		Np := R.Len()
		pCheck, ddrCheck, ddsCheck := make([]float64, Np), make([]float64, Np), make([]float64, Np)
		for i, rVal := range R.DataP {
			sVal := S.DataP[i]
			ddrCheck[i], ddsCheck[i] = GradSimplex2DPTerm(rVal, sVal, ii, jj)
			pCheck[i] = Simplex2DPTerm(rVal, sVal, ii, jj)
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
			//var ddr, dds float64
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
	{
		N := 1
		R, S := NodesEpsilon(N - 1)
		rt := NewRTElement(N, R, S)
		fmt.Println(rt.V[0].Print("V0"))
		fmt.Println(rt.V[1].Print("V1"))
		fmt.Println(rt.Div.Print("Div"))
	}
	if true { // Check Divergence for polynomial vector fields of order < N against analytical solution
		Nend := 8
		for N := 1; N < Nend; N++ {
			R, S := NodesEpsilon(N - 1)
			rt := NewRTElement(N, R, S)
			for cOrder := 0; cOrder <= N; cOrder++ {
				fmt.Printf("Check Order = %d, ", cOrder)
				s1, s2, divCheck := checkSolution(rt, cOrder)
				sp := rt.ProjectFunctionOntoBasis(s1, s2)
				sm := utils.NewMatrix(rt.Np, 1, sp)
				divM := rt.Div.Mul(sm)
				//fmt.Println(divM.Print("divM"))
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
		N := 7
		NRT := N + 1
		R, S := NodesEpsilon(N)
		//Nint := R.Len()
		rt := NewRTElement(NRT, R, S)
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
			chart2d.CrossGlyph, chart2d.Solid, utils.GetColor(utils.White)); err != nil {
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
