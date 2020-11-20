package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"

	"github.com/stretchr/testify/assert"
)

func TestRTElement(t *testing.T) {
	{
		// Check term-wise orthogonal 2D polynomial basis
		N := 2
		R, S := NodesEpsilon(N - 1)
		a, b := RStoAB(R, S)
		ii, jj := 1, 1
		p := Simplex2DP(a, b, ii, jj)
		ddr, dds := GradSimplex2DP(a, b, ii, jj)
		Np := R.Len()
		pCheck, ddrCheck, ddsCheck := make([]float64, Np), make([]float64, Np), make([]float64, Np)
		for i, rVal := range R.Data() {
			sVal := S.Data()[i]
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
			r := rt.R.Data()[i]
			s := rt.S.Data()[i]
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
				minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, divM.Data(), divCheck)
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
		//RR := rt.R.Subset(Nint, rt.R.Len()-1)
		//SS := rt.S.Subset(Nint, rt.S.Len()-1)
		//fmt.Println(RR.Transpose().Print("RR"))
		//fmt.Println(SS.Transpose().Print("SS"))
		s1, s2 := make([]float64, rt.R.Len()), make([]float64, rt.R.Len())
		for i := range rt.R.Data() {
			s1[i] = 1
			s2[i] = 1
		}
		if plot {
			chart := PlotTestTri(true)
			points := utils.ArraysToPoints(rt.R.Data(), rt.S.Data())
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
			r := rt.Rm.Data()[i]
			s := rt.Sm.Data()[i]
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
				div := Dr.Mul(s1m).Add(Ds.Mul(s2m)).Data()
				minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, div, divCheck)
				assert.True(t, near(minerrInt, 0.0, 0.00001))
				assert.True(t, near(maxerrInt, 0.0, 0.00001))
				assert.True(t, near(minerrEdge, 0.0, 0.00001))
				assert.True(t, near(maxerrEdge, 0.0, 0.00001))
			}
		}
	}
*/
