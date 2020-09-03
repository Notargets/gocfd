package DG2D

import (
	"fmt"
	"image/color"
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"

	"github.com/stretchr/testify/assert"
)

func TestRTElement(t *testing.T) {
	{ // RT0 Validation
		oosr2 := 1. / math.Sqrt(2)
		R, S := NodesEpsilon(0)
		rt := NewRTElement(0, R, S)
		// Vandermonde Matrices, one for each of r,s directions hand calculated for the RT0 case
		// Note: The Vandermonde matrix is defined as: V_i_j = Psi_j(X_i), j is column number
		checkV1 := utils.NewMatrix(3, 3, []float64{
			oosr2, -.5, .5,
			0, -1, 0,
			oosr2, -.5, .5,
		})
		checkV2 := utils.NewMatrix(3, 3, []float64{
			oosr2, .5, -.5,
			oosr2, .5, -.5,
			0, 0, -1,
		})
		assert.True(t, nearVec(checkV1.Data(), rt.V1.Data(), 0.000001))
		assert.True(t, nearVec(checkV2.Data(), rt.V2.Data(), 0.000001))
		for j := range rt.R.Data() {
			r, s := rt.R.AtVec(j), rt.S.AtVec(j)
			p1, p2 := rt.EvaluatePolynomial(j, r, s)
			//			fmt.Printf("poly[%d] at (%8.5f, %8.5f) = [%8.5f, %8.5f]\n", j, r, s, p1, p2)
			switch j {
			case 0: // Edge 1
				assert.True(t, near(oosr2, p1, 0.00001))
				assert.True(t, near(oosr2, p2, 0.00001))
			case 1: // Edge 2
				assert.True(t, near(-1, p1, 0.00001))
				assert.True(t, near(0.5, p2, 0.00001))
			case 2: // Edge 3
				assert.True(t, near(0.5, p1, 0.00001))
				assert.True(t, near(-1, p2, 0.00001))
			}
		}
	}
	{ // RT 1 Validation
		N := 0
		NRT := N + 1
		R, S := NodesEpsilon(N)
		rt := NewRTElement(NRT, R, S)
		//assert.Equal(t, (NRT+1)*(NRT+3), rt.R.Len())
		//assert.Equal(t, (NRT+1)*(NRT+3), rt.S.Len())
		assert.Equal(t, rt.Npm, rt.R.Len())
		assert.Equal(t, rt.Npm, rt.S.Len())
		/*
			Reconstruct the single coefficient matrix to compare with the Matlab solution
		*/
		Np := (NRT + 1) * (NRT + 3)
		Ainv := utils.NewMatrix(Np, Np)
		NGroup1 := (NRT+1)*NRT/2 + (NRT + 1) // The first set of polynomial terms
		for i := 0; i < NGroup1; i++ {
			rowG1 := rt.A1.Row(i)
			Ainv.M.SetRow(i, rowG1.Data())
			rowG2 := rt.A2.Row(i + NGroup1)
			Ainv.M.SetRow(i+NGroup1, rowG2.Data())
		}
		for i := 2 * NGroup1; i < Np; i++ {
			rowG3 := rt.A1.Row(i)
			Ainv.M.SetRow(i, rowG3.Data())
		}
		// Matlab solution
		CheckAinv := utils.NewMatrix(Np, Np, []float64{
			0.75, -0.75, 0.3536, 0.3536, 0.4045, -0.1545, -0.4045, 0.1545,
			-0.75, -1.5, 1.279, 0.4886, 0.25, 0.25, -0.9635, 0.7135,
			-0.75, -1.5, 0.135, 0.9256, 0.559, -0.559, -0.6545, -0.09549,
			-0.75, 0.75, 0.3536, 0.3536, -0.4045, 0.1545, 0.4045, -0.1545,
			-1.5, -0.75, 0.9256, 0.135, -0.6545, -0.09549, 0.559, -0.559,
			-1.5, -0.75, 0.4886, 1.279, -0.9635, 0.7135, 0.25, 0.25,
			-1.5, -0.75, 0.9256, 0.135, -0.6545, -0.09549, -0.559, 0.559,
			-0.75, -1.5, 0.135, 0.9256, -0.559, 0.559, -0.6545, -0.09549,
		})
		assert.True(t, nearVec(CheckAinv.Data(), Ainv.Data(), 0.001))

		// Verify Vandermonde matrices against Matlab solution
		CheckV1 := utils.NewMatrix(Np, Np, []float64{
			1.0, 0, 0, 0, 0, 0, 0, 0,
			1.0, 0, 0, 0, 0, 0, 0, 0,
			0.6, -0.6, 1.023, 0, 0.2472, 0.07639, -0.5236, 0.6472,
			0.6, -0.6, 0, 0.3909, 0.5236, -0.6472, -0.2472, -0.07639,
			0, 0, 0, 0, -1.0, 0, 0, 0,
			0, 0, 0, 0, 0, -1.0, 0, 0,
			1.2, 0.6, -0.108, -0.3496, -0.6472, 0.5236, 0.2764, 0,
			1.2, 0.6, 0.9153, -0.7405, 0.07639, 0.2472, 0, 0.7236,
		})
		CheckV2 := utils.NewMatrix(Np, Np, []float64{
			0, 1.0, 0, 0, 0, 0, 0, 0,
			0, 1.0, 0, 0, 0, 0, 0, 0,
			-0.6, 0.6, 0.3909, 0, -0.2472, -0.07639, 0.5236, -0.6472,
			-0.6, 0.6, 0, 1.023, -0.5236, 0.6472, 0.2472, 0.07639,
			0.6, 1.2, -0.3496, -0.108, 0.2764, 0, -0.6472, 0.5236,
			0.6, 1.2, -0.7405, 0.9153, 0, 0.7236, 0.07639, 0.2472,
			0, 0, 0, 0, 0, 0, -1.0, 0,
			0, 0, 0, 0, 0, 0, 0, -1.0,
		})
		if false {
			CheckV1, CheckV2 = CombineBasis(rt.N, CheckV1, CheckV2)
		}
		assert.True(t, nearVec(CheckV1.Data(), rt.V1.Data(), 0.001))
		assert.True(t, nearVec(CheckV2.Data(), rt.V2.Data(), 0.001))
		// Validate derivative matrices against Matlab solution
		CheckDr1 := utils.NewMatrix(Np, Np, []float64{
			0.5, -0.5, 0.6171, 0.09003, 0.8727, 0.1273, -0.3727, 0.3727,
			0.5, -0.5, 0.6171, 0.09003, 0.8727, 0.1273, -0.3727, 0.3727,
			-1.756, -1.5, 2.047, 0.1954, -0.08541, -0.08541, -1.171, 1.256,
			0.2562, -1.5, 0.5117, 0.7818, 0.5854, 0.5854, -0.7562, 0.1708,
			2.585, 0.6708, -0.6325, -0.1954, 1.809, 0.191, 0.4472, -0.3618,
			1.915, -0.6708, -0.5117, 0.6325, 1.309, 0.691, -0.1382, -0.4472,
			1.342, 0.6708, 0.3162, -0.5578, 1.394, -0.2236, 0.191, 0.309,
			-1.342, -0.6708, 1.972, -0.3162, 0.2236, -0.3944, -0.809, 1.309,
		})
		CheckDr2 := utils.NewMatrix(Np, Np, []float64{
			-1.0, -0.5, 0.6171, 0.09003, -0.4363, -0.06366, 0.7454, -0.7454,
			-1.0, -0.5, 0.6171, 0.09003, -0.4363, -0.06366, 0.7454, -0.7454,
			-0.8292, -0.4146, 0.5117, 0.07465, -0.3618, -0.05279, 0.809, -0.809,
			-2.171, -1.085, 1.34, 0.1954, -0.9472, -0.1382, 0.309, -0.309,
			-0.8292, -0.4146, 0.5117, 0.07465, -0.3618, -0.05279, 0.809, -0.809,
			-2.171, -1.085, 1.34, 0.1954, -0.9472, -0.1382, 0.309, -0.309,
			0, 0, 0, 0, 0, 0, 1.118, -1.118,
			0, 0, 0, 0, 0, 0, 1.118, -1.118,
		})
		if false {
			CheckDr1, CheckDr2 = CombineBasis(rt.N, CheckDr1, CheckDr2)
		}
		assert.True(t, nearVec(CheckDr1.Data(), rt.Dr1.Data(), 0.001))
		assert.True(t, nearVec(CheckDr2.Data(), rt.Dr2.Data(), 0.001))
		CheckDs1 := utils.NewMatrix(Np, Np, []float64{
			-0.5, -1.0, 0.09003, 0.6171, 0.7454, -0.7454, -0.4363, -0.06366,
			-0.5, -1.0, 0.09003, 0.6171, 0.7454, -0.7454, -0.4363, -0.06366,
			-1.085, -2.171, 0.1954, 1.34, 0.309, -0.309, -0.9472, -0.1382,
			-0.4146, -0.8292, 0.07465, 0.5117, 0.809, -0.809, -0.3618, -0.05279,
			0, 0, 0, 0, 1.118, -1.118, 0, 0,
			0, 0, 0, 0, 1.118, -1.118, 0, 0,
			-0.4146, -0.8292, 0.07465, 0.5117, 0.809, -0.809, -0.3618, -0.05279,
			-1.085, -2.171, 0.1954, 1.34, 0.309, -0.309, -0.9472, -0.1382,
		})
		CheckDs2 := utils.NewMatrix(Np, Np, []float64{
			-0.5, 0.5, 0.09003, 0.6171, -0.3727, 0.3727, 0.8727, 0.1273,
			-0.5, 0.5, 0.09003, 0.6171, -0.3727, 0.3727, 0.8727, 0.1273,
			-1.5, 0.2562, 0.7818, 0.5117, -0.7562, 0.1708, 0.5854, 0.5854,
			-1.5, -1.756, 0.1954, 2.047, -1.171, 1.256, -0.08541, -0.08541,
			0.6708, 1.342, -0.5578, 0.3162, 0.191, 0.309, 1.394, -0.2236,
			-0.6708, -1.342, -0.3162, 1.972, -0.809, 1.309, 0.2236, -0.3944,
			0.6708, 2.585, -0.1954, -0.6325, 0.4472, -0.3618, 1.809, 0.191,
			-0.6708, 1.915, 0.6325, -0.5117, -0.1382, -0.4472, 1.309, 0.691,
		})
		if false {
			CheckDs1, CheckDs2 = CombineBasis(rt.N, CheckDs1, CheckDs2)
		}
		assert.True(t, nearVec(CheckDs1.Data(), rt.Ds1.Data(), 0.001))
		assert.True(t, nearVec(CheckDs2.Data(), rt.Ds2.Data(), 0.001))
	}
	plot := false
	if plot {
		N := 6
		NRT := N + 1
		R, S := NodesEpsilon(N)
		rt := NewRTElement(NRT, R, S)
		s1, s2 := make([]float64, rt.R.Len()), make([]float64, rt.R.Len())
		for i := range rt.R.Data() {
			/*
				s1[i] = math.Sin(rt.S.Data()[i]*math.Pi) / 5
				s2[i] = math.Sin(rt.R.Data()[i]*math.Pi) / 5
			*/
			s1[i] = 1
			s2[i] = 1
		}
		s1, s2 = rt.ProjectFunctionOntoBasis2(s1, s2)

		if plot {
			chart := PlotTestTri(true)
			points := arraysToPoints(rt.R.Data(), rt.S.Data())
			f := arraysToVector(s1, s2, 0.1)
			_ = chart.AddVectors("test function", points, f, chart2d.Solid, getColor(green))
			sleepForever()
		}
	}

	// Check Divergence
	{
		N := 1
		R, S := NodesEpsilon(N - 1)
		rt := NewRTElement(N, R, S)
		divCheck := make([]float64, rt.Npm)
		F1, F2 := make([]float64, rt.Npm), make([]float64, rt.Npm)
		fmt.Println(rt.A1.Print("A1"))
		fmt.Println(rt.A2.Print("A2"))
		for i := range rt.R.Data() {
			F1[i], F2[i] = 1., 1.
			r := rt.R.Data()[i]
			s := rt.S.Data()[i]
			coeffs1 := rt.A1.Col(i).Data()
			coeffs2 := rt.A2.Col(i).Data()
			a2 := coeffs1[1]
			a7, a8 := coeffs1[6], coeffs1[7]
			a5 := coeffs2[4]
			// Manually calculated divergence for RT1 element
			divCheck[i] = F1[i]*(a2+2*a7*r+a8*s) + F2[i]*(a5+a7*r+2*a8*s)
			fmt.Printf("i, a2,a5,a7,a8 = %d, %8.5f, %8.5f, %8.5f, %8.5f\n", i, a2, a5, a7, a8)
		}
		div := rt.Divergence(F1, F2)
		fmt.Printf("div= %v\n", div)
		fmt.Printf("divCheck = %v\n", divCheck)
		// TODO: Fix divergence calculation to match the check
		//assert.True(t, nearVec(div, divCheck, 0.00001))
	}

	if false {
		Nend := 8
		for N := 1; N < Nend; N++ {
			R, S := NodesEpsilon(N - 1)
			rt := NewRTElement(N, R, S)
			Npm := (N + 6) * (N + 1) / 2
			s1, s2 := make([]float64, Npm), make([]float64, Npm)
			for i := 0; i < Npm; i++ {
				r := rt.R.Data()[i]
				s := rt.S.Data()[i]
				s1[i] = r * r
				s2[i] = s * s
			}
			s1, s2 = rt.ProjectFunctionOntoBasis2(s1, s2)
			div := rt.Divergence(s1, s2)
			// Restrict divergence to internal points
			errors := make([]float64, Npm)
			for i := 0; i < Npm; i++ {
				r := rt.R.Data()[i]
				s := rt.S.Data()[i]
				// d/dR
				ddr := 2 * r
				// d/dS
				dds := 2 * s
				divCheck := ddr + dds
				error := div[i] - divCheck
				errors[i] = error
			}
			minerrInt, maxerrInt := errors[0], errors[0]
			Nint := N * (N + 1) / 2
			minerrEdge, maxerrEdge := errors[Nint], errors[Nint]
			for i := 0; i < Nint; i++ {
				errAbs := math.Abs(errors[i])
				if minerrInt > errAbs {
					minerrInt = errAbs
				}
				if maxerrInt < errAbs {
					maxerrInt = errAbs
				}
			}
			for i := Nint; i < Npm; i++ {
				errAbs := math.Abs(errors[i])
				if minerrEdge > errAbs {
					minerrEdge = errAbs
				}
				if maxerrEdge < errAbs {
					maxerrEdge = errAbs
				}
			}
			fmt.Printf("Order = %d, ", N)
			fmt.Printf("Min, Max Int Err = %8.5f, %8.5f, Min, Max Edge Err = %8.5f, %8.5f\n", minerrInt, maxerrInt, minerrEdge, maxerrEdge)
			//fmt.Printf("Errors = %v\n", errors)
		}
	}
}

func arraysToVector(r1, r2 []float64, scaleO ...float64) (g [][2]float64) {
	var (
		scale float64 = 1
	)
	g = make([][2]float64, len(r1))
	if len(scaleO) > 0 {
		scale = scaleO[0]
	}
	for i := range r1 {
		g[i][0] = r1[i] * scale
		g[i][1] = r2[i] * scale
	}
	return
}

func arraysToPoints(r1, r2 []float64) (points []graphics2D.Point) {
	points = make([]graphics2D.Point, len(r1))
	for i := range r1 {
		points[i].X[0] = float32(r1[i])
		points[i].X[1] = float32(r2[i])
	}
	return
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
		if err := chart.AddTriMesh("TriMesh", points, trimesh,
			chart2d.CrossGlyph, chart2d.Solid, getColor(white)); err != nil {
			panic("unable to add graph series")
		}
	}
	return
}

type ColorName uint8

const (
	white ColorName = iota
	blue
	red
	green
	black
)

func getColor(name ColorName) (c color.RGBA) {
	switch name {
	case white:
		c = color.RGBA{
			R: 255,
			G: 255,
			B: 255,
			A: 0,
		}
	case blue:
		c = color.RGBA{
			R: 50,
			G: 0,
			B: 255,
			A: 0,
		}
	case red:
		c = color.RGBA{
			R: 255,
			G: 0,
			B: 50,
			A: 0,
		}
	case green:
		c = color.RGBA{
			R: 25,
			G: 255,
			B: 25,
			A: 0,
		}
	case black:
		c = color.RGBA{
			R: 0,
			G: 0,
			B: 0,
			A: 0,
		}
	}
	return
}
