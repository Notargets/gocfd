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
			fmt.Printf("poly[%d] at (%8.5f, %8.5f) = [%8.5f, %8.5f]\n", j, r, s, p1, p2)
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
		nr, _ := rt.V1.Dims()
		assert.Equal(t, (NRT+1)*(NRT+3), nr)
		assert.Equal(t, (NRT+1)*(NRT+3), rt.R.Len())
		assert.Equal(t, (NRT+1)*(NRT+3), rt.S.Len())
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
		fmt.Println(Ainv.Print("Ainv"))
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
		s1, s2 = rt.ProjectFunctionOntoBasis(s1, s2)

		if plot {
			chart := PlotTestTri(true)
			points := arraysToPoints(rt.R.Data(), rt.S.Data())
			f := arraysToVector(s1, s2, 0.1)
			_ = chart.AddVectors("test function", points, f, chart2d.Solid, getColor(green))
			sleepForever()
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
