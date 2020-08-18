package DG2D

import (
	"fmt"
	"image/color"
	"math"
	"testing"

	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"

	"github.com/stretchr/testify/assert"
)

func TestRTElement(t *testing.T) {
	{
		plot := false
		N := 2
		NRT := N + 1
		R, S := NodesEpsilon(N)
		rt := NewRTElement(NRT, R, S)
		assert.NotNil(t, rt.V1)
		//fmt.Println(rt.V1.Print("V1"))
		//fmt.Println(rt.V2.Print("V2"))
		nr, _ := rt.V1.Dims()
		assert.Equal(t, (NRT+1)*(NRT+3), nr)
		assert.Equal(t, (NRT+1)*(NRT+3), rt.R.Len())
		assert.Equal(t, (NRT+1)*(NRT+3), rt.S.Len())

		s1, s2 := make([]float64, rt.R.Len()), make([]float64, rt.R.Len())
		for i := range rt.R.Data() {
			s1[i] = math.Sin(rt.S.Data()[i]*math.Pi) / 5
			s2[i] = math.Sin(rt.R.Data()[i]*math.Pi) / 5
		}
		s1, s2 = rt.ProjectFunctionOntoBasis(s1, s2)
		//s1, s2 = rt.RebuildFunctionFromBasis(s1, s2)
		for i := range rt.R.Data() {
			fmt.Printf("f(%8.3f,%8.3f)= %8.3f,%8.3f\n", rt.R.AtVec(i), rt.S.AtVec(i), s1[i], s2[i])
		}
		var f1, f2 float64
		f1, f2 = rt.Interpolate(-0.5, -0.5, s1, s2)
		fmt.Printf("f1, f2 = %v, %v\n", f1, f2)
		f1, f2 = rt.Interpolate(0.25, 0.25, s1, s2)
		fmt.Printf("f1, f2 = %v, %v\n", f1, f2)

		if plot {
			chart := PlotTestTri(rt.R.Data(), rt.S.Data(), false)
			points := arraysToPoints(rt.R.Data(), rt.S.Data())
			f := arraysToVector(s1, s2)
			_ = chart.AddVectors("test function", points, f, chart2d.Solid, getColor(blue))
			sleepForever()
		}
	}
}

func arraysToVector(r1, r2 []float64) (g [][2]float64) {
	g = make([][2]float64, len(r1))
	for i := range r1 {
		g[i][0] = r1[i]
		g[i][1] = r2[i]
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

func PlotTestTri(x, y []float64, plotGeom bool) (chart *chart2d.Chart2D) {
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
