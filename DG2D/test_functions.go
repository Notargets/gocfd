package DG2D

import (
	"math"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/utils"
)

type ScalarTestField interface {
	P(r, s float64, P int) (val float64)
	Gradient(r, s float64, P int) (grad [2]float64)
}

type VectorTestField interface {
	F(r, s float64, P int) (f1, f2 float64)
	Divergence(r, s float64, P int) (div float64)
	String() string
}

type SinCosVectorField struct{}

func (scf SinCosVectorField) String() string { return "Sin / Cos Field" }

func (scf SinCosVectorField) F(r, s float64, P int) (f1, f2 float64) {
	var (
		Pi = math.Pi
	)
	f1, f2 = math.Sin(Pi*(r+1)), math.Cos(Pi*(s+1))
	return
}

func (scf SinCosVectorField) Divergence(r, s float64, P int) (div float64) {
	var (
		Pi = math.Pi
	)
	div = Pi * (math.Cos(Pi*(r+1)) - math.Sin(Pi*(s+1)))
	return
}

type PolyVectorField struct{}

func (lpf PolyVectorField) String() string { return "R + S Permutation Field" }

func (lpf PolyVectorField) F(r, s float64, P int) (f1, f2 float64) {
	var (
		p = float64(P)
	)
	f1, f2 = math.Pow(r+s+10, p), math.Pow(10*(r+s), p)
	return
}

func (lpf PolyVectorField) Divergence(r, s float64, P int) (div float64) {
	var (
		p = float64(P)
	)
	if P > 0 {
		div = p * (math.Pow(r+s+10, p-1) + 10*math.Pow(10*(r+s), p-1))
	}
	return
}

type PolyVectorField3 struct{}

func (lpf PolyVectorField3) String() string { return "[s+10,10r]^p Zero Divergence Field" }

func (lpf PolyVectorField3) F(r, s float64, P int) (f1, f2 float64) {
	var (
		p = float64(P)
	)
	f1, f2 = math.Pow(s+10, p), math.Pow(10*r, p)
	return
}

func (lpf PolyVectorField3) Divergence(r, s float64, P int) (div float64) {
	return
}

type PolyVectorField2 struct{}

func (lpf PolyVectorField2) String() string { return "[r,s]^p Simple Field" }

func (lpf PolyVectorField2) F(r, s float64, P int) (f1, f2 float64) {
	var (
		p = float64(P)
	)
	f1, f2 = math.Pow(r, p), math.Pow(s, p)
	return
}

func (lpf PolyVectorField2) Divergence(r, s float64, P int) (div float64) {
	var (
		p = float64(P)
	)
	if P > 0 {
		div = p * (math.Pow(r, p-1) + math.Pow(s, p-1))
	}
	return
}

type PolyScalarField struct{}

func (lpf PolyScalarField) P(r, s float64, P int) (val float64) {
	var (
		p = float64(P)
	)
	val = math.Pow(10*r+s+10, p)
	return
}

func (lpf PolyScalarField) Gradient(r, s float64, P int) (grad [2]float64) {
	var (
		p = float64(P)
	)
	if P > 0 {
		grad = [2]float64{
			10. * p * math.Pow(10*r+s+10, p-1.),
			p * math.Pow(10*r+s+10, p-1.),
		}
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
		if err := chart.AddTriMesh("TriMesh", trimesh,
			chart2d.CrossGlyph, 0.1, chart2d.Solid,
			utils.GetColor(utils.Black)); err != nil {
			panic("unable to add graph series")
		}
	}
	return
}
