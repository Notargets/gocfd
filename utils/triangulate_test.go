package utils

import (
	"testing"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
)

func TestTriangulate(t *testing.T) {
	var (
		points  []graphics2D.Point
		triMesh graphics2D.TriMesh
		K       = 1
	)

	points = make([]graphics2D.Point, 3)
	points[0].X[0], points[0].X[1] = -1, -1
	points[1].X[0], points[1].X[1] = 1, -1
	points[2].X[0], points[2].X[1] = -1, 1

	triMesh.Triangles = make([]graphics2D.Triangle, K)
	triMesh.Geometry = points
	triMesh.Triangles[0].Nodes[0] = 0
	triMesh.Triangles[0].Nodes[1] = 1
	triMesh.Triangles[0].Nodes[2] = 2

	plot := false
	if plot {
		plotTriangles(triMesh)
	}
	return
}

func plotTriangles(triMesh graphics2D.TriMesh) {
	var (
		chart *chart2d.Chart2D
	)
	colorMap := utils2.NewColorMap(0, 1, 1)
	box := graphics2D.NewBoundingBox(triMesh.GetGeometry())
	chart = chart2d.NewChart2D(1024, 1024, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
	chart.AddColorMap(colorMap)
	go chart.Plot()

	if err := chart.AddTriMesh("TriMesh", triMesh.Geometry, triMesh,
		chart2d.CrossGlyph, chart2d.Solid, GetColor(White)); err != nil {
		panic("unable to add graph series")
	}
	SleepForever()
}
