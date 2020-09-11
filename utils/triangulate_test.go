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

	//N := 7
	//R, S := DG2D.NodesEpsilon(N - 1)
	//rt := DG2D.NewRTElement(N, R, S)
	R := []float64{-0.9600, 0.9201, -0.9600, -0.7366, 0.4731, -0.7366, -0.3333, -0.0297, -0.9405, -0.0297, 0.7358, -0.9517, -0.7841, -0.7841, -0.9517, 0.7358, 0.4017, -0.9434, -0.4583, -0.4583, -0.9434, 0.4017, 0.0733, -0.7064, -0.3669, -0.3669, -0.7064, 0.0733, -0.9600, 0.9201, -0.9600, -0.7366, 0.4731, -0.7366, -0.3333, -0.0297, -0.9405, -0.0297, 0.7358, -0.9517, -0.7841, -0.7841, -0.9517, 0.7358, 0.4017, -0.9434, -0.4583, -0.4583, -0.9434, 0.4017, 0.0733, -0.7064, -0.3669, -0.3669, -0.7064, 0.0733, 0.9195, 0.7388, 0.4779, 0.1653, -0.1653, -0.4779, -0.7388, -0.9195, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195}
	S := []float64{-0.9600, -0.9600, 0.9201, -0.7366, -0.7366, 0.4731, -0.3333, -0.9405, -0.0297, -0.0297, -0.9517, 0.7358, 0.7358, -0.9517, -0.7841, -0.7841, -0.9434, 0.4017, 0.4017, -0.9434, -0.4583, -0.4583, -0.7064, 0.0733, 0.0733, -0.7064, -0.3669, -0.3669, -0.9600, -0.9600, 0.9201, -0.7366, -0.7366, 0.4731, -0.3333, -0.9405, -0.0297, -0.0297, -0.9517, 0.7358, 0.7358, -0.9517, -0.7841, -0.7841, -0.9434, 0.4017, 0.4017, -0.9434, -0.4583, -0.4583, -0.7064, 0.0733, 0.0733, -0.7064, -0.3669, -0.3669, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000}

	points = make([]graphics2D.Point, 3+len(R))
	points[0].X[0], points[0].X[1] = -1, -1
	points[1].X[0], points[1].X[1] = 1, -1
	points[2].X[0], points[2].X[1] = -1, 1
	for i := 0; i < len(R); i++ {
		points[i+3].X[0], points[i+3].X[1] = float32(R[i]), float32(S[i])
	}

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
