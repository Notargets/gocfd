package geometry2D

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
)

func TestTriangulate(t *testing.T) {
	{ //Test Legalize Edge
		R := []float64{-1, 1, -1} // Vertices
		S := []float64{-1, -1, 1}
		tm := NewTriMesh(R, S)

		tri := &Tri{}
		tri.AddEdge(NewEdge([2]int{0, 1}, true))
		e2 := NewEdge([2]int{1, 2}, true)
		tri.AddEdge(e2)
		tri.AddEdge(NewEdge([2]int{2, 0}, true))
		tm.AddBoundingTriangle(tri)

		tm.AddPoint(-0.33, -0.33)

		plot := false
		if plot {
			plotTriangles(tm.ToGraphMesh())
			utils.SleepFor(100000)
		}
	}
	if false { //Test Legal Edge test
		R := []float64{-0.9600, 0.9201, -0.9600, -0.7366, 0.4731, -0.7366, -0.3333, -0.0297, -0.9405, -0.0297, 0.7358, -0.9517, -0.7841, -0.7841, -0.9517, 0.7358, 0.4017, -0.9434, -0.4583, -0.4583, -0.9434, 0.4017, 0.0733, -0.7064, -0.3669, -0.3669, -0.7064, 0.0733, -0.9600, 0.9201, -0.9600, -0.7366, 0.4731, -0.7366, -0.3333, -0.0297, -0.9405, -0.0297, 0.7358, -0.9517, -0.7841, -0.7841, -0.9517, 0.7358, 0.4017, -0.9434, -0.4583, -0.4583, -0.9434, 0.4017, 0.0733, -0.7064, -0.3669, -0.3669, -0.7064, 0.0733, 0.9195, 0.7388, 0.4779, 0.1653, -0.1653, -0.4779, -0.7388, -0.9195, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195}
		S := []float64{-0.9600, -0.9600, 0.9201, -0.7366, -0.7366, 0.4731, -0.3333, -0.9405, -0.0297, -0.0297, -0.9517, 0.7358, 0.7358, -0.9517, -0.7841, -0.7841, -0.9434, 0.4017, 0.4017, -0.9434, -0.4583, -0.4583, -0.7064, 0.0733, 0.0733, -0.7064, -0.3669, -0.3669, -0.9600, -0.9600, 0.9201, -0.7366, -0.7366, 0.4731, -0.3333, -0.9405, -0.0297, -0.0297, -0.9517, 0.7358, 0.7358, -0.9517, -0.7841, -0.7841, -0.9434, 0.4017, 0.4017, -0.9434, -0.4583, -0.4583, -0.7064, 0.0733, 0.0733, -0.7064, -0.3669, -0.3669, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000}
		R = append(R, -1, 1, -1) // Vertices
		S = append(S, -1, -1, 1)
		R = append(R, -0.9999999) // Almost at vertex, but inside
		S = append(S, -1)
		testCheck := make([]bool, len(R))
		for i := range testCheck {
			testCheck[i] = true
		}
		testCheck[len(R)-4] = false // vertices will be legal
		testCheck[len(R)-3] = false
		testCheck[len(R)-2] = false
		// Check expected results, most points are inside the circumscribing circle, so will be "illegal"
		for i, r := range R {
			s := S[i]
			fmt.Printf("point[%8.5f,%8.5f] ", r, s)
			assert.Equal(t, IsIllegalEdge(r, s, -1, -1, 1, -1, -1, 1), testCheck[i])
			if testCheck[i] {
				// testCheck is true, so edge is illegal, inside
				fmt.Printf("is inside, illegal (ccw), ")
			} else {
				fmt.Printf("is outside, legal (ccw), ")
			}
			// Test the opposite direction of the base triangle
			assert.Equal(t, IsIllegalEdge(r, s, -1, 1, 1, -1, -1, -1), testCheck[i])
			if testCheck[i] {
				// testCheck is true, so edge is illegal, inside
				fmt.Printf("is inside, illegal\n")
			} else {
				fmt.Printf("is outside, legal\n")
			}
		}
	}

	if false {

		N := 7
		Ninterior := N * (N + 1) / 2
		//R, S := DG2D.NodesEpsilon(N - 1)
		//rt := DG2D.NewRTElement(N, R, S)
		RR := []float64{-0.9600, 0.9201, -0.9600, -0.7366, 0.4731, -0.7366, -0.3333, -0.0297, -0.9405, -0.0297, 0.7358, -0.9517, -0.7841, -0.7841, -0.9517, 0.7358, 0.4017, -0.9434, -0.4583, -0.4583, -0.9434, 0.4017, 0.0733, -0.7064, -0.3669, -0.3669, -0.7064, 0.0733, -0.9600, 0.9201, -0.9600, -0.7366, 0.4731, -0.7366, -0.3333, -0.0297, -0.9405, -0.0297, 0.7358, -0.9517, -0.7841, -0.7841, -0.9517, 0.7358, 0.4017, -0.9434, -0.4583, -0.4583, -0.9434, 0.4017, 0.0733, -0.7064, -0.3669, -0.3669, -0.7064, 0.0733, 0.9195, 0.7388, 0.4779, 0.1653, -0.1653, -0.4779, -0.7388, -0.9195, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195}
		SS := []float64{-0.9600, -0.9600, 0.9201, -0.7366, -0.7366, 0.4731, -0.3333, -0.9405, -0.0297, -0.0297, -0.9517, 0.7358, 0.7358, -0.9517, -0.7841, -0.7841, -0.9434, 0.4017, 0.4017, -0.9434, -0.4583, -0.4583, -0.7064, 0.0733, 0.0733, -0.7064, -0.3669, -0.3669, -0.9600, -0.9600, 0.9201, -0.7366, -0.7366, 0.4731, -0.3333, -0.9405, -0.0297, -0.0297, -0.9517, 0.7358, 0.7358, -0.9517, -0.7841, -0.7841, -0.9434, 0.4017, 0.4017, -0.9434, -0.4583, -0.4583, -0.7064, 0.0733, 0.0733, -0.7064, -0.3669, -0.3669, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000}
		R := []float64{-1, 1, -1}
		S := []float64{-1, -1, 1}
		// Strip off redundant geometry2D
		for i := Ninterior; i < len(RR); i++ {
			R = append(R, RR[i])
			S = append(S, SS[i])
		}

		points := utils.ArraysToPoints(R, S)

		triMesh := graphics2D.TriMesh{
			BaseGeometryClass: graphics2D.BaseGeometryClass{
				Geometry: points,
			},
			Triangles: []graphics2D.Triangle{},
		}

		tri := graphics2D.Triangle{Nodes: [3]int32{0, 1, 2}}
		triMesh.Triangles = append(triMesh.Triangles, tri)
		var (
			chart *chart2d.Chart2D
			delay int
			plot  = false
		)
		if plot {
			delay = 500
			chart = plotTriangles(triMesh)
			utils.SleepFor(delay)
		}
		var ii int
		for i := 3; i < len(points); i++ {
			flipped, tri1, tri2 := LegalizeEdge(i, tri, R, S)
			if flipped {
				fmt.Println("tri was flipped")
				triMesh.Triangles[ii] = tri1
			}
			triMesh.Triangles = append(triMesh.Triangles, tri2)
			ii++
			tri = triMesh.Triangles[ii]
			if plot {
				updateTriMesh(chart, triMesh)
				utils.SleepFor(delay)
			}
			fmt.Printf("Point[%8.5f,%8.5f]\n", points[i].X[0], points[i].X[1])
		}
	}
}

func plotTriangles(triMesh graphics2D.TriMesh) (chart *chart2d.Chart2D) {
	colorMap := utils2.NewColorMap(0, 1, 1)
	box := graphics2D.NewBoundingBox(triMesh.GetGeometry())
	chart = chart2d.NewChart2D(1024, 1024, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
	chart.AddColorMap(colorMap)
	go chart.Plot()

	updateTriMesh(chart, triMesh)
	return
}

func updateTriMesh(chart *chart2d.Chart2D, triMesh graphics2D.TriMesh) {
	if err := chart.AddTriMesh("TriMesh", triMesh.Geometry, triMesh,
		chart2d.CrossGlyph, chart2d.Solid, utils.GetColor(utils.White)); err != nil {
		panic("unable to add graph series")
	}
}
