package utils

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
)

func TestTriangulate(t *testing.T) {
	{ //Test Legal Edge test
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

	N := 7
	Ninterior := N * (N + 1) / 2
	//R, S := DG2D.NodesEpsilon(N - 1)
	//rt := DG2D.NewRTElement(N, R, S)
	RR := []float64{-0.9600, 0.9201, -0.9600, -0.7366, 0.4731, -0.7366, -0.3333, -0.0297, -0.9405, -0.0297, 0.7358, -0.9517, -0.7841, -0.7841, -0.9517, 0.7358, 0.4017, -0.9434, -0.4583, -0.4583, -0.9434, 0.4017, 0.0733, -0.7064, -0.3669, -0.3669, -0.7064, 0.0733, -0.9600, 0.9201, -0.9600, -0.7366, 0.4731, -0.7366, -0.3333, -0.0297, -0.9405, -0.0297, 0.7358, -0.9517, -0.7841, -0.7841, -0.9517, 0.7358, 0.4017, -0.9434, -0.4583, -0.4583, -0.9434, 0.4017, 0.0733, -0.7064, -0.3669, -0.3669, -0.7064, 0.0733, 0.9195, 0.7388, 0.4779, 0.1653, -0.1653, -0.4779, -0.7388, -0.9195, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195}
	SS := []float64{-0.9600, -0.9600, 0.9201, -0.7366, -0.7366, 0.4731, -0.3333, -0.9405, -0.0297, -0.0297, -0.9517, 0.7358, 0.7358, -0.9517, -0.7841, -0.7841, -0.9434, 0.4017, 0.4017, -0.9434, -0.4583, -0.4583, -0.7064, 0.0733, 0.0733, -0.7064, -0.3669, -0.3669, -0.9600, -0.9600, 0.9201, -0.7366, -0.7366, 0.4731, -0.3333, -0.9405, -0.0297, -0.0297, -0.9517, 0.7358, 0.7358, -0.9517, -0.7841, -0.7841, -0.9434, 0.4017, 0.4017, -0.9434, -0.4583, -0.4583, -0.7064, 0.0733, 0.0733, -0.7064, -0.3669, -0.3669, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195, -0.9195, -0.7388, -0.4779, -0.1653, 0.1653, 0.4779, 0.7388, 0.9195, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000, -1.0000}
	R := []float64{-1, 1, -1}
	S := []float64{-1, -1, 1}
	// Strip off redundant geometry
	for i := Ninterior; i < len(RR); i++ {
		R = append(R, RR[i])
		S = append(S, SS[i])
	}

	points := ArraysToPoints(R, S)

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
		delay = 400
		chart = plotTriangles(triMesh)
		SleepFor(delay)
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
			SleepFor(delay)
		}
		fmt.Printf("Point[%8.5f,%8.5f]\n", points[i].X[0], points[i].X[1])
	}

	if plot {
		plotTriangles(triMesh)
		SleepFor(delay)
	}
	return
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
		chart2d.CrossGlyph, chart2d.Solid, GetColor(White)); err != nil {
		panic("unable to add graph series")
	}
}

func IsIllegalEdge(prX, prY, piX, piY, pjX, pjY, pkX, pkY float64) bool {
	/*
		pr is a new point for candidate triangle pi-pj-pr
		pi-pj is a shared edge between pi-pj-pk and pi-pj-pr
		if pr lies inside the circle defined by pi-pj-pk:
			- The edge pi-pj should be swapped with pr-pk to make two new triangles:
				pi-pr-pk and pj-pk-pr
	*/
	inCircle := func(ax, ay, bx, by, cx, cy, dx, dy float64) (inside bool) {
		// Calculate handedness, counter-clockwise is (positive) and clockwise is (negative)
		signBit := math.Signbit((bx-ax)*(cy-ay) - (cx-ax)*(by-ay))
		ax_ := ax - dx
		ay_ := ay - dy
		bx_ := bx - dx
		by_ := by - dy
		cx_ := cx - dx
		cy_ := cy - dy
		det := (ax_*ax_+ay_*ay_)*(bx_*cy_-cx_*by_) -
			(bx_*bx_+by_*by_)*(ax_*cy_-cx_*ay_) +
			(cx_*cx_+cy_*cy_)*(ax_*by_-bx_*ay_)
		if signBit {
			return det < 0
		} else {
			return det > 0
		}
	}
	return inCircle(piX, piY, pjX, pjY, pkX, pkY, prX, prY)
}

func LegalizeEdge(index int, tri graphics2D.Triangle, X, Y []float64) (flipped bool, triOut1, triOut2 graphics2D.Triangle) {
	var (
		prX, prY      = X[index], Y[index]
		p1x, p2x, p3x = X[tri.Nodes[0]], X[tri.Nodes[1]], X[tri.Nodes[2]]
		p1y, p2y, p3y = Y[tri.Nodes[0]], Y[tri.Nodes[1]], Y[tri.Nodes[2]]
	)
	if IsIllegalEdge(prX, prY, p1x, p1y, p2x, p2y, p3x, p3y) {
		flipped = true
		triOut1.Nodes[0] = tri.Nodes[0]
		triOut1.Nodes[1] = tri.Nodes[1]
		triOut1.Nodes[2] = int32(index)

		triOut2.Nodes[0] = tri.Nodes[1]
		triOut2.Nodes[1] = tri.Nodes[2]
		triOut2.Nodes[2] = int32(index)
	} else {
		triOut1.Nodes = tri.Nodes
		triOut2.Nodes[0] = tri.Nodes[0]
		triOut2.Nodes[1] = tri.Nodes[2]
		triOut2.Nodes[2] = int32(index)
	}
	return
}
