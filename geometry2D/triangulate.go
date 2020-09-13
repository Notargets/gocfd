package geometry2D

import (
	"math"

	graphics2D "github.com/notargets/avs/geometry"
)

type Point struct {
	X [2]float64
}

type Edge struct {
	Tris []*Tri // Associated triangles
}

type Tri struct {
	Verts [3]int
	Edges [3]int
}

type TriMesh struct {
	Tris   []*Tri
	Edges  []*Edge
	Points []Point
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

func NewTriMesh(X, Y []float64) (tris *TriMesh) {
	pts := make([]Point, len(X))
	for i, x := range X {
		y := Y[i]
		pts[i].X[0] = x
		pts[i].X[1] = y
	}
	tris = &TriMesh{
		Tris:   nil,
		Edges:  nil,
		Points: pts,
	}
	return
}

func (tm *TriMesh) LegalizeEdge(e *Edge, testPtI int) {
	// The edge must have two triangles associated with it for the test to make sense
	if len(e.Tris) < 2 {
		return
	}
	// Find which of the two triangles does not contain the test point
	triNotContainingPt := func(ptI int, tris []*Tri) (tri *Tri) {
		for _, tri = range tris {
			var contains bool
			for _, vert := range tri.Verts {
				if vert == ptI {
					contains = true
					break
				}
			}
			if !contains {
				return
			}
		}
		return
	}
	baseTri := triNotContainingPt(testPtI, e.Tris)
	_ = baseTri
	return
}

func (tm *TriMesh) ToGraphMesh() (trisOut graphics2D.TriMesh) {
	pts := make([]graphics2D.Point, len(tm.Points))
	for i, pt := range tm.Points {
		pts[i].X[0] = float32(pt.X[0])
		pts[i].X[1] = float32(pt.X[1])
	}
	tris := make([]graphics2D.Triangle, len(tm.Tris))
	for i, tri := range tm.Tris {
		tris[i].Nodes[0] = int32(tri.Verts[0])
		tris[i].Nodes[1] = int32(tri.Verts[1])
		tris[i].Nodes[2] = int32(tri.Verts[2])
	}
	trisOut = graphics2D.TriMesh{
		BaseGeometryClass: graphics2D.BaseGeometryClass{
			Geometry: pts,
		},
		Triangles:  tris,
		Attributes: nil,
	}
	return
}
