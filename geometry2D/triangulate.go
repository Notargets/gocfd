package geometry2D

import (
	"github.com/notargets/avs/geometry"
	"github.com/pradeep-pyro/triangle"
)

func ConvertEdgePlusInternal(EdgeX, EdgeY,
	X, Y []float64) (triXY []float32, constrainedEdge [][2]int32) {
	var (
		EdgeLen = len(EdgeX)
	)
	triXY = make([]float32, 2*(EdgeLen+len(X)))
	for i := 0; i < EdgeLen; i++ {
		triXY[2*i+0] = float32(EdgeX[i])
		triXY[2*i+1] = float32(EdgeY[i])
	}
	for i := EdgeLen; i < len(X)+EdgeLen; i++ {
		triXY[2*i+0] = float32(X[i-EdgeLen])
		triXY[2*i+1] = float32(Y[i-EdgeLen])
	}
	constrainedEdge = make([][2]int32, EdgeLen)
	for i := 0; i < EdgeLen; i++ {
		ii := int32(i)
		ip := ii + 1
		if i == EdgeLen-1 {
			ip = 0
		}
		constrainedEdge[i] = [2]int32{ii, ip}
	}
	return
}

func TriangulateTriangle(EdgeX, EdgeY, X, Y []float64) (tm geometry.TriMesh) {
	var (
		EdgeLen = len(EdgeX)
	)
	// The input EdgeX and EdgeY should be a contiguous array of points
	// defining the edges of the triangle
	triXY, constrainedEdge := ConvertEdgePlusInternal(EdgeX, EdgeY, X, Y)
	input := triangle.NewTriangulateIO()

	points := make([][2]float64, len(X)+EdgeLen)
	for i := 0; i < len(triXY)/2; i++ {
		points[i][0] = float64(triXY[2*i+0])
		points[i][1] = float64(triXY[2*i+1])
	}
	input.SetPoints(points)

	// Define constrained edges (segments)
	input.SetSegments(constrainedEdge)

	// Configure triangulation options
	// options := "pzQ" // p: PSLG, z: zero indexing, Q: quiet
	options := triangle.NewOptions()
	options.ConformingDelaunay = false // Disable conforming Delaunay (
	// which adds points)
	options.SegmentSplitting = triangle.NoSplitting // Prevent segment splitting
	options.MaxSteinerPoints = 0                    // Explicitly prohibit any Steiner points
	// options.Area = 0.0                              // Disable area constraint
	options.Angle = 0.0 // Disable angle constraint

	// Run triangulation
	// output := triangle.Triangulate(input, options, true)
	output := triangle.Triangulate(input, options, false)

	tris := output.Triangles()
	tris64 := make([][3]int64, len(tris))
	for i, tri := range tris {
		for n := 0; n < 3; n++ {
			tris64[i][n] = int64(tri[n])
		}
	}
	tm = geometry.NewTriMesh(triXY, tris64)

	return
}
