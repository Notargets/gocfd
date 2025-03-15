package DG2D

import (
	"fmt"
	"math"
)

type EdgeSegmentFluxIndex struct {
	interiorPtsIndex []int     // The interior point supporting each segment
	segmentArea      []float64 // The area for each segment in R,S space
	edgeTris         [][3]int  // used for debugging
}

func (dfr *DFR2D) getEdgeSegmentFluxIndex() (efi *EdgeSegmentFluxIndex) {
	// The output of this function is a list of point indices within the
	// interior of the Solution element that map to segments along the edge
	// boundary of the element.
	//
	// For the P+1 = NpEdge points on each of the three edges for an element,
	// we define P = NpEdge - 1 "interior" segments.
	// This excludes the two line segments adjacent to the corners of the
	// element.
	//
	// Each segment has one supporting point in the interior of the solution
	// element that is used to compute one side of the Riemann problem for
	// each segment. Each segment's area is used to compute an average of the
	// flux over the two segments bordering an edge point.
	// The two "outside" flux edge points adjoining the triangle corners for
	// each of the three edges have only one supporting segment,
	// so they don't use an area average weighting.
	efi = &EdgeSegmentFluxIndex{
		interiorPtsIndex: make([]int, 3*(dfr.FluxElement.NpEdge-1)),
		segmentArea:      make([]float64, 3*(dfr.FluxElement.NpEdge-1)),
	}
	gmRT := dfr.TriangulateRTElement()

	Np_RTBoundary := 3 * (1 + dfr.FluxElement.NpEdge)

	// fmt.Println("TriVerts = ", gmRT.TriVerts)

	// Index of edge segment pairs
	var targetPairs [][2]int
	targetPairs = make([][2]int, Np_RTBoundary-1)
	for i := 0; i < Np_RTBoundary-1; i++ { // Edges of boundary
		targetPairs[i] = [2]int{i, i + 1}
	}

	// Filter out "corner" triangles, those composed of only edge vertices
	var triList [][3]int
	for _, tri := range gmRT.TriVerts {
		edgeBound := int64(Np_RTBoundary)
		if !(tri[0] < edgeBound && tri[1] < edgeBound && tri[2] < edgeBound) {
			// ! Corner tri
			triList = append(triList, [3]int{int(tri[0]), int(tri[1]), int(tri[2])})
		}
	}

	// Build an index from vertex to the list of triangle indices that include that vertex.
	vertexToTris := make(map[int][]int)
	for idx, tri := range triList {
		for _, v := range tri {
			vertexToTris[v] = append(vertexToTris[v], idx)
		}
	}

	// For each target pair, intersect the lists of triangle indices for each vertex.
	triMap := make([][][3]int, len(targetPairs))
	for idx, pair := range targetPairs {
		a, b := pair[0], pair[1]
		trisA := vertexToTris[a]
		trisB := vertexToTris[b]

		// Build a set from one of the lists (choose the smaller list for efficiency).
		set := make(map[int]struct{})
		if len(trisA) < len(trisB) {
			for _, triIdx := range trisA {
				set[triIdx] = struct{}{}
			}
			for _, triIdx := range trisB {
				if _, found := set[triIdx]; found {
					triMap[idx] = append(triMap[idx], triList[triIdx])
				}
			}
		} else {
			for _, triIdx := range trisB {
				set[triIdx] = struct{}{}
			}
			for _, triIdx := range trisA {
				if _, found := set[triIdx]; found {
					triMap[idx] = append(triMap[idx], triList[triIdx])
				}
			}
		}
	}
	convEdgePtIndex := func(GI int) (I int) {
		var (
			NpEdge = dfr.FluxElement.NpEdge
			NpInt  = dfr.FluxElement.NpInt
		)
		switch {
		case GI <= NpEdge:
			I = GI - 1 + 2*NpInt
		case GI > NpEdge && GI <= 2*NpEdge+1:
			I = GI - 2 + 2*NpInt
		case GI > 2*NpEdge+1 && GI <= 3*NpEdge+2:
			I = GI - 3 + 2*NpInt
		}
		return
	}
	lineLength := func(r1, s1, r2, s2 float64) (a float64) {
		a = math.Sqrt((r2-r1)*(r2-r1) + (s2-s1)*(s2-s1))
		return
	}
	// Compose the segment supporting tris
	// var edgeTris [][3]int
	var (
		sk           int
		RFlux, SFlux = dfr.FluxElement.R.DataP, dfr.FluxElement.S.DataP
	)
	for _, tris := range triMap {
		for _, tri := range tris {
			// fmt.Println("Tri: ", tri)
			efi.edgeTris = append(efi.edgeTris, tri)
			var areaPts [2]int
			var areaPtNum int
			for i := 0; i < 3; i++ {
				if tri[i] >= Np_RTBoundary {
					// This is the support point - convert to interior index
					efi.interiorPtsIndex[sk] = tri[i] - Np_RTBoundary
				} else {
					areaPts[areaPtNum] = convEdgePtIndex(tri[i])
					// NpEdge := dfr.FluxElement.NpEdge
					// fmt.Println("NpEdge, edge pts index, converted to RT: ",
					// 	NpEdge,
					// 	tri[i],
					// 	convEdgePtIndex(tri[i]))
					areaPtNum++
				}
			}
			if areaPtNum != 2 {
				fmt.Println("aPtNum, Np_RTBoundary = ", areaPtNum, Np_RTBoundary)
				panic("didn't find a triangle")
			}
			a1, a2 := areaPts[0], areaPts[1]
			r1, s1, r2, s2 := RFlux[a1], SFlux[a1], RFlux[a2], SFlux[a2]
			efi.segmentArea[sk] = lineLength(r1, s1, r2, s2)
			sk++
		}
	}
	return
}
