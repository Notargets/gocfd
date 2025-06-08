package Euler2D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG2D"
	"github.com/stretchr/testify/assert"
)

func TestEdges2(t *testing.T) {
	// This is mostly a placeholder
	var (
		ip     = ipDefault
		Q_Face [4]utils.Matrix
	)
	meshFile := "../../DG2D/test_data/test_tris_9.neu"
	// edges := make(EdgeKeySlice, len(dfr.Tris.Edges))
	// es := edges.NewEdgeStorage(dfr)
	c := NewEuler(ip, meshFile, 1, false, false)
	NpEdge := c.DFR.FluxElement.NpEdge
	for n := 0; n < 4; n++ {
		Q_Face[n] = utils.NewMatrix(3*NpEdge, c.DFR.K)
		for k := 0; k < c.DFR.K; k++ {
			for i := 0; i < 3*NpEdge; i++ {
				Q_Face[n].Set(i, k, math.Pow(10., float64(k)))
			}
		}
	}
	putEdgeValues(c, Q_Face)
	// lines := make(map[color.RGBA][]float32)
	// DG2D.AddLine(-1, -1, -1, 1, utils2.BLUE, lines)
	// DG2D.AddLine(-1, 1, 1, -1, utils2.BLUE, lines)
	// DG2D.AddLine(1, -1, -1, -1, utils2.BLUE, lines)
	// fmt.Printf("N: %d\n", P)
	// if testing.Verbose() {
	// 	DG2D.PlotLinesAndText(lines, nil)
	// }
}

func putEdgeValues(c *Euler, Q_Face [4]utils.Matrix) {
	var (
		Nedge       = c.DFR.FluxElement.NpEdge
		edgeQValues = make([][4]float64, Nedge)
		pm          = c.Partitions
		edgeKeys    = c.SortedEdgeKeys[0]
	)
	for _, en := range edgeKeys {
		e := c.DFR.Tris.Edges[en]
		var (
			// kLGlobal             = int(e.ConnectedTris[0])
			// kL, KmaxL, myThreadL = pm.GetLocalK(int(e.ConnectedTris[0]))
			kL, KmaxL, _ = pm.GetLocalK(int(e.ConnectedTris[0]))
			edgeNumberL  = int(e.ConnectedTriEdgeNumber[0])
			// normalL              = c.GetFaceNormal(kLGlobal, edgeNumberL)
			shiftL = edgeNumberL * Nedge
		)
		// Store solution for this edge - use left side only per Cockburn and Shu's algorithm for Laplacian, alternate flux sides
		for i := 0; i < Nedge; i++ {
			indL := kL + (i+shiftL)*KmaxL
			for n := 0; n < 4; n++ {
				edgeQValues[i][n] = Q_Face[n].DataP[indL]
			}
		}
		c.EdgeStore.PutEdgeValues(en, EdgeQValues, edgeQValues)
	}
}

func TestEdges(t *testing.T) {
	dfr := DG2D.NewDFR2D(1, false, "../../DG2D/test_data/test_tris_9.neu")
	assert.Equal(t, len(dfr.Tris.Edges), 19)
	edges := make(EdgeKeySlice, len(dfr.Tris.Edges))
	var i int
	for key := range dfr.Tris.Edges {
		edges[i] = key
		i++
	}
	edges.Sort()
	// fmt.Printf("len(Edges) = %d, Edges = %v\n", len(edges), edges)
	l := make([]int, len(edges))
	for i, e := range edges {
		// fmt.Printf("vertex[edge[%d]]=%d\n", i, e.GetVertices(false)[1])
		l[i] = e.GetVertices(false)[1]
	}
	assert.Equal(t, []int{1, 2, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 9}, l)
	edges2 := EdgeKeySliceSortLeft(edges)
	edges2.Sort()
	for i, e := range edges2 {
		// v := e.GetVertices(false)
		// fmt.Printf("vertex2[edge[%d]]=[%d,%d]\n", i, v[0], v[1])
		l[i] = e.GetVertices(false)[0]
	}
	assert.Equal(t, []int{0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 7, 8}, l)
}
