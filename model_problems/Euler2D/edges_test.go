package Euler2D

import (
	"fmt"
	"image/color"
	"testing"

	"github.com/notargets/gocfd/DG2D"
	"github.com/stretchr/testify/assert"

	"github.com/notargets/avs/utils"
)

func TestEdges2(t *testing.T) {
	// TODO: Test whether Q values in edge storage are correctly stored in
	//  left and right forms (they kinda can't be rn) for each edge
	var (
		P = 2
	)
	dfr := DG2D.NewDFR2D(P, false, "../../DG2D/test_data/test_tris_9.neu")
	// edges := make(EdgeKeySlice, len(dfr.Tris.Edges))
	// es := edges.NewEdgeStorage(dfr)
	_ = dfr
	lines := make(map[color.RGBA][]float32)
	DG2D.AddLine(-1, -1, -1, 1, utils.BLUE, lines)
	DG2D.AddLine(-1, 1, 1, -1, utils.BLUE, lines)
	DG2D.AddLine(1, -1, -1, -1, utils.BLUE, lines)
	fmt.Printf("Order: %d\n", P)
	if testing.Verbose() {
		DG2D.PlotLinesAndText(lines, nil)
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
