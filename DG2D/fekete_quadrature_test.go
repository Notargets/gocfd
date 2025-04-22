package DG2D

import (
	"fmt"
	"image/color"
	"testing"

	"github.com/notargets/avs/utils"
)

func TestNodes2DFekete(t *testing.T) {
	var (
		N = 2
	)
	r, s, _ := Nodes2DFekete(N)
	var lines map[color.RGBA][]float32
	lines = make(map[color.RGBA][]float32)
	var interiorCount, extCount int
	var vertexCount int
	for i, rr := range r {
		ss := s[i]
		AddCrossHairs([]float32{float32(rr), float32(ss)}, utils.WHITE, lines)
		if ss == -1 || rr == -1 ||
			ss == 1 || rr == 1 ||
			rr+ss == 0 {
			fmt.Printf("Ext Point: [%.2f,%.2f]\n", rr, ss)
			extCount++
		}
		if (ss == -1 && rr == -1) ||
			(ss == 1 && rr == -1) ||
			(ss == -1 && rr == 1) {
			vertexCount++
		}
	}
	totalCount := len(r)
	interiorCount = totalCount - extCount
	singleEdge := (totalCount - interiorCount - vertexCount) / 3
	fmt.Printf("Order: %d\n", N)
	fmt.Printf("Exterior count: %d\n", extCount)
	fmt.Printf("Vertex count: %d\n", vertexCount)
	fmt.Printf("Interior count: %d\n", interiorCount)
	fmt.Printf("Single Edge count: %d\n", singleEdge)
	fmt.Printf("Total count: %d\n", totalCount)
	if testing.Verbose() {
		PlotLinesAndText(lines, nil)
	}
}
