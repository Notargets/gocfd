package DG2D

import (
	"fmt"
	"image/color"
	"testing"

	"github.com/notargets/avs/utils"
)

func TestNodes2DUniform(t *testing.T) {
	var (
		P = 2
		N = P * 6
	)
	rs := UniformRS(N)
	pts := getInteriorPoints(rs)
	lines := make(map[color.RGBA][]float32)
	AddLine(-1, -1, -1, 1, utils.BLUE, lines)
	AddLine(-1, 1, 1, -1, utils.BLUE, lines)
	AddLine(1, -1, -1, -1, utils.BLUE, lines)
	for _, pt := range pts {
		AddCrossHairs([]float32{float32(pt.R), float32(pt.S)}, utils.WHITE, lines)
	}
	fmt.Printf("Order: %d\n", N)
	if testing.Verbose() {
		PlotLinesAndText(lines, nil)
	}
}

func TestNodes2DOptimized(t *testing.T) {
	var (
		P = 3
		N = P * 6
	)
	pts, err := OptimizeRTNodesCont(P)
	if err != nil {
		t.Fatal(err)
	}
	lines := make(map[color.RGBA][]float32)
	AddLine(-1, -1, -1, 1, utils.BLUE, lines)
	AddLine(-1, 1, 1, -1, utils.BLUE, lines)
	AddLine(1, -1, -1, -1, utils.BLUE, lines)
	for _, pt := range pts {
		AddCrossHairs([]float32{float32(pt.R), float32(pt.S)}, utils.WHITE, lines)
	}
	fmt.Printf("Order: %d\n", N)
	if testing.Verbose() {
		PlotLinesAndText(lines, nil)
	}
}
