package DG2D

import (
	"fmt"
	"image/color"
	"testing"

	utils2 "github.com/notargets/avs/utils"
)

func TestPlotEdgeTriangulation(t *testing.T) {
	var (
		N = 2
	)
	if !testing.Verbose() {
		return
	}
	dfr := NewDFR2D(N, false,
		"test_data/test_10tris_centered.neu")
	efi := dfr.getEdgeSegmentFluxIndex()

	fmt.Println("2*NpInt, NpEdge: ", 2*dfr.FluxElement.NpInt, dfr.FluxElement.NpEdge)
	fmt.Println("Segment Indices: ", efi.interiorPtsIndex)
	fmt.Println("Segment Areas: ", efi.segmentArea)

	gmRT := dfr.TriangulateRTElement()
	Np_RTBoundary := 3 * (1 + dfr.FluxElement.NpEdge)
	lines := make(map[color.RGBA][]float32)
	for i := 0; i < Np_RTBoundary; i++ {
		ip := i + 1
		if i == Np_RTBoundary-1 {
			ip = 0
		}
		lines[utils2.RED] = append(lines[utils2.RED],
			gmRT.XY[2*i], gmRT.XY[2*i+1],
			gmRT.XY[2*ip], gmRT.XY[2*ip+1])
	}
	for _, tri := range efi.edgeTris {
		i1, i2, i3 := tri[0], tri[1], tri[2]
		lines[utils2.WHITE] = append(lines[utils2.WHITE],
			gmRT.XY[2*i1], gmRT.XY[2*i1+1],
			gmRT.XY[2*i2], gmRT.XY[2*i2+1],
			gmRT.XY[2*i2], gmRT.XY[2*i2+1],
			gmRT.XY[2*i3], gmRT.XY[2*i3+1],
			gmRT.XY[2*i3], gmRT.XY[2*i3+1],
			gmRT.XY[2*i1], gmRT.XY[2*i1+1])
	}
	var xyCross []float32
	for i := 0; i < dfr.SolutionElement.Np; i++ {
		x, y := dfr.SolutionElement.R.DataP[i], dfr.SolutionElement.S.DataP[i]
		xyCross = append(xyCross, float32(x), float32(y))
	}
	addCrossHairs(xyCross, utils2.GREEN, lines)
	PlotLines(lines)
}
