package Euler2D

import (
	"fmt"
	"image/color"
	"testing"

	"github.com/notargets/gocfd/DG2D"

	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/InputParameters"
	"github.com/notargets/gocfd/utils"
)

func TestProjection(t *testing.T) {
	var (
		N = 2
	)
	if !testing.Verbose() {
		return
	}
	ip := &InputParameters.InputParameters2D{
		Title:           "",
		CFL:             2.5,
		FluxType:        "Roe",
		InitType:        "Freestream",
		PolynomialOrder: N,
		FinalTime:       0,
		Minf:            2.5,
		Gamma:           1.4,
		Alpha:           0,
		BCs:             nil,
		Limiter:         "PerssonC0",
		Kappa:           3,
		PlotFields:      nil,
	}
	meshFile := "../../DG2D/test_data/test_10tris_centered.neu"
	c := NewEuler(ip, meshFile, 1, false, false)
	// dfr := NewDFR2D(N, false,
	// 	"test_data/test_10tris_centered.neu")
	dfr := c.DFR
	efi := dfr.EdgeSegmentIndex

	fmt.Println("2*NpInt, NpEdge: ", 2*dfr.FluxElement.NpInt, dfr.FluxElement.NpEdge)
	fmt.Println("Segment Indices: ", efi.InteriorPtsIndex)
	fmt.Println("Segment Areas: ", efi.SegmentArea)

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
	for _, tri := range efi.EdgeTris {
		i1, i2, i3 := tri[0], tri[1], tri[2]
		lines[utils2.WHITE] = append(lines[utils2.WHITE],
			gmRT.XY[2*i1], gmRT.XY[2*i1+1],
			gmRT.XY[2*i2], gmRT.XY[2*i2+1],
			gmRT.XY[2*i2], gmRT.XY[2*i2+1],
			gmRT.XY[2*i3], gmRT.XY[2*i3+1],
			gmRT.XY[2*i3], gmRT.XY[2*i3+1],
			gmRT.XY[2*i1], gmRT.XY[2*i1+1])
	}

	var (
		Np      = dfr.SolutionElement.Np
		NpInt   = dfr.FluxElement.NpInt
		NpEdge  = dfr.FluxElement.NpEdge
		text    []DG2D.RenderText
		xyCross []float32
	)
	Q := utils.NewMatrix(Np, 1)
	QEdge := utils.NewMatrix(NpEdge*3, 1)
	for i := 0; i < Np; i++ {
		Q.DataP[i] = float64(i)
		x, y := dfr.SolutionElement.R.DataP[i], dfr.SolutionElement.S.DataP[i]
		xyCross = append(xyCross, float32(x), float32(y))
		text = append(text, DG2D.RenderText{
			Color: utils2.WHITE,
			Text:  fmt.Sprintf("%.2f", Q.DataP[i]),
			Pitch: 36,
			X:     float32(x),
			Y:     float32(y),
		})
	}
	c.ProjectEdges(Q, QEdge, []int{0})
	offset := 2 * NpInt
	for i := 0; i < 3*NpEdge; i++ {
		x, y := dfr.FluxElement.R.DataP[i+offset], dfr.FluxElement.S.DataP[i+offset]
		text = append(text, DG2D.RenderText{
			Color: utils2.WHITE,
			Text:  fmt.Sprintf("%.2f", QEdge.DataP[i]),
			Pitch: 36,
			X:     float32(x),
			Y:     float32(y),
		})
	}
	DG2D.AddCrossHairs(xyCross, utils2.GREEN, lines)
	DG2D.PlotLinesAndText(lines, text)
}
