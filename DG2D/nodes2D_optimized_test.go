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

func TestNodes2DOptimizedBaseline(t *testing.T) {
	var (
		P     = 4
		Alpha = 0.7
	)
	cond, leb := AnalyzeRTOutput(P, append(MakeEdgePointsFromDist(
		gaussLegendreNodes(P+1)), UniformRSAlpha(P-1, Alpha)...))
	fmt.Printf("Uniform Points with GL Edges:\n P=%d: cond(V)=%.3e, "+
		"Edge Lebesgue≈%.3f\n", P, cond, leb)
	R, S := NodesEpsilon(P - 1)
	cond, leb = AnalyzeRTOutput(P,
		AddEdgePointsToInterior(gaussLegendreNodes(P+1), R, S))
	fmt.Printf("Epsilon Points with GL Edges:\n P=%d: cond(V)=%.3e, "+
		"Edge Lebesgue≈%.3f\n", P, cond, leb)
	cond, leb = AnalyzeRTOutput(P,
		AddEdgePointsToInterior(GetOptimizedEdgePointsEpsilon(P), R, S))
	fmt.Printf("Epsilon Points with Opt Edges:\n P=%d: cond(V)=%.3e, "+
		"Edge Lebesgue≈%.3f\n", P, cond, leb)
	// pts := append(MakeEdgePointsFromDist(gaussLegendreNodes(P+1)),
	// 	UniformRSAlpha(P-1, Alpha)...)
	// pts := AddEdgePointsToInterior(GetOptimizedEdgePointsEpsilon(P), R, S)
	// pts := AddEdgePointsToInterior(P, gaussLegendreNodes(P+1), R, S)
	ptsInt := WilliamsShunnJameson(P - 1)
	R, S = MakeRSFromPoints(ptsInt)
	// pts := AddEdgePointsToInterior(gaussLegendreNodes(P+1), R, S)
	// pts := AddEdgePointsToInterior(GetOptimizedEdgePointsEpsilon(P), R, S)
	pts := AddEdgePointsToInterior(GetOptimizedEdgePointsEpsilon(P), R, S)
	cond, leb = AnalyzeRTOutput(P, pts)
	fmt.Printf("Williams, Shunn, Jameson Points GL Edges:\n P=%d: cond("+
		"V)=%.3e, "+
		"Edge Lebesgue≈%.3f\n", P, cond, leb)

	lines := make(map[color.RGBA][]float32)
	AddLine(-1, -1, -1, 1, utils.BLUE, lines)
	AddLine(-1, 1, 1, -1, utils.BLUE, lines)
	AddLine(1, -1, -1, -1, utils.BLUE, lines)
	for _, pt := range pts {
		AddCrossHairs([]float32{float32(pt.R), float32(pt.S)}, utils.WHITE, lines)
	}
	if testing.Verbose() {
		PlotLinesAndText(lines, nil)
	}
}
