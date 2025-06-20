package DG2D

import (
	"fmt"
	"image/color"
	"testing"

	"github.com/notargets/avs/utils"
)

func _TestNodes2DUniform(t *testing.T) {
	var (
		P = 2
		N = P * 6
	)
	rs := UniformRS(N)
	pts := GetInteriorPoints(rs)
	lines := make(map[color.RGBA][]float32)
	AddLine(-1, -1, -1, 1, utils.BLUE, lines)
	AddLine(-1, 1, 1, -1, utils.BLUE, lines)
	AddLine(1, -1, -1, -1, utils.BLUE, lines)
	for _, pt := range pts {
		AddCrossHairs([]float32{float32(pt.R), float32(pt.S)}, utils.WHITE, lines)
	}
	fmt.Printf("N: %d\n", N)
	if testing.Verbose() {
		PlotLinesAndText(lines, nil)
	}
}

func _TestNodes2DOptimizedBaseline(t *testing.T) {
	var (
		P     = 7
		Alpha = 0.7
	)
	cond, leb := AnalyzeRTOutput(P, append(MakeEdgePointsFromDist(
		gaussLegendreNodes(P+1)), UniformRSAlpha(P-1, Alpha)...))
	fmt.Printf("Uniform Points with GL Edges:\n P=%d: cond(V)=%.3e, "+
		"Edge Lebesgue≈%.3f\n", P, cond, leb)
	R, S := MakeRSFromPoints(WilliamsShunnJameson(P - 1))
	cond, leb = AnalyzeRTOutput(P,
		AddEdgePointsToInterior(gaussLegendreNodes(P+1), R, S))
	fmt.Printf("Epsilon Points with GL Edges:\n P=%d: cond(V)=%.3e, "+
		"Edge Lebesgue≈%.3f\n", P, cond, leb)
	cond, leb = AnalyzeRTOutput(P,
		AddEdgePointsToInterior(GetOptimizedEdgePoints(P), R, S))
	fmt.Printf("Epsilon Points with Opt Edges:\n P=%d: cond(V)=%.3e, "+
		"Edge Lebesgue≈%.3f\n", P, cond, leb)
	// pts := append(MakeEdgePointsFromDist(gaussLegendreNodes(P+1)),
	// 	UniformRSAlpha(P-1, Alpha)...)
	// pts := AddEdgePointsToInterior(GetOptimizedEdgePoints(P), R, S)
	// pts := AddEdgePointsToInterior(P, gaussLegendreNodes(P+1), R, S)
	R, S = MakeRSFromPoints(WilliamsShunnJameson(P - 1))
	// pts := AddEdgePointsToInterior(gaussLegendreNodes(P+1), R, S)
	// pts := AddEdgePointsToInterior(GetOptimizedEdgePoints(P), R, S)
	pts := AddEdgePointsToInterior(GetOptimizedEdgePoints(P), R, S)
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
