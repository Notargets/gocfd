package DG2D

import (
	"fmt"
	"image/color"
	"log"
	"sort"
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
		P     = 3
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
	pts := append(MakeEdgePointsFromDist(gaussLegendreNodes(P+1)),
		UniformRSAlpha(P-1, Alpha)...)
	// pts := AddEdgePointsToInterior(GetOptimizedEdgePointsEpsilon(P), R, S)
	// pts := AddEdgePointsToInterior(P, gaussLegendreNodes(P+1), R, S)
	R, S = MakeRSFromPoints(pts)
	R.Transpose().Print("R")
	S.Transpose().Print("S")

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

func TestNodes2DOptimized3(t *testing.T) {
	var (
		P = 2
	)
	res, err := TuneRT(P)
	if err != nil {
		t.Fatal(err)
	}
	SortInPlaceByLebEdge(res)
	for _, r := range res {
		r.Print()
	}
}

func SortInPlaceByLebEdge(results []Result) {
	sort.Slice(results, func(i, j int) bool {
		return results[i].LebEdge < results[j].LebEdge
	})
}

type Result struct {
	Eps, Alpha, Gamma float64
	Cond, LebEdge     float64
}

func (r Result) Print() {
	fmt.Printf("Eps=%.2f, Alpha=%.0e, Gamma=%.1f, Cond=%.3e, EdgeLeb=%.3f\n",
		r.Eps, r.Alpha, r.Gamma, r.Cond, r.LebEdge)
}

func TuneRT(P int) ([]Result, error) {
	epsilons := []float64{0.05, 0.10, 0.15}  // λ_min floor
	alphas := []float64{1e2, 1e3, 1e4}       // hinge strengths
	gammas := []float64{0.0, 0.1, 1.0, 10.0} // cond weights

	var out []Result
	for _, eps := range epsilons {
		for _, alpha := range alphas {
			for _, gamma := range gammas {
				pts, err := OptimizeRTWithParams(P, eps, alpha, gamma)
				if err != nil {
					return nil, err
				}
				c, l := AnalyzeRTOutput(P, pts)
				out = append(out, Result{Eps: eps, Alpha: alpha, Gamma: gamma, Cond: c, LebEdge: l})
				// fmt.Printf("eps=%.2f α=%.0e γ=%.1f → cond=%.3e, edge‐Leb=%.3f\n",
				// 	eps, alpha, gamma, c, l)
			}
		}
	}
	return out, nil
}
func TestNodes2DOptimized2(t *testing.T) {
	var (
		P = 2
	)
	// 1) build the 3×(P+1)=9 Gauss–Legendre edge points exactly
	gl := gaussLegendreNodes(P + 1) // length=3
	edges := make([]Point, 0, 3*len(gl))
	for _, t := range gl {
		edges = append(edges,
			Point{R: t, S: -1}, // bottom
			Point{R: -1, S: t}, // left
			Point{R: t, S: -t}, // hypotenuse
		)
	}

	// 2) define the “baseline” RT‐P=2 interior points in barycentric form
	bary := [][3]float64{
		{0.5, 0.25, 0.25},
		{0.25, 0.5, 0.25},
		{0.25, 0.25, 0.5},
	}
	interior := make([]Point, len(bary))
	for i, b := range bary {
		r := 2*b[1] - 1
		s := 2*b[2] - 1
		interior[i] = Point{R: r, S: s}
	}

	// 3) baseline full set: edges then interior
	baseline := append(edges, interior...)

	// 4) your optimized result from before:
	optimized, err := OptimizeRTWithParams(P, 1.e-1, 1e3, 1.0)
	if err != nil {
		log.Fatal(err)
	}

	// 5) analyze both
	cond0, leb0 := AnalyzeRTOutput(P, baseline)
	cond1, leb1 := AnalyzeRTOutput(P, optimized)

	fmt.Printf("P=%d baseline: cond≈%.3f, edge‐Leb≈%.3f\n", P, cond0, leb0)
	fmt.Printf("P=%d optimized: cond≈%.3f, edge‐Leb≈%.3f\n", P, cond1, leb1)
}
