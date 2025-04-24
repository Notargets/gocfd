package DG2D

import (
	"fmt"
	"math"
	"testing"
)

// truncateMantissa truncates x so that only the top p bits of its 52-bit mantissa are retained.
func round(x float64, p uint) float64 {
	I := 1.
	for i := 0; i < int(p); i++ {
		I *= 10
	}
	return math.Round(x*I) / I
}

func TestRound(t *testing.T) {
	f := -0.3849001817874452
	ft := round(f, 7)
	fmt.Printf("f, ft := %15.10f, %15.10f\n", f, ft)
}

func TestEdgeOptimization(t *testing.T) {
	// This only needs to be done once and the output placed into the edge
	// distributions in raviart_thomas_element.go
	if false {
		var (
			NMin = 1
			NMax = 1
			tol  = 1.e-5
		)
		for N := NMin; N <= NMax; N++ {
			fmt.Printf("RT Order %d\n", N+1)
			R, S := NodesEpsilon(N)
			SolutionBasis := NewJacobiBasis2D(N, R, S, 0, 0)
			epd := OptimizePointDistribution(N, SolutionBasis)
			//		fmt.Println(epd.RBottom)
			if testing.Verbose() {
				fmt.Printf("Begin/End Lebesque Constants: %.2f,%.2f\n",
					epd.InitialLebesque, epd.FinalLebesque)
				for i, r := range epd.RBottom {
					if math.Abs(r) < tol {
						epd.RBottom[i] = 0.
					} else {
						epd.RBottom[i] = round(r, 8)
					}
				}
				fmt.Printf("R = %#v\n", epd.RBottom)
			}
		}
	}
}
