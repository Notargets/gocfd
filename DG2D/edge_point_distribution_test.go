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
	if false {
		var (
			NMax = 7
			tol  = 1.e-5
		)
		for N := 0; N <= NMax; N++ {
			fmt.Printf("Order %d\n", N)
			dfr := NewDFR2D(N, false)
			epd := dfr.OptimizePointDistribution()
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
