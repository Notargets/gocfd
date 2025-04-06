package DG2D

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestEdgeOptimization(t *testing.T) {
	var (
		NMax = 7
	)
	for N := 2; N <= NMax; N++ {
		dfr := NewDFR2D(N, false)
		epd := dfr.OptimizePointDistribution()
		if testing.Verbose() {
			fmt.Printf("Begin/End Lebesque Constants: %.2f,%.2f\n",
				epd.InitialLebesque, epd.FinalLebesque)
		}
		switch N {
		case 2:
			assert.InDeltaf(t, 1.50, epd.InitialLebesque, 0.01, "")
			assert.InDeltaf(t, 1.42, epd.FinalLebesque, 0.01, "")
		case 3:
			assert.InDeltaf(t, 1.64, epd.InitialLebesque, 0.01, "")
			assert.InDeltaf(t, 1.56, epd.FinalLebesque, 0.01, "")
		case 4:
			assert.InDeltaf(t, 1.78, epd.InitialLebesque, 0.01, "")
			assert.InDeltaf(t, 1.67, epd.FinalLebesque, 0.01, "")
		case 5:
			assert.InDeltaf(t, 4.01, epd.InitialLebesque, 0.01, "")
			assert.InDeltaf(t, 1.78, epd.FinalLebesque, 0.01, "")
		case 6:
			assert.InDeltaf(t, 4.33, epd.InitialLebesque, 0.01, "")
			assert.InDeltaf(t, 1.86, epd.FinalLebesque, 0.01, "")
		case 7:
			assert.InDeltaf(t, 4.48, epd.InitialLebesque, 0.01, "")
			assert.InDeltaf(t, 1.93, epd.FinalLebesque, 0.01, "")
		}
		// fmt.Println("Bottom edge points: R =", epd.RBottom, " S =", epd.SBottom)
		// fmt.Println("Left edge points: R =", epd.RLeft, " S =", epd.SLeft)
		// fmt.Println("Hypotenuse edge points: R =", epd.RHyp, " S =", epd.SHyp)
		// fmt.Printf("Condition numbers:\n")
		// fmt.Printf("  Bottom edge: %.6e\n", epd.CondBottom)
		// fmt.Printf("  Left edge: %.6e\n", epd.CondLeft)
		// fmt.Printf("  Hypotenuse: %.6e\n", epd.CondHyp)
	}
}
