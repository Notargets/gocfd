package DG2D

import (
	"testing"
)

func TestRomeroJamesonRTBasis_ComposePolyBasis(t *testing.T) {
	for P := 1; P <= 6; P++ {
		Np := (P + 1) * (P + 3)
		var JJ int
		for i := 0; i <= P+3; i++ { // Allowing up to P+1
			for j := 0; j <= (P + 3 - i); j++ {
				JJ++
			}
		}
		t.Logf("P = %d, Np = %d, JJ = %d\n", P, Np, JJ)
	}

}
