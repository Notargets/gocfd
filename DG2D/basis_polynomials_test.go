package DG2D

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG1D"
)

func TestNewLagrangePolynomial1D(t *testing.T) {
	P := 2
	R := utils.NewVector(P+1, DG1D.LegendreZeros(P))

	lp1d := NewLagrangePolynomial1D(R, P, 0, 0)
	// Check the Lagrange property
	for j := 0; j < lp1d.Np; j++ {
		for i := 0; i < lp1d.Np; i++ {
			r := R.AtVec(i)
			psi := lp1d.getPolynomial(r, j)
			fmt.Printf("psi[%d,%d,%f]=%f\n", i, j, r, psi)
			if i == j {
				assert.InDeltaf(t, psi, 1, 0.000001, "")
			} else {
				assert.InDeltaf(t, psi, 0, 0.000001, "")
			}
		}
	}
}
