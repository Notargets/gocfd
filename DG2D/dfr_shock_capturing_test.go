package DG2D

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/utils"
	"github.com/stretchr/testify/assert"
)

func TestShockFinder(t *testing.T) {
	var (
		NMin, NMax = 2, 4
	)
	meshFile := "test_data/test_10tris_centered.neu"
	for N := NMin; N <= NMax; N++ {
		t.Logf("Order: %d", N)
		dfr := NewDFR2D(N, false, meshFile)
		var (
			Kmax  = dfr.K
			Np    = dfr.SolutionElement.Np
			Kappa = 3.
		)
		// c := NewEuler(ip, meshFile, 1, false, false)
		// SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM12, c.Q[0])
		var Q [4]utils.Matrix
		for n := 0; n < 4; n++ {
			Q[n] = utils.NewMatrix(Np, Kmax)
		}
		SetTestFieldQ(dfr, NORMALSHOCKTESTM2, Q)
		F, T := false, true
		shockedElements := []bool{F, F, T, T, F, F, F, T, T, F}
		// SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM5, c.Q[0])
		// SetTestFieldQ(c.dfr, DG2D.FIXEDVORTEXTEST, c.Q[0])
		Dens := Q[0]
		// scratch := c.ShockFinder.Qalt.DataP
		ShockFinder := dfr.NewAliasShockFinder(Kappa)
		scratch := make([]float64, Np)
		for k := 0; k < Kmax; k++ {
			for i := 0; i < dfr.SolutionElement.Np; i++ {
				ind := k + i*Kmax
				scratch[i] = Dens.DataP[ind]
			}
			sigma := ShockFinder.ShockIndicator(scratch)
			fmt.Printf("ShockFinder Sigma[%d] = %f\n", k, sigma)
			hasShock := ShockFinder.ElementHasShock(scratch)
			assert.Equal(t, hasShock, shockedElements[k])
		}
	}
}
