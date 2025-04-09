package DG2D

import (
	"testing"

	"github.com/notargets/gocfd/utils"
)

func TestPlotProjection(t *testing.T) {
	var (
		N = 1
	)
	if !testing.Verbose() {
		return
	}
	// DFR := NewDFR2D(N, false, "test_data/test_10tris_centered.
	// neu")
	meshFile := "test_data/test_10tris_centered.neu"

	dfr := NewDFR2D(N, false, meshFile)
	var Q [4]utils.Matrix
	for n := 0; n < 4; n++ {
		Q[n] = utils.NewMatrix(dfr.SolutionElement.Np, dfr.K)
	}
	// SetTestFieldQ(dfr, RADIAL1TEST, Q)
	// SetTestFieldQ(dfr, NORMALSHOCKTESTM12, Q)
	// SetTestFieldQ(dfr, NORMALSHOCKTESTM18, Q)
	// SetTestFieldQ(dfr, NORMALSHOCKTESTM2, Q)
	SetTestFieldQ(dfr, NORMALSHOCKTESTM5, Q)
	// SetTestFieldQ(dfr, FIXEDVORTEXTEST, Q)
	// SetTestFieldQ(dfr, INTEGERTEST, Q)
	var QInterp [4]utils.Matrix
	for n := 0; n < 4; n++ {
		QInterp[n] = dfr.GraphInterp.Mul(Q[n])
	}
	PlotField(QInterp[0].Transpose().DataP, dfr.GraphMesh, 0.0, 0.0)
}
