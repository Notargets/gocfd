package Euler2D

import (
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/InputParameters"
)

func TestPlotEntropyInterpolation(t *testing.T) {
	var (
		N = 7
	)
	if !testing.Verbose() {
		return
	}
	// DFR := NewDFR2D(N, false, "../../DG2D/test_data/test_10tris_centered.
	// neu")
	meshFile := "../../DG2D/test_data/test_10tris_centered.neu"
	ip := &InputParameters.InputParameters2D{
		Title:             "",
		CFL:               2.5,
		FluxType:          "Roe",
		InitType:          "Freestream",
		PolynomialOrder:   N,
		FinalTime:         10,
		Minf:              2,
		Gamma:             1.4,
		Alpha:             0,
		BCs:               nil,
		LocalTimeStepping: true,
		MaxIterations:     10000,
		ImplicitSolver:    false,
		Limiter:           "PerssonC0",
		Kappa:             3,
	}
	c := NewEuler(ip, meshFile, 1, false, false)
	// DG2D.SetTestFieldQ(c.DFR, DG2D.RADIAL1TEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM12, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM18, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM2, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.NORMALSHOCKTESTM5, c.Q[0])
	DG2D.SetTestFieldQ(c.DFR, DG2D.FIXEDVORTEXTEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.DFR, DG2D.INTEGERTEST, c.Q[0])
	// fMGraph := c.DFR.GraphInterp.Mul(c.Q[0][3]).Transpose()
	SwitchToEntropyVariables(c.Q[0], c.FSFar.Gamma)
	var QInterp [4]utils.Matrix
	for n := 0; n < 4; n++ {
		QInterp[n] = c.DFR.GraphInterp.Mul(c.Q[0][n])
	}
	SwitchToConservedVariables(QInterp, c.FSFar.Gamma)
	DG2D.PlotField(QInterp[0].Transpose().DataP, c.DFR.GraphMesh, 0.0, 0.0)
}
