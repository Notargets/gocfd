package Euler2D

import (
	"fmt"
	"os"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/InputParameters"
)

func TestPlotVariousFields2(t *testing.T) {
	var (
		N = 7
	)
	if !testing.Verbose() {
		return
	}
	// dfr := NewDFR2D(N, false, "../../DG2D/test_data/test_10tris_centered.
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
	// SetTestFieldQ(c.dfr, DG2D.RADIAL1TEST, c.Q[0])
	// SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM12, c.Q[0])
	SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM2, c.Q[0])
	// SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM5, c.Q[0])
	// SetTestFieldQ(c.dfr, DG2D.FIXEDVORTEXTEST, c.Q[0])
	// SetTestFieldQ(c.dfr, DG2D.INTEGERTEST, c.Q[0])

	// c.Q[0][0].Print("Density")
	// c.Q[0][1].Print("Umom")
	// c.Q[0][2].Print("Vmom")
	// c.Q[0][3].Print("E")
	// os.Exit(1)

	Dens := c.Q[0][0]
	Kmax := c.dfr.K
	scratch := c.ShockFinder.Qalt.DataP
	for k := 0; k < Kmax; k++ {
		for i := 0; i < c.dfr.SolutionElement.Np; i++ {
			ind := k + i*Kmax
			scratch[i] = Dens.DataP[ind]
		}
		sigma := c.ShockFinder.ShockIndicator(scratch)
		fmt.Printf("ShockFinder Sigma[%d] = %f\n", k, sigma)
	}
	os.Exit(1)
	fMGraph := c.dfr.GraphInterp.Mul(c.Q[0][3]).Transpose()
	DG2D.PlotField(fMGraph.DataP, c.dfr.GraphMesh, 0.0, 0.0)
}

func SetTestFieldQ(dfr *DG2D.DFR2D, tf DG2D.TestField, Q [4]utils.Matrix) {
	var (
		Np, Kmax = Q[0].Dims()
	)
	X, Y := dfr.SolutionX.Transpose().DataP, dfr.SolutionY.Transpose().DataP
	field := DG2D.SetTestField(X, Y, tf)
	for n := 0; n < 4; n++ {
		for k := 0; k < Kmax; k++ {
			for i := 0; i < Np; i++ {
				Q[n].Set(i, k, field[i+k*Np+n*Np*Kmax])
			}
		}
	}
	return
}
