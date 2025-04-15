package Euler2D

import (
	"testing"

	"github.com/notargets/gocfd/DG2D"

	"github.com/notargets/gocfd/InputParameters"

	"github.com/notargets/gocfd/utils"
)

func TestPlotProjection(t *testing.T) {
	var (
		N = 4
	)
	if !testing.Verbose() {
		return
	}
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
	dfr := c.DFR
	// gp.MassMatrix.Print("MM")
	// gp.Minv.Print("Minv")
	var Q, Q_Face, Q_Face_P0 [4]utils.Matrix
	for n := 0; n < 4; n++ {
		Q[n] = utils.NewMatrix(dfr.SolutionElement.Np, dfr.K)
		Q_Face[n] = utils.NewMatrix(dfr.FluxElement.NpEdge*3, dfr.K)
		Q_Face_P0[n] = utils.NewMatrix(dfr.FluxElement.NpEdge*3, dfr.K)
	}
	// SetTestFieldQ(dfr, RADIAL1TEST, Q)
	// SetTestFieldQ(dfr, NORMALSHOCKTESTM12, Q)
	// SetTestFieldQ(dfr, NORMALSHOCKTESTM18, Q)
	// SetTestFieldQ(dfr, NORMALSHOCKTESTM2, Q)
	DG2D.SetTestFieldQ(dfr, DG2D.NORMALSHOCKTESTM5, Q)
	// SetTestFieldQ(dfr, FIXEDVORTEXTEST, Q)
	// SetTestFieldQ(dfr, INTEGERTEST, Q)

	NpEdge := dfr.FluxElement.NpEdge
	var QInterp [4]utils.Matrix
	for n := 0; n < 4; n++ {
		QInterp[n] = dfr.GraphInterp.Mul(Q[n])
	}
	n := 0
	c.ShockFinder.UpdateShockedCells(Q[0])
	c.InterpolateSolutionToShockedEdges(c.ShockFinder, Q_Face, Q_Face_P0)
	for _, k := range c.ShockFinder.ShockCells.Cells() {
		// QInterp has 1 extra points per edge
		for nEdge := 0; nEdge < 3; nEdge++ {
			off1 := nEdge * NpEdge
			off2 := nEdge * (NpEdge + 1)
			for i := 0; i < NpEdge; i++ {
				QInterp[n].Set(off2+i+1, k, Q_Face[n].At(off1+i, k))
			}
		}
	}
	// Replace vertices for plotting
	dfr.AverageGraphFieldVertices(QInterp[n])
	// }
	DG2D.PlotField(QInterp[0].Transpose().DataP, dfr.GraphMesh, 0.0, 0.0)
}
