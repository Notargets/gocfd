package Euler2D

import (
	"fmt"
	"testing"

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
	// DG2D.SetTestFieldQ(c.dfr, DG2D.RADIAL1TEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM12, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM18, c.Q[0])
	DG2D.SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM2, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM5, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.FIXEDVORTEXTEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.INTEGERTEST, c.Q[0])

	// c.Q[0][0].String("Density")
	// c.Q[0][1].String("Umom")
	// c.Q[0][2].String("Vmom")
	// c.Q[0][3].String("E")
	// os.Exit(1)

	// Dens := c.Q[0][0]
	// Kmax := c.dfr.K
	// scratch := c.ShockFinder.Qalt.DataP
	// for k := 0; k < Kmax; k++ {
	// 	for i := 0; i < c.dfr.SolutionElement.Np; i++ {
	// 		ind := k + i*Kmax
	// 		scratch[i] = Dens.DataP[ind]
	// 	}
	// 	sigma := c.ShockFinder.ShockIndicator(scratch)
	// 	fmt.Printf("ShockFinder Sigma[%d] = %f\n", k, sigma)
	// }
	fieldNum := 3
	d1, d2 := c.Q[0][fieldNum].Dims()
	fmt.Printf("Q dims before interpolation = %d:%d\n", d1, d2)
	// fMGraph := c.dfr.GraphInterp.Mul(c.Q[0][3]).Transpose()
	QInterp := c.dfr.GraphInterp.Mul(c.Q[0][fieldNum])
	d1, d2 = QInterp.Dims()
	fmt.Printf("Q dims after interpolation = %d:%d\n", d1, d2)
	// For testing, zero out the edge and vertex values
	// for k := 0; k < Kmax; k++ {
	// 	for i := 0; i < 3*(c.dfr.FluxElement.NpEdge+1); i++ {
	// 		QInterp.Set(i, k, 0.)
	// 	}
	// }
	c.GetFirstOrderEdgeProjection_ForGraphing(c.Q[0][fieldNum], &QInterp)
	// QInterp.String("QInterp")
	DG2D.PlotField(QInterp.Transpose().DataP, c.dfr.GraphMesh, 0.0, 0.0)
}
