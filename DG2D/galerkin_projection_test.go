package DG2D

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/utils"
)

func TestPlotProjection(t *testing.T) {
	var (
		N  = 7
		NP = 0
	)
	if !testing.Verbose() {
		return
	}
	// DFR := NewDFR2D(N, false, "test_data/test_10tris_centered.
	// neu")
	meshFile := "test_data/test_10tris_centered.neu"

	dfr := NewDFR2D(N, false, meshFile)
	Rn, Sn := MakeRSFromPoints(WilliamsShunnJameson(NP))
	nb := NewJacobiBasis2D(NP, Rn, Sn)
	gp := NewGalerkinProjection(dfr.SolutionBasis, nb)
	_ = gp
	// gp.MassMatrix.Print("MM")
	// gp.Minv.Print("Minv")
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
		// Replace vertices for plotting
		dfr.AverageGraphFieldVertices(QInterp[n])
	}
	ShockFinder := dfr.NewAliasShockFinder(3)
	// Check for shocked elements, project those and interpolate new edge values
	NpInt := dfr.SolutionElement.Np
	Uh := utils.NewMatrix(NpInt, 1)
	GraphR, GraphS := dfr.GetRSForGraphMesh()
	NpEdge := dfr.FluxElement.NpEdge
	NpEdgeGraph := 3 * (NpEdge + 1)
	grE := GraphR.Subset(0, NpEdgeGraph-1)
	gsE := GraphS.Subset(0, NpEdgeGraph-1)
	IProj := gp.GetProjectedInterpolationMatrix(grE, gsE)
	for k := 0; k < dfr.K; k++ {
		for i := 0; i < NpInt; i++ {
			Uh.DataP[i] = Q[0].At(i, k)
		}
		if ShockFinder.ElementHasShock(Uh.DataP) {
			fmt.Println("Shock found at K = ", k)
			vals := IProj.Mul(Uh)
			for i := 0; i < NpEdgeGraph; i++ {
				QInterp[0].Set(i, k, vals.DataP[i])
			}
		}
	}
	dfr.AverageGraphFieldVertices(QInterp[0])
	PlotField(QInterp[0].Transpose().DataP, dfr.GraphMesh, 0.0, 0.0)
}
