package Euler2D

import (
	"fmt"
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
	// DG2D.SetTestFieldQ(c.dfr, DG2D.RADIAL1TEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM12, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM2, c.Q[0])
	DG2D.SetTestFieldQ(c.dfr, DG2D.NORMALSHOCKTESTM5, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.FIXEDVORTEXTEST, c.Q[0])
	// DG2D.SetTestFieldQ(c.dfr, DG2D.INTEGERTEST, c.Q[0])

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
	fieldNum := 3
	d1, d2 := c.Q[0][fieldNum].Dims()
	fmt.Printf("Q dims before interpolation = %d:%d\n", d1, d2)
	// fMGraph := c.dfr.GraphInterp.Mul(c.Q[0][3]).Transpose()
	QInterp := c.dfr.GraphInterp.Mul(c.Q[0][fieldNum])
	d1, d2 = QInterp.Dims()
	fmt.Printf("Q dims after interpolation = %d:%d\n", d1, d2)
	// For testing, zero out the edge and vertex values
	for k := 0; k < Kmax; k++ {
		for i := 0; i < 3*(c.dfr.FluxElement.NpEdge+1); i++ {
			QInterp.Set(i, k, 0.)
		}
	}
	c.FirstOrderEdgeProjection_ForGraphing(c.Q[0][fieldNum], &QInterp)
	// QInterp.Print("QInterp")
	DG2D.PlotField(QInterp.Transpose().DataP, c.dfr.GraphMesh, 0.0, 0.0)
}

func (c *Euler) FirstOrderEdgeProjection_ForGraphing(Q utils.Matrix,
	QGraph *utils.Matrix) {
	var (
		dfr         = c.dfr
		NpInt, KMax = Q.Dims()
		NpEdge      = dfr.FluxElement.NpEdge
		NpGraph     = 3*(1+NpEdge) + NpInt
		efi         = dfr.GetEdgeSegmentFluxIndex()
	)
	if QGraph.IsEmpty() {
		QGraphP := utils.NewMatrix(NpGraph, KMax)
		*QGraph = QGraphP
	}
	for k := 0; k < KMax; k++ {
		// There are NpEdge-1 interior points supporting reconstruction of
		// NpEdge-1 sub-segments on each of the three edges
		// Here we will use the two adjoining corner segments to construct
		// the vertex value and we'll average segments to create interior
		// node values

		// Below in reconstructed efi coordinates:
		// Nseg = NpE-1,  Nefi = 3Nseg
		Nseg := NpEdge - 1
		Nefi := 3 * Nseg
		if len(efi.InteriorPtsIndex) != Nefi {
			panic("edge segment point count is not correct")
		}
		// We'll do the vertices first
		// NpE == NpEdge, below in GraphMesh coordinates
		// v1,   Edge1,   v2,        Edge2,         v3         Edge3
		//  0,  1->NpE,  NpE+1, (NpE+2)->2*NpE+1, 2*NpE+2, 2*NpE+3->3*NpE+2

		// Efi coordinates:
		//      v1              v2                 v3
		// 0+(Nefi-1)/2  ((Nseg-1)+Nseg)/2  (2Nseg-1)+2Nseg/2
		fmt.Println("NpEdge, Nseg, Nefi, NpInt = ", NpEdge, Nseg, Nefi, NpInt)
		fmt.Println("InteriorPtsIndex = ", efi.InteriorPtsIndex)
		v1a, v1b := efi.InteriorPtsIndex[0], efi.InteriorPtsIndex[Nefi-1]
		fmt.Println("v1a/b = ", v1a, v1b)
		QGraph.Set(0, k, 0.5*(Q.At(v1a, k)+Q.At(v1b, k)))
		v2a, v2b := efi.InteriorPtsIndex[Nseg-1], efi.InteriorPtsIndex[Nseg]
		fmt.Println("v2a/b = ", v2a, v2b)
		QGraph.Set(NpEdge+1, k, 0.5*(Q.At(v2a, k)+Q.At(v2b, k)))
		v3a, v3b := efi.InteriorPtsIndex[2*Nseg-1], efi.InteriorPtsIndex[2*Nseg]
		fmt.Println("v3a/b = ", v3a, v3b)
		QGraph.Set(2*NpEdge+2, k, 0.5*(Q.At(v3a, k)+Q.At(v3b, k)))

		// Beginning/End Edge Points (0-based), inclusive:
		// Edge1			Edge2				Edge3
		// b:0 e:Nseg-1		b:Nseg e:2*Nseg-1	b:2*Nseg e:3*Nseg-1 or e:Nefi-1

		// Below are loop indices in efi coordinates (non-inclusive)
		// Edge1         Edge2              Edge3
		// 0->Nseg       Nseg->2Nseg        2Nseg->Nefi

		// NpE == NpEdge, below in GraphMesh coordinates (ranges non-inclusive)
		// v1,   Edge1,   v2,        Edge2,         v3         Edge3
		//  0,  1->NpE+1, NpE+1, (NpE+2)->2*NpE+2, 2*NpE+2, 2*NpE+3->3*NpE+3
		if true {
			var skEdge, skSeg int
			for n := 0; n < 3; n++ { // Each edge
				// Beginning point of range, excluding vertex
				// QGraph.Set(n*NpEdge+n+1, k, Q.At(efi.InteriorPtsIndex[skSeg], k))
				skEdge = 1 + n*(NpEdge+1) // Skip the vertex
				sleft := efi.InteriorPtsIndex[skSeg]
				QGraph.Set(skEdge, k, Q.At(sleft, k))
				skEdge++
				// Interior range
				for i := 1; i < NpEdge-1; i++ { // Averaging segment values
					sright := efi.InteriorPtsIndex[skSeg+1]
					// QGraph.Set(i+n*NpEdge+n+1, k, 0.5*(Q.At(sleft, k)+Q.At(sright,k)))
					QGraph.Set(skEdge, k, 0.5*(Q.At(sleft, k)+Q.At(sright, k)))
					skSeg++
					skEdge++
					sleft = sright
				}
				// End point of range, excluding vertex
				// QGraph.Set((n+1)*NpEdge+n, k, Q.At(efi.InteriorPtsIndex[skSeg], k))
				QGraph.Set(skEdge, k, Q.At(sleft, k))
				skEdge++
				skSeg++
			}
			fmt.Println("skEdge, skSeg = ", skEdge, skSeg)
		}
		// if sk != Nefi {
		// 	fmt.Printf("Interior point edge count: %d, sk = %d\n", Nefi, sk)
		// 	panic("edge segment point count is not correct")
		// }
	}
}
