package Euler2D

import (
	"fmt"
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

func TestNeighborNotify(t *testing.T) {
	var (
		N = 2
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

	nen := c.NewNeighborNotifier()

	NpEdge := dfr.FluxElement.NpEdge
	var QInterp [4]utils.Matrix
	for n := 0; n < 4; n++ {
		QInterp[n] = dfr.GraphInterp.Mul(Q[n])
	}
	n := 0
	c.ShockFinder.UpdateShockedCells(Q[0])
	c.InterpolateSolutionToShockedEdges(c.ShockFinder, Q_Face, Q_Face_P0)
	for _, k := range c.ShockFinder.ShockCells.Cells() {
		nen.PostNotification(0, k)
		// QInterp has 1 extra points per edge
		for nEdge := 0; nEdge < 3; nEdge++ {
			off1 := nEdge * NpEdge
			off2 := nEdge * (NpEdge + 1)
			for i := 0; i < NpEdge; i++ {
				QInterp[n].Set(off2+i+1, k, Q_Face[n].At(off1+i, k))
			}
		}
	}
	fmt.Println("Here 1")
	for _, mp := range nen.mb.PostMsgQs[0] {
		for _, msg := range mp.Cells() {
			fmt.Println(*msg)
		}
	}
	nen.DeliverNotifications(0)
	fmt.Println("Here 2")
	fmt.Println(nen.ReadNotifications(0))
	// Replace vertices for plotting
	// dfr.AverageGraphFieldVertices(QInterp[n])
	// }
	// DG2D.PlotField(QInterp[0].Transpose().DataP, dfr.GraphMesh, 0.0, 0.0)
}

type NeighborMsg struct {
	KNeighbor, KRemote int
}

type NeighborNotifier struct {
	c  *Euler
	mb *utils.MailBox[*NeighborMsg]
}

func (c *Euler) NewNeighborNotifier() (nen *NeighborNotifier) {
	nen = &NeighborNotifier{
		c:  c,
		mb: utils.NewMailBox[*NeighborMsg](c.Partitions.ParallelDegree),
	}
	return
}

func (nen *NeighborNotifier) PostNotification(myThread, myKLocal int) {
	var (
		c         = nen.c
		myK       = c.Partitions.GetGlobalK(myKLocal, myThread)
		Neighbors = c.DFR.Tris.EtoE[myK]
	)
	// The message to each neighbor is simply to connect to this node via it's
	// local EtoE list of element->face mappings.
	for _, nbrK := range Neighbors { // For each neighboring tri
		targetThread, _, _ := c.Partitions.GetBucket(nbrK)
		nen.mb.PostMessage(myThread, targetThread, &NeighborMsg{
			nbrK,
			myK,
		})
		break
	}
}

func (nen *NeighborNotifier) DeliverNotifications(myThread int) {
	// Must be called in myThread before receivers can receive notifications
	nen.mb.DeliverMyMessages(myThread)
}

func (nen *NeighborNotifier) ReadNotifications(myThread int) (faces [][2]int) {
	var (
		c = nen.c
	)
	// Must be called after waiting for the DeliverNotifications is done in a
	// sync
	nen.mb.ReceiveMyMessages(myThread)
	for _, msg := range nen.mb.ReceiveMsgQs[myThread].Cells() {
		myK := msg.KNeighbor
		myKLocal, _, _ := c.Partitions.GetLocalK(myK)
		remoteK := msg.KRemote
		for nEdge, nbr := range c.DFR.Tris.EtoE[myK] {
			if nbr == remoteK {
				faces = append(faces, [2]int{nEdge, myKLocal})
				break
			}
		}
	}
	return
}
