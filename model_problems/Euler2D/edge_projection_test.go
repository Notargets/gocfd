package Euler2D

import (
	"fmt"
	"testing"

	"github.com/stretchr/testify/assert"

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
	fmt.Println(c.DFR.Tris.EtoE)
	nen := c.NewNeighborNotifier()
	for k := 0; k < c.DFR.K; k++ {
		nen.PostNotification(0, k)
	}

	// for _, mp := range nen.mb.PostMsgQs[0] {
	// 	for _, msg := range mp.Cells() {
	// 		fmt.Printf("Message: Neighbor/From %d/%d\n", msg.KNeighbor,
	// 			msg.KRemote)
	// 	}
	// }
	nen.DeliverNotifications(0)
	for _, conn := range nen.ReadNotifications(0) {
		remoteFace := conn[0]
		remoteK := conn[1]
		localK := conn[2]
		// fmt.Printf("Edge[%d] Element[%d] ConnectsTo[%d]\n",
		// 	remoteFace, remoteK, localK)
		assert.Equal(t, c.DFR.Tris.EtoE[remoteK][remoteFace], localK)
	}
	// fmt.Println(nen.ReadNotifications(0))
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
	// fmt.Println("myK (global) for myKLocal", myK, myKLocal)
	// fmt.Println("Neighbors: ", Neighbors)
	// The message to each neighbor is simply to connect to this node via it's
	// local EtoE list of element->face mappings.
	for _, nbrK := range Neighbors { // For each neighboring tri
		if nbrK == -1 {
			continue
		}
		targetThread, _, _ := c.Partitions.GetBucket(nbrK)
		// fmt.Println("targetThread, nbrk = ", targetThread, nbrK)
		nen.mb.PostMessage(myThread, targetThread, &NeighborMsg{
			nbrK,
			myK,
		})
	}
}

func (nen *NeighborNotifier) DeliverNotifications(myThread int) {
	// Must be called in myThread before receivers can receive notifications
	nen.mb.DeliverMyMessages(myThread)
}

func (nen *NeighborNotifier) ReadNotifications(myThread int) (faces [][3]int) {
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
		for nEdge, neighbor := range c.DFR.Tris.EtoE[myK] {
			if neighbor == remoteK {
				faces = append(faces, [3]int{nEdge, myKLocal, neighbor})
				break
			}
		}
	}
	nen.mb.ClearMyMessages(myThread)
	return
}
