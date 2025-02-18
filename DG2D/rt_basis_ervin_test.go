package DG2D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"
)

func TestErvinBasis(t *testing.T) {
	var PMin, PMax int
	PMin = 1
	PMax = 1
	DivergencePolynomialField_RT1vsRTK_at1_Test(t, PMin, PMax)
}

func DivergencePolynomialField_RT1vsRTK_at1_Test(t *testing.T, PMin, PMax int) {
	BasisType := 0 // manual RT1
	// BasisType := 1 // RTK logic

	t.Log("Begin Divergence Test")
	// P := 1
	PStart := PMin
	PEnd := PMax
	for P := PStart; P <= PEnd; P++ {
		// PFieldStart := 0
		// PFieldEnd := P
		t.Logf("---------------------------------------------\n")
		t.Logf("Checking Divergence for RT%d\n", P)
		t.Logf("---------------------------------------------\n")
		rt := SetupRTTest(P)

		// Test the manual RT1 versus using the RTK logic
		e := NewErvinRTBasis(rt.P, rt.R, rt.S)
		edgeLocation := 2 * e.NpInt
		switch BasisType {
		case 0:
			g1 := rt.R.DataP[edgeLocation]
			g2 := rt.R.DataP[edgeLocation+1]
			e.Phi = e.ComposePhiRT1(g1, g2)
		// case 2:
		// 	g1 := rt.R.DataP[edgeLocation]
		// 	g2 := rt.R.DataP[edgeLocation+1]
		// 	g3 := rt.R.DataP[edgeLocation+2]
		// 	e.Phi = e.ComposePhiRT2(g1, g2, g3)
		case 1:
			e.InteriorPolyKBasis = NewJacobiBasis2D(e.P-1,
				rt.R.Copy().Subset(0, e.NpInt-1),
				rt.S.Copy().Subset(0, e.NpInt-1),
				0, 0)
			e.Phi = e.ComposePhiRTK(rt.R.DataP[edgeLocation : edgeLocation+e.NpEdge])
		}
		rt.Phi = e.Phi

		rt.V = rt.ComposeV(rt.Phi)
		rt.VInv = rt.V.InverseWithCheck()
		rt.Div = rt.ComputeDivergenceMatrix()
		rt.V.Print("V (Psi)")
		rt.Div.Print("Div (Psi)")

	}
}

func SetupRTTest(P int) (rt *RTElement) {
	var (
		Np     = (P + 1) * (P + 3)
		NpInt  = P * (P + 1) / 2 // Number of interior points is same as the 2D scalar space one order lesser
		NpEdge = P + 1           // Each edge is P+1 nodes
		oosr2  = 0.5 * math.Sqrt2
	)
	rt = &RTElement{
		P:          P,
		Np:         Np,
		NpInt:      NpInt,
		NpEdge:     NpEdge,
		V:          utils.NewMatrix(Np, Np),
		Div:        utils.NewMatrix(Np, Np),
		DOFVectors: make([]*ConstantVector, Np),
	}

	if P > 0 {
		if P < 9 {
			rt.RInt, rt.SInt = NodesEpsilon(P - 1)
		} else {
			rt.RInt, rt.SInt = XYtoRS(Nodes2D(P - 1))
			// Weak approach to pull back equidistant distro from the edges
			rt.RInt.Scale(0.93)
			rt.SInt.Scale(0.93)
			panic("This distribution is broken - TODO: fix\n")
		}
	}

	// Construct the unit vectors for the DOFs

	for i := 0; i < NpInt; i++ {
		rt.DOFVectors[i] = NewConstantVector(1, 0)
		rt.DOFVectors[i+NpInt] = NewConstantVector(0, 1)
	}

	offset := 2 * NpInt
	for i := 0; i < NpEdge; i++ {
		rt.DOFVectors[offset+i] = NewConstantVector(0, -1)
		rt.DOFVectors[offset+i+NpEdge] = NewConstantVector(oosr2, oosr2)
		rt.DOFVectors[offset+i+2*NpEdge] = NewConstantVector(-1, 0)
	}

	rt.R, rt.S = rt.ExtendGeomToRT(rt.RInt, rt.SInt)
	return
}
