package DG2D

import (
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

func TestRTElement_CalculateBasis(t *testing.T) {
	{
		var rtb *RTBasis2DSimplex
		var Rint, Sint utils.Vector
		P := 1
		rtb = &RTBasis2DSimplex{
			P:      P,
			Np:     (P + 1) * (P + 3),
			NpInt:  P * (P + 1) / 2, // Number of interior points is same as the 2D scalar space one order lesser
			NpEdge: P + 1,           // Each edge is P+1 nodes
		}
		if rtb.P > 0 {
			Rint, Sint = NodesEpsilon(rtb.P - 1)
		}
		rtb.R, rtb.S = rtb.ExtendGeomToRT(Rint, Sint)
		RGauss := DG1D.LegendreZeros(rtb.P)
		//for i := 0; i < rtb.Np; i++ {
		i := 0
		r, s := rtb.R.DataP[i], rtb.S.DataP[i]

		// Core basis
		p0, p1 := rtb.getCoreBasisTerm(e1, r, s)
		assert.InDeltaf(t, 0.471, p0, 0.001, "")
		assert.InDeltaf(t, 0.471, p1, 0.001, "")
		p0, p1 = rtb.getCoreBasisTerm(e2, r, s)
		assert.InDeltaf(t, -0.667, p0, 0.001, "")
		assert.InDeltaf(t, 0.333, p1, 0.001, "")
		p0, p1 = rtb.getCoreBasisTerm(e3, r, s)
		assert.InDeltaf(t, 0.333, p0, 0.001, "")
		assert.InDeltaf(t, -0.667, p1, 0.001, "")

		// Lagrange 1D basis
		l1xi := rtb.Lagrange1DPoly(r, RGauss, 0, RDir)
		l1eta := rtb.Lagrange1DPoly(s, RGauss, 0, SDir)
		assert.InDeltaf(t, 0.789, l1xi, 0.001, "")
		assert.InDeltaf(t, 0.789, l1eta, 0.001, "")
		l2xi := rtb.Lagrange1DPoly(r, RGauss, 1, RDir)
		l2eta := rtb.Lagrange1DPoly(s, RGauss, 1, SDir)
		assert.InDeltaf(t, 0.211, l2xi, 0.001, "")
		assert.InDeltaf(t, 0.211, l2eta, 0.001, "")

		// Core basis derivatives
		p0, p1 = rtb.getCoreBasisTerm(e1, r, s, Dr)
		assert.InDeltaf(t, 0.707, p0, 0.001, "")
		assert.InDeltaf(t, 0.000, p1, 0.001, "")
		p0, p1 = rtb.getCoreBasisTerm(e2, r, s, Dr)
		assert.InDeltaf(t, 0.500, p0, 0.001, "")
		assert.InDeltaf(t, 0.000, p1, 0.001, "")
		p0, p1 = rtb.getCoreBasisTerm(e3, r, s, Dr)
		assert.InDeltaf(t, 0.500, p0, 0.001, "")
		assert.InDeltaf(t, 0.000, p1, 0.001, "")
		p0, p1 = rtb.getCoreBasisTerm(e1, r, s, Ds)
		assert.InDeltaf(t, 0.000, p0, 0.001, "")
		assert.InDeltaf(t, 0.707, p1, 0.001, "")
		p0, p1 = rtb.getCoreBasisTerm(e2, r, s, Ds)
		assert.InDeltaf(t, 0.000, p0, 0.001, "")
		assert.InDeltaf(t, 0.500, p1, 0.001, "")
		p0, p1 = rtb.getCoreBasisTerm(e3, r, s, Ds)
		assert.InDeltaf(t, 0.000, p0, 0.001, "")
		assert.InDeltaf(t, 0.500, p1, 0.001, "")

		// Lagrange 1D basis derivatives
		l1xiDr := rtb.Lagrange1DPoly(r, RGauss, 0, RDir, Dr)
		l1etaDr := rtb.Lagrange1DPoly(s, RGauss, 0, SDir, Dr)
		assert.InDeltaf(t, -0.433, l1xiDr, 0.001, "")
		assert.InDeltaf(t, 0.000, l1etaDr, 0.001, "")
		l1xiDs := rtb.Lagrange1DPoly(r, RGauss, 0, RDir, Ds)
		l1etaDs := rtb.Lagrange1DPoly(s, RGauss, 0, SDir, Ds)
		assert.InDeltaf(t, 0.000, l1xiDs, 0.001, "")
		assert.InDeltaf(t, -0.433, l1etaDs, 0.001, "")

		l2xiDr := rtb.Lagrange1DPoly(r, RGauss, 1, RDir, Dr)
		l2etaDr := rtb.Lagrange1DPoly(s, RGauss, 1, SDir, Dr)
		assert.InDeltaf(t, 0.433, l2xiDr, 0.001, "")
		assert.InDeltaf(t, 0.000, l2etaDr, 0.001, "")
		l2xiDs := rtb.Lagrange1DPoly(r, RGauss, 1, RDir, Ds)
		l2etaDs := rtb.Lagrange1DPoly(s, RGauss, 1, SDir, Ds)
		assert.InDeltaf(t, 0.000, l2xiDs, 0.001, "")
		assert.InDeltaf(t, 0.433, l2etaDs, 0.001, "")
	}
	if true {
		rtb := NewRTBasis2DSimplex(1, false)
		assert.Equal(t, rtb.P, 1)
	}
}
