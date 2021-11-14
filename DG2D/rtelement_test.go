package DG2D

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

func TestRTElement_CalculateBasis(t *testing.T) {
	{
		P := 1
		rtb := NewRTBasis2DSimplex(P, false)
		// Lagrange 1D polynomial
		{
			g1 := 0.5 - math.Sqrt(3)/6
			g2 := 0.5 + math.Sqrt(3)/6
			g1 = 2*g1 - 1
			g2 = 2*g2 - 1
			l1 := func(t float64) (p float64) {
				p = (t - g2) / (g1 - g2)
				return
			}
			l2 := func(t float64) (p float64) {
				p = (t - g1) / (g2 - g1)
				return
			}
			RGauss := DG1D.LegendreZeros(P)
			for _, r := range RGauss {
				assert.InDeltaf(t, l1(r), rtb.Lagrange1DPoly(r, RGauss, 0, RDir), 0.00001, "")
				assert.InDeltaf(t, l2(r), rtb.Lagrange1DPoly(r, RGauss, 1, RDir), 0.00001, "")
				assert.InDeltaf(t, 0, rtb.Lagrange1DPoly(r, RGauss, 0, RDir, Ds), 0.00001, "")
				assert.InDeltaf(t, 0, rtb.Lagrange1DPoly(r, RGauss, 1, RDir, Ds), 0.00001, "")
				// Check derivative
				assert.InDeltaf(t, 1.0/(g1-g2), rtb.Lagrange1DPoly(r, RGauss, 0, RDir, Dr), 0.00001, "")
				assert.InDeltaf(t, 1.0/(g2-g1), rtb.Lagrange1DPoly(r, RGauss, 1, RDir, Dr), 0.00001, "")
			}
		}
		// Linear 2D Polynomial
		{
			j := 0
			assert.InDeltaf(t, 0, rtb.LinearPoly(0, 0, j), 0.00001, "")
			assert.InDeltaf(t, 1, rtb.LinearPoly(-1, -1, j), 0.00001, "")
			assert.InDeltaf(t, 0, rtb.LinearPoly(-1, 1, j), 0.00001, "")
			assert.InDeltaf(t, 0, rtb.LinearPoly(1, -1, j), 0.00001, "")
			assert.InDeltaf(t, -1, rtb.LinearPoly(1, 1, j), 0.00001, "")
			assert.InDeltaf(t, 1./3, rtb.LinearPoly(-1./3, -1./3, j), 0.00001, "")

			j = 1
			assert.InDeltaf(t, 0.5, rtb.LinearPoly(0, 0, j), 0.00001, "")
			assert.InDeltaf(t, 0, rtb.LinearPoly(-1, -1, j), 0.00001, "")
			assert.InDeltaf(t, 0, rtb.LinearPoly(-1, 1, j), 0.00001, "")
			assert.InDeltaf(t, 1, rtb.LinearPoly(1, -1, j), 0.00001, "")
			assert.InDeltaf(t, 1, rtb.LinearPoly(1, 1, j), 0.00001, "")
			assert.InDeltaf(t, 1./3, rtb.LinearPoly(-1./3, -1./3, j), 0.00001, "")

			j = 2
			assert.InDeltaf(t, 0.5, rtb.LinearPoly(0, 0, j), 0.00001, "")
			assert.InDeltaf(t, 0, rtb.LinearPoly(-1, -1, j), 0.00001, "")
			assert.InDeltaf(t, 1, rtb.LinearPoly(-1, 1, j), 0.00001, "")
			assert.InDeltaf(t, 0, rtb.LinearPoly(1, -1, j), 0.00001, "")
			assert.InDeltaf(t, 1, rtb.LinearPoly(1, 1, j), 0.00001, "")
			assert.InDeltaf(t, 1./3, rtb.LinearPoly(-1./3, -1./3, j), 0.00001, "")

			// Dr
			assert.InDeltaf(t, -0.5, rtb.LinearPoly(0, 0, 0, Dr), 0.00001, "")
			assert.InDeltaf(t, 0.5, rtb.LinearPoly(0, 0, 1, Dr), 0.00001, "")
			assert.InDeltaf(t, 0.0, rtb.LinearPoly(0, 0, 2, Dr), 0.00001, "")

			// Ds
			assert.InDeltaf(t, -0.5, rtb.LinearPoly(0, 0, 0, Ds), 0.00001, "")
			assert.InDeltaf(t, 0.0, rtb.LinearPoly(0, 0, 1, Ds), 0.00001, "")
			assert.InDeltaf(t, 0.5, rtb.LinearPoly(0, 0, 2, Ds), 0.00001, "")
		}
	}
	if false {
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
	if false {
		P := 1
		rtb := NewRTBasis2DSimplex(P, false)
		assert.Equal(t, rtb.P, P)
		rtb.V[0].Print("V0")
		rtb.V[1].Print("V1")
		rtb.DivInt.Print("DivInt")
		R, S := NodesEpsilon(P - 1)
		rt := NewRTElement(R, S, P, false)
		rt.V[0].Print("V0 rt")
		rt.V[1].Print("V1 rt")
		rt.DivInt.Print("DivInt rt")
	}
}
