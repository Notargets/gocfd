package DG2D

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/DG1D"
)

func TestRTElement_CalculateBasis(t *testing.T) {
	{
		P := 1
		rtb := NewRTBasis2DSimplex(P)
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
		// Core basis functions e1-e5
		{
			var r, s, p0, p1 float64
			sr2 := math.Sqrt(2)
			// e1
			r, s = -1., -1.
			p0, p1 = rtb.getCoreBasisTerm(e1, r, s)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			r, s = 1., 1.
			p0, p1 = rtb.getCoreBasisTerm(e1, r, s)
			assert.InDeltaf(t, sr2, p0, 0.000001, "")
			assert.InDeltaf(t, sr2, p1, 0.000001, "")
			r, s = 0., 0.
			p0, p1 = rtb.getCoreBasisTerm(e1, r, s)
			assert.InDeltaf(t, 0.5*sr2, p0, 0.000001, "")
			assert.InDeltaf(t, 0.5*sr2, p1, 0.000001, "")

			// e2
			r, s = -1., -1.
			p0, p1 = rtb.getCoreBasisTerm(e2, r, s)
			assert.InDeltaf(t, -1, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			r, s = 1., 1.
			p0, p1 = rtb.getCoreBasisTerm(e2, r, s)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 1, p1, 0.000001, "")
			r, s = 0., 0.
			p0, p1 = rtb.getCoreBasisTerm(e2, r, s)
			assert.InDeltaf(t, -.5, p0, 0.000001, "")
			assert.InDeltaf(t, .5, p1, 0.000001, "")

			r, s = -1., -1.
			p0, p1 = rtb.getCoreBasisTerm(e3, r, s)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, -1, p1, 0.000001, "")
			r, s = 1., 1.
			p0, p1 = rtb.getCoreBasisTerm(e3, r, s)
			assert.InDeltaf(t, 1, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			r, s = 0., 0.
			p0, p1 = rtb.getCoreBasisTerm(e3, r, s)
			assert.InDeltaf(t, .5, p0, 0.000001, "")
			assert.InDeltaf(t, -.5, p1, 0.000001, "")

			r, s = -1., -1.
			p0, p1 = rtb.getCoreBasisTerm(e4, r, s)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			r, s = 1., 1.
			p0, p1 = rtb.getCoreBasisTerm(e4, r, s)
			assert.InDeltaf(t, 1, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			r, s = 0., 0.
			p0, p1 = rtb.getCoreBasisTerm(e4, r, s)
			assert.InDeltaf(t, .25, p0, 0.000001, "")
			assert.InDeltaf(t, -.25, p1, 0.000001, "")

			r, s = -1., -1.
			p0, p1 = rtb.getCoreBasisTerm(e5, r, s)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			r, s = 1., 1.
			p0, p1 = rtb.getCoreBasisTerm(e5, r, s)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 1, p1, 0.000001, "")
			r, s = 0., 0.
			p0, p1 = rtb.getCoreBasisTerm(e5, r, s)
			assert.InDeltaf(t, -.25, p0, 0.000001, "")
			assert.InDeltaf(t, .25, p1, 0.000001, "")

			r, s = -1, -1
			p0, p1 = rtb.getCoreBasisTerm(e1, r, s, Dr)
			assert.InDeltaf(t, 0.5*sr2, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			p0, p1 = rtb.getCoreBasisTerm(e2, r, s, Dr)
			assert.InDeltaf(t, 0.5, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			p0, p1 = rtb.getCoreBasisTerm(e3, r, s, Dr)
			assert.InDeltaf(t, 0.5, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			p0, p1 = rtb.getCoreBasisTerm(e4, r, s, Dr)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
			p0, p1 = rtb.getCoreBasisTerm(e5, r, s, Dr)
			assert.InDeltaf(t, -.5, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")

			p0, p1 = rtb.getCoreBasisTerm(e1, r, s, Ds)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 0.5*sr2, p1, 0.000001, "")
			p0, p1 = rtb.getCoreBasisTerm(e2, r, s, Ds)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 0.5, p1, 0.000001, "")
			p0, p1 = rtb.getCoreBasisTerm(e3, r, s, Ds)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 0.5, p1, 0.000001, "")
			p0, p1 = rtb.getCoreBasisTerm(e4, r, s, Ds)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, -0.5, p1, 0.000001, "")
			p0, p1 = rtb.getCoreBasisTerm(e5, r, s, Ds)
			assert.InDeltaf(t, 0, p0, 0.000001, "")
			assert.InDeltaf(t, 0, p1, 0.000001, "")
		}
	}
	if true {
		P := 2
		rtb := NewRTBasis2DSimplex(P)
		assert.Equal(t, rtb.P, P)
		rt := NewRTElement(P)
		assert.InDeltaSlicef(t, rt.V[0].DataP, rtb.V[0].DataP, 0.000001, "")
		assert.InDeltaSlicef(t, rt.V[1].DataP, rtb.V[1].DataP, 0.000001, "")
		rtb.DivInt.Print("DivInt")
		rt.DivInt.Print("DivInt rt")
		// assert.InDeltaSlicef(t, rt.DivInt.DataP, rtb.DivInt.DataP, 0.000001, "")
		// rtb.Div.Print("Div")
		// rt.Div.Print("Div rt")
	}
}
