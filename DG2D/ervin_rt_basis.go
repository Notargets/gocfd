package DG2D

import (
	"math"
)

type ErvinRTBasis struct {
	P int
}

func NewErvinRTBasis(P int) (e *ErvinRTBasis) {
	e = &ErvinRTBasis{
		P: P,
	}
	return
}

func (e *ErvinRTBasis) ComposePhi(t []float64) (phi []BasisPolynomialTerm) {
	// Compose the basis for the RT1 element based on Ervin
	var (
		// Get the two edge points for use in the lagrange terms for RT1
		sr2      = math.Sqrt2
		conv     = func(r float64) float64 { return (r + 1.) / 2. }
		dconvdrs = 1. / 2.
		Dot      = func(v1, v2 [2]float64) float64 { return v1[0]*v2[0] + v1[1]*v2[1] }
		Scale    = func(v [2]float64, scale float64) [2]float64 {
			return [2]float64{v[0] * scale, v[1] * scale}
		}
		Sum = func(v [2]float64) float64 { return v[0] + v[1] }
		e1v = func(r, s float64) (v [2]float64) {
			// Bottom edge
			xi, eta := conv(r), conv(s)
			v = [2]float64{xi, eta - 1.}
			// v = [2]float64{xi - 0.5, eta - 1.}
			return
		}
		e2v = func(r, s float64) (v [2]float64) {
			// Hypotenuse
			xi, eta := conv(r), conv(s)
			v = [2]float64{sr2 * xi, sr2 * eta}
			return
		}
		e3v = func(r, s float64) (v [2]float64) {
			// Left edge
			xi, eta := conv(r), conv(s)
			v = [2]float64{xi - 1., eta}
			// v = [2]float64{xi - 1., eta - 0.5}
			return
		}
		e4v = func(r, s float64) (v [2]float64) {
			xi, eta := conv(r), conv(s)
			v = [2]float64{eta * xi, eta * (eta - 1.)}
			return
		}
		e5v = func(r, s float64) (v [2]float64) {
			xi, eta := conv(r), conv(s)
			v = [2]float64{xi * (xi - 1.), eta * xi}
			return
		}
		E4Vector = BasisVectorStruct{
			// Interior vector E4
			Eval: func(r, s float64) [2]float64 { return e4v(r, s) },
			Dot: func(r, s float64, f [2]float64) float64 {
				return Dot(e4v(r, s), f)
			},
			Project: func(r, s float64, scale float64) [2]float64 {
				return Scale(e4v(r, s), scale)
			},
			Divergence: func(r, s float64) (div float64) {
				div = (3.*s + 1.) / 4.
				return
			},
			Sum: func(r, s float64) float64 { return Sum(e4v(r, s)) },
		}
		E5Vector = BasisVectorStruct{
			// Interior vector E5
			Eval: func(r, s float64) [2]float64 { return e5v(r, s) },
			Dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e5v(r, s), f)
			},
			Project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e5v(r, s), scale)
			},
			Divergence: func(r, s float64) (div float64) {
				div = (3.*r + 1.) / 4.
				return
			},
			Sum: func(r, s float64) float64 { return Sum(e5v(r, s)) },
		}
		Edge1Vector = BasisVectorStruct{
			// Bottom edge
			Eval: func(r, s float64) (v [2]float64) { return e1v(r, s) },
			Dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e1v(r, s), f)
			},
			Project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e1v(r, s), scale)
			},
			Divergence: func(r, s float64) (div float64) { return 1. },
			Sum:        func(r, s float64) float64 { return Sum(e1v(r, s)) },
		}
		Edge2Vector = BasisVectorStruct{
			// Hypotenuse
			Eval: func(r, s float64) (v [2]float64) { return e2v(r, s) },
			Dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e2v(r, s), f)
			},
			Project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e2v(r, s), scale)
			},
			Divergence: func(r, s float64) (div float64) { return sr2 },
			Sum:        func(r, s float64) float64 { return Sum(e2v(r, s)) },
		}
		Edge3Vector = BasisVectorStruct{
			// Left edge
			Eval: func(r, s float64) (v [2]float64) { return e3v(r, s) },
			Dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e3v(r, s), f)
			},
			Project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e3v(r, s), scale)
			},
			Divergence: func(r, s float64) (div float64) { return 1. },
			Sum:        func(r, s float64) float64 { return Sum(e3v(r, s)) },
		}
		l1Func = func(tt float64) float64 {
			var (
				g1, g2 = conv(t[0]), conv(t[1])
			)
			return (tt - g2) / (g1 - g2)
		}
		l1Deriv = func() float64 {
			var (
				g1, g2 = conv(t[0]), conv(t[1])
			)
			return dconvdrs / (g1 - g2)
		}
		l2Func = func(tt float64) float64 {
			var (
				g1, g2 = conv(t[0]), conv(t[1])
			)
			return (tt - g1) / (g2 - g1)
		}
		l2Deriv = func() float64 { return -l1Deriv() }
		l1xi    = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(conv(r)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{l1Deriv(), 0}
			},
			OrderOfTerm: 1,
		}
		l1etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(conv(s)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, l1Deriv()}
			},
			OrderOfTerm: 1,
		}
		l1etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(conv(-s)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, -l1Deriv()}
			},
			OrderOfTerm: 1,
		}
		l2xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(conv(r)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{l2Deriv(), 0}
			},
			OrderOfTerm: 1,
		}
		l2etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(conv(s)) },
			Gradient: func(r, s float64) [2]float64 {
				// return [2]float64{-l2Deriv(), l2Deriv()}
				return [2]float64{0, l2Deriv()}
			},
			OrderOfTerm: 1,
		}
		l2etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(conv(-s)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, -l2Deriv()}
			},
			OrderOfTerm: 1,
		}
		q1Func = func(tt float64) float64 {
			var (
				g1_q, g2_q, g3_q = conv(t[0]), conv(t[1]), conv(t[2])
			)
			return (tt - g2_q) * (tt - g3_q) / ((g1_q - g2_q) * (g1_q - g3_q))
		}
		q1Deriv = func(tt float64) (deriv float64) {
			var (
				g1_q, g2_q, g3_q = conv(t[0]), conv(t[1]), conv(t[2])
			)
			// return dconvdrs / (g1 - g2)
			deriv = ((tt - g2_q) + (tt - g3_q)) / ((g1_q - g2_q) * (g1_q - g3_q))
			deriv += dconvdrs
			return
		}
		q1xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q1Func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{q1Deriv(r), 0}
				return
			},
			OrderOfTerm: 2,
		}
		q1etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q1Func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, q1Deriv(s)}
				return
			},
			OrderOfTerm: 2,
		}
		q2Func = func(tt float64) float64 {
			var (
				g1_q, g2_q, g3_q = conv(t[0]), conv(t[1]), conv(t[2])
			)
			return (tt - g1_q) * (tt - g3_q) / ((g2_q - g1_q) * (g2_q - g3_q))
		}
		q2Deriv = func(tt float64) (deriv float64) {
			var (
				g1_q, g2_q, g3_q = conv(t[0]), conv(t[1]), conv(t[2])
			)
			deriv = ((tt - g1_q) + (tt - g3_q)) / ((g2_q - g1_q) * (g2_q - g3_q))
			deriv += dconvdrs
			return
		}
		q2xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q2Func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{q2Deriv(r), 0}
				return
			},
			OrderOfTerm: 2,
		}
		q2etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q2Func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, q2Deriv(s)}
				return
			},
			OrderOfTerm: 2,
		}
		q3Func = func(tt float64) float64 {
			var (
				g1_q, g2_q, g3_q = conv(t[0]), conv(t[1]), conv(t[2])
			)
			return (tt - g1_q) * (tt - g2_q) / ((g3_q - g1_q) * (g3_q - g2_q))
		}
		q3Deriv = func(tt float64) (deriv float64) {
			var (
				g1_q, g2_q, g3_q = conv(t[0]), conv(t[1]), conv(t[2])
			)
			deriv = ((tt - g1_q) + (tt - g2_q)) / ((g3_q - g1_q) * (g3_q - g2_q))
			deriv += dconvdrs
			return
		}
		q3xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q3Func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{q3Deriv(r), 0}
				return
			},
			OrderOfTerm: 2,
		}
		q3etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q3Func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, q3Deriv(s)}
				return
			},
			OrderOfTerm: 2,
		}
		q1etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q1Func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -q1Deriv(s)}
				return
			},
			OrderOfTerm: 2,
		}
		q2etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q2Func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -q2Deriv(s)}
				return
			},
			OrderOfTerm: 2,
		}
		q3etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return q3Func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -q3Deriv(s)}
				return
			},
			OrderOfTerm: 2,
		}
		e45mult1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return 1. - conv(r) - conv(s)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{-0.5, -0.5}
			},
			OrderOfTerm: 1,
		}
		e45mult2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return conv(r)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{0.5, 0}
			},
			OrderOfTerm: 1,
		}
		e45mult3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return conv(s)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{0, 0.5}
			},
			OrderOfTerm: 1,
		}
		constant = BasisPolynomialMultiplier{
			Eval:        func(r, s float64) float64 { return 1. },
			Gradient:    func(r, s float64) (grad [2]float64) { return },
			OrderOfTerm: 0,
		}
	)
	switch e.P {
	case 1:
		phi = []BasisPolynomialTerm{
			{PolyMultiplier: constant, BasisVector: E4Vector},
			{PolyMultiplier: constant, BasisVector: E5Vector},
			{PolyMultiplier: l1xi, BasisVector: Edge1Vector},
			{PolyMultiplier: l2xi, BasisVector: Edge1Vector},
			{PolyMultiplier: l1etaEdge2, BasisVector: Edge2Vector},
			{PolyMultiplier: l2etaEdge2, BasisVector: Edge2Vector},
			{PolyMultiplier: l2etaEdge3, BasisVector: Edge3Vector},
			{PolyMultiplier: l1etaEdge3, BasisVector: Edge3Vector},
		}
	case 2:
		phi = []BasisPolynomialTerm{
			{PolyMultiplier: e45mult1, BasisVector: E4Vector},
			{PolyMultiplier: e45mult2, BasisVector: E4Vector},
			{PolyMultiplier: e45mult3, BasisVector: E4Vector},
			{PolyMultiplier: e45mult1, BasisVector: E5Vector},
			{PolyMultiplier: e45mult2, BasisVector: E5Vector},
			{PolyMultiplier: e45mult3, BasisVector: E5Vector},
			{PolyMultiplier: q1xi, BasisVector: Edge1Vector},
			{PolyMultiplier: q2xi, BasisVector: Edge1Vector},
			{PolyMultiplier: q3xi, BasisVector: Edge1Vector},
			{PolyMultiplier: q1etaEdge2, BasisVector: Edge2Vector},
			{PolyMultiplier: q2etaEdge2, BasisVector: Edge2Vector},
			{PolyMultiplier: q3etaEdge2, BasisVector: Edge2Vector},
			{PolyMultiplier: q3etaEdge3, BasisVector: Edge3Vector},
			{PolyMultiplier: q2etaEdge3, BasisVector: Edge3Vector},
			{PolyMultiplier: q1etaEdge3, BasisVector: Edge3Vector},
		}
	default:
		panic("Order unsupported for Ervin basis methods")
	}

	return
}
