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
		sr2   = math.Sqrt2
		conv  = func(r float64) float64 { return (r + 1.) / 2. }
		Dot   = func(v1, v2 [2]float64) float64 { return v1[0]*v2[0] + v1[1]*v2[1] }
		Scale = func(v [2]float64, scale float64) [2]float64 {
			return [2]float64{v[0] * scale, v[1] * scale}
		}
		Sum = func(v [2]float64) float64 { return v[0] + v[1] }
		// normalize = func(v [2]float64) (normed [2]float64) {
		// 	norm := 1. / math.Sqrt(v[0]*v[0]+v[1]*v[1])
		// 	normed = [2]float64{v[0] * norm, v[1] * norm}
		// 	return
		// }
		e1v = func(r, s float64) (v [2]float64) {
			// Bottom edge
			xi, eta := conv(r), conv(s)
			// v = [2]float64{xi, eta - 1.}
			v = [2]float64{xi - 0.5, eta - 1.}
			// v = [2]float64{0, -1}
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
			// v = [2]float64{xi - 1., eta}
			v = [2]float64{xi - 1., eta - 0.5}
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
			Divergence: func(r, s float64) (div float64) { return },
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
			Divergence: func(r, s float64) (div float64) { return },
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
			Divergence: func(r, s float64) (div float64) { return },
			Sum:        func(r, s float64) float64 { return Sum(e3v(r, s)) },
		}
		l1xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1, g2 = t[0], t[1]
				)
				return (r - g2) / (g1 - g2)
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1, g2 = t[0], t[1]
				)
				div = 1 / (g1 - g2)
				return
			},
			OrderOfTerm: 1,
		}
		l1eta = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1, g2 = t[0], t[1]
				)
				return (s - g2) / (g1 - g2)
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1, g2 = t[0], t[1]
				)
				div = 1 / (g1 - g2)
				return
			},
			OrderOfTerm: 1,
		}
		l2xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1, g2 = t[0], t[1]
				)
				return (r - g1) / (g2 - g1)
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1, g2 = t[0], t[1]
				)
				div = 1 / (g2 - g1)
				return
			},
			OrderOfTerm: 1,
		}
		l2eta = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1, g2 = t[0], t[1]
				)
				return (s - g1) / (g2 - g1)
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1, g2 = t[0], t[1]
				)
				div = 1 / (g2 - g1)
				return
			},
			OrderOfTerm: 1,
		}
		q1xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				return (r - g2_q) * (r - g3_q) / ((g1_q - g2_q) * (g1_q - g3_q))
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				div = ((r - g2_q) + (r - g3_q)) / ((g1_q - g2_q) * (g1_q - g3_q))
				return
			},
			OrderOfTerm: 2,
		}
		q1eta = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				return (s - g2_q) * (s - g3_q) / ((g1_q - g2_q) * (g1_q - g3_q))
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				div = ((s - g2_q) + (s - g3_q)) / ((g1_q - g2_q) * (g1_q - g3_q))
				return
			},
			OrderOfTerm: 2,
		}
		q2xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				return (r - g1_q) * (r - g3_q) / ((g2_q - g1_q) * (g2_q - g3_q))
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				div = ((r - g1_q) + (r - g3_q)) / ((g2_q - g1_q) * (g2_q - g3_q))
				return
			},
			OrderOfTerm: 2,
		}
		q2eta = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				return (s - g1_q) * (s - g3_q) / ((g2_q - g1_q) * (g2_q - g3_q))
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				div = ((s - g1_q) + (s - g3_q)) / ((g2_q - g1_q) * (g2_q - g3_q))
				return
			},
			OrderOfTerm: 2,
		}
		q3xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				return (r - g1_q) * (r - g2_q) / ((g3_q - g1_q) * (g3_q - g2_q))
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				div = ((r - g1_q) + (r - g2_q)) / ((g3_q - g1_q) * (g3_q - g2_q))
				return
			},
			OrderOfTerm: 2,
		}
		q3eta = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				return (s - g1_q) * (s - g2_q) / ((g3_q - g1_q) * (g3_q - g2_q))
			},
			Divergence: func(r, s float64) (div float64) {
				var (
					g1_q, g2_q, g3_q = t[0], t[1], t[2]
				)
				div = ((s - g1_q) + (s - g2_q)) / ((g3_q - g1_q) * (g3_q - g2_q))
				return
			},
			OrderOfTerm: 2,
		}
		e45mult1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return 1. - conv(r) - conv(s)
			},
			Divergence: func(r, s float64) (div float64) {
				return -1.
			},
			OrderOfTerm: 1,
		}
		e45mult2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return conv(r)
			},
			Divergence: func(r, s float64) (div float64) {
				return 0.5
			},
			OrderOfTerm: 1,
		}
		e45mult3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return conv(s)
			},
			Divergence: func(r, s float64) (div float64) {
				return 0.5
			},
			OrderOfTerm: 1,
		}
		constant = BasisPolynomialMultiplier{
			Eval:        func(r, s float64) float64 { return 1. },
			Divergence:  func(r, s float64) (div float64) { return 0. },
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
			{PolyMultiplier: l1eta, BasisVector: Edge2Vector},
			{PolyMultiplier: l2eta, BasisVector: Edge2Vector},
			{PolyMultiplier: l2eta, BasisVector: Edge3Vector},
			{PolyMultiplier: l1eta, BasisVector: Edge3Vector},
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
			{PolyMultiplier: q1eta, BasisVector: Edge2Vector},
			{PolyMultiplier: q2eta, BasisVector: Edge2Vector},
			{PolyMultiplier: q3eta, BasisVector: Edge2Vector},
			{PolyMultiplier: q3eta, BasisVector: Edge3Vector},
			{PolyMultiplier: q2eta, BasisVector: Edge3Vector},
			{PolyMultiplier: q1eta, BasisVector: Edge3Vector},
		}
	default:
		panic("Order unsupported for Ervin basis methods")
	}
	return
}
