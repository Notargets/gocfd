package DG2D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type ErvinRTBasis struct {
	P                  int
	Np, NpInt, NpEdge  int
	E4Vector           BasisVectorStruct
	E5Vector           BasisVectorStruct
	Edge1Vector        BasisVectorStruct
	Edge2Vector        BasisVectorStruct
	Edge3Vector        BasisVectorStruct
	Phi                []BasisPolynomialTerm
	InteriorPolyKBasis *JacobiBasis2D
}

func NewErvinRTBasis(P int, R, S utils.Vector) (e *ErvinRTBasis) {
	var (
		sr2   = math.Sqrt2
		Dot   = func(v1, v2 [2]float64) float64 { return v1[0]*v2[0] + v1[1]*v2[1] }
		Scale = func(v [2]float64, scale float64) [2]float64 {
			return [2]float64{v[0] * scale, v[1] * scale}
		}
		e1v = func(r, s float64) (v [2]float64) {
			// Bottom edge
			xi, eta := e.conv(r), e.conv(s)
			v = [2]float64{xi, eta - 1.}
			// v = [2]float64{xi - 0.5, eta - 1.}
			return
		}
		e2v = func(r, s float64) (v [2]float64) {
			// Hypotenuse
			xi, eta := e.conv(r), e.conv(s)
			v = [2]float64{sr2 * xi, sr2 * eta}
			return
		}
		e3v = func(r, s float64) (v [2]float64) {
			// Left edge
			xi, eta := e.conv(r), e.conv(s)
			v = [2]float64{xi - 1., eta}
			// v = [2]float64{xi - 1., eta - 0.5}
			return
		}
		e4v = func(r, s float64) (v [2]float64) {
			xi, eta := e.conv(r), e.conv(s)
			v = [2]float64{eta * xi, eta * (eta - 1.)}
			return
		}
		e5v = func(r, s float64) (v [2]float64) {
			xi, eta := e.conv(r), e.conv(s)
			v = [2]float64{xi * (xi - 1.), eta * xi}
			return
		}
	)

	e = &ErvinRTBasis{
		P:      P,
		Np:     (P + 1) * (P + 3),
		NpInt:  P * (P + 1) / 2,
		NpEdge: P + 1,
		E4Vector: BasisVectorStruct{
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
		},
		E5Vector: BasisVectorStruct{
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
		},
		Edge1Vector: BasisVectorStruct{
			// Bottom edge
			Eval: func(r, s float64) (v [2]float64) { return e1v(r, s) },
			Dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e1v(r, s), f)
			},
			Project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e1v(r, s), scale)
			},
			Divergence: func(r, s float64) (div float64) { return 1. },
		},
		Edge2Vector: BasisVectorStruct{
			// Hypotenuse
			Eval: func(r, s float64) (v [2]float64) { return e2v(r, s) },
			Dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e2v(r, s), f)
			},
			Project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e2v(r, s), scale)
			},
			Divergence: func(r, s float64) (div float64) { return sr2 },
		},
		Edge3Vector: BasisVectorStruct{
			// Left edge
			Eval: func(r, s float64) (v [2]float64) { return e3v(r, s) },
			Dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e3v(r, s), f)
			},
			Project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e3v(r, s), scale)
			},
			Divergence: func(r, s float64) (div float64) { return 1. },
		},
	}
	edgeLocation := 2 * e.NpInt
	switch e.P {
	case 1:
		// e.g1, e.g2 = conv(SingleEdgePts[0]), conv(SingleEdgePts[1])
		// SingleEdgePoints are the P+1 points along each edge, the same for all
		// rt.R.DataP[edgeLocation:edgeLocation+rt.NpEdge])
		g1 := R.DataP[edgeLocation]
		g2 := R.DataP[edgeLocation+1]
		e.Phi = e.ComposePhiRT1(g1, g2)
	case 2:
		g1 := R.DataP[edgeLocation]
		g2 := R.DataP[edgeLocation+1]
		g3 := R.DataP[edgeLocation+2]
		e.Phi = e.ComposePhiRT2(g1, g2, g3)
	default:
		e.InteriorPolyKBasis = NewJacobiBasis2D(e.P-1, R, S, 0, 0)
	}
	return
}

func (e *ErvinRTBasis) ComposePhiRTK() (phi []BasisPolynomialTerm) {
	// phi = []BasisPolynomialTerm{
	// 	{PolyMultiplier: constant, BasisVector: e.E4Vector},
	// 	{PolyMultiplier: constant, BasisVector: e.E5Vector},
	// 	{PolyMultiplier: l1xiEdge1, BasisVector: e.Edge1Vector},
	// 	{PolyMultiplier: l2xiEdge1, BasisVector: e.Edge1Vector},
	// 	{PolyMultiplier: l1etaEdge2, BasisVector: e.Edge2Vector},
	// 	{PolyMultiplier: l2etaEdge2, BasisVector: e.Edge2Vector},
	// 	{PolyMultiplier: l2etaEdge3, BasisVector: e.Edge3Vector},
	// 	{PolyMultiplier: l1etaEdge3, BasisVector: e.Edge3Vector},
	// }
	return
}

func (e *ErvinRTBasis) ComposePhiRT1(g1, g2 float64) (phi []BasisPolynomialTerm) {
	// Compose the basis for the RT1 element based on Ervin
	var (
		constant = BasisPolynomialMultiplier{
			Eval:     func(r, s float64) float64 { return 1. },
			Gradient: func(r, s float64) (grad [2]float64) { return },
		}

		l1Func = func(t float64) float64 {
			return (t - g2) / (g1 - g2)
		}
		l1Deriv = func() float64 {
			return e.dconvDrDs() / (g1 - g2)
		}

		l2Func = func(t float64) float64 {
			return (t - g1) / (g2 - g1)
		}
		l2Deriv = func() float64 { return -l1Deriv() }

		l1xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(e.conv(r)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{l1Deriv(), 0}
			},
		}
		l2xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(e.conv(r)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{l2Deriv(), 0}
			},
		}

		l1etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(e.conv(s)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, l1Deriv()}
			},
		}
		l2etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(e.conv(s)) },
			Gradient: func(r, s float64) [2]float64 {
				// return [2]float64{-l2Deriv(), l2Deriv()}
				return [2]float64{0, l2Deriv()}
			},
		}

		l1etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(e.conv(-s)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, -l1Deriv()}
			},
		}
		l2etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(e.conv(-s)) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, -l2Deriv()}
			},
		}
	)
	phi = []BasisPolynomialTerm{
		{PolyMultiplier: constant, BasisVector: e.E4Vector},
		{PolyMultiplier: constant, BasisVector: e.E5Vector},
		{PolyMultiplier: l1xiEdge1, BasisVector: e.Edge1Vector},
		{PolyMultiplier: l2xiEdge1, BasisVector: e.Edge1Vector},
		{PolyMultiplier: l1etaEdge2, BasisVector: e.Edge2Vector},
		{PolyMultiplier: l2etaEdge2, BasisVector: e.Edge2Vector},
		{PolyMultiplier: l2etaEdge3, BasisVector: e.Edge3Vector},
		{PolyMultiplier: l1etaEdge3, BasisVector: e.Edge3Vector},
	}

	return
}

func (e *ErvinRTBasis) conv(r float64) float64 { return (r + 1.) / 2. }

func (e *ErvinRTBasis) dconvDrDs() float64 { return 1. / 2. }

func (e *ErvinRTBasis) lagrange1D(r, s float64, t []float64, j int) (val float64) {
	var (
		div float64
	)
	for i := 0; i < e.NpEdge; i++ {

	}
	_ = div
	return
}

func (e *ErvinRTBasis) ComposePhiRT2(g1, g2, g3 float64) (phi []BasisPolynomialTerm) {
	// Compose the basis for the RT1 element based on Ervin
	var (
		// Get the two edge points for use in the lagrange terms for RT1
		div1 = 1. / ((g1 - g2) * (g1 - g3))
		div2 = 1. / ((g2 - g1) * (g2 - g3))
		div3 = 1. / ((g3 - g1) * (g3 - g2))

		l1func = func(t float64) float64 {
			return (t - g2) * (t - g3) * div1
		}
		l2Func = func(t float64) float64 {
			return (t - g1) * (t - g3) * div2
		}
		l3Func = func(t float64) float64 {
			return (t - g1) * (t - g2) * div3
		}

		l1Deriv = func(t float64) (deriv float64) {
			deriv = ((t - g2) + (t - g3)) * div1
			deriv *= e.dconvDrDs()
			return
		}
		l2Deriv = func(t float64) (deriv float64) {
			deriv = ((t - g1) + (t - g3)) * div2
			deriv *= e.dconvDrDs()
			return
		}
		l3Deriv = func(t float64) (deriv float64) {
			deriv = ((t - g1) + (t - g2)) * div3
			deriv *= e.dconvDrDs()
			return
		}

		l1xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(e.conv(r)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l1Deriv(e.conv(r)), 0}
				return
			},
		}
		l2xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(e.conv(r)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l2Deriv(e.conv(r)), 0}
				return
			},
		}
		l3xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(e.conv(r)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l3Deriv(e.conv(r)), 0}
				return
			},
		}

		l1etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(e.conv(s)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l1Deriv(e.conv(s))}
				return
			},
		}
		l2etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(e.conv(s)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l2Deriv(e.conv(s))}
				return
			},
		}
		l3etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(e.conv(s)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l3Deriv(e.conv(s))}
				return
			},
		}

		l1etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(e.conv(-s)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l1Deriv(e.conv(-s))}
				return
			},
		}
		l2etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(e.conv(-s)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l2Deriv(e.conv(-s))}
				return
			},
		}
		l3etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(e.conv(-s)) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l3Deriv(e.conv(-s))}
				return
			},
		}

		e45mult1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return 1. - e.conv(r) - e.conv(s)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{-0.5, -0.5}
			},
		}
		e45mult2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return e.conv(r)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{0.5, 0}
			},
		}
		e45mult3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return e.conv(s)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{0, 0.5}
			},
		}
	)
	phi = []BasisPolynomialTerm{
		{PolyMultiplier: e45mult1, BasisVector: e.E4Vector},
		{PolyMultiplier: e45mult2, BasisVector: e.E4Vector},
		{PolyMultiplier: e45mult3, BasisVector: e.E4Vector},
		{PolyMultiplier: e45mult1, BasisVector: e.E5Vector},
		{PolyMultiplier: e45mult2, BasisVector: e.E5Vector},
		{PolyMultiplier: e45mult3, BasisVector: e.E5Vector},
		{PolyMultiplier: l1xiEdge1, BasisVector: e.Edge1Vector},
		{PolyMultiplier: l2xiEdge1, BasisVector: e.Edge1Vector},
		{PolyMultiplier: l3xiEdge1, BasisVector: e.Edge1Vector},
		{PolyMultiplier: l1etaEdge2, BasisVector: e.Edge2Vector},
		{PolyMultiplier: l2etaEdge2, BasisVector: e.Edge2Vector},
		{PolyMultiplier: l3etaEdge2, BasisVector: e.Edge2Vector},
		{PolyMultiplier: l3etaEdge3, BasisVector: e.Edge3Vector},
		{PolyMultiplier: l2etaEdge3, BasisVector: e.Edge3Vector},
		{PolyMultiplier: l1etaEdge3, BasisVector: e.Edge3Vector},
	}
	return
}
