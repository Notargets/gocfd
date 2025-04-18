package DG2D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type ErvinRTBasis struct {
	P                  int
	Np, NpInt, NpEdge  int
	E4Vector           BaseVector
	E5Vector           BaseVector
	Edge1Vector        BaseVector
	Edge2Vector        BaseVector
	Edge3Vector        BaseVector
	Phi                []VectorI
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
		E4Vector: BaseVector{
			// Interior vector E4
			eval: func(r, s float64) [2]float64 { return e4v(r, s) },
			dot: func(r, s float64, f [2]float64) float64 {
				return Dot(e4v(r, s), f)
			},
			project: func(r, s float64, scale float64) [2]float64 {
				return Scale(e4v(r, s), scale)
			},
			divergence: func(r, s float64) (div float64) {
				// eta + eta -1 + eta = 3 eta -1 = 3*(S+1)/2 -1 = (3s+3)/2 -2/2
				// = (3s+1)/2
				div = e.dconvDrDs() * (3.*s + 1.) / 2.
				return
			},
		},
		E5Vector: BaseVector{
			// Interior vector E5
			eval: func(r, s float64) [2]float64 { return e5v(r, s) },
			dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e5v(r, s), f)
			},
			project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e5v(r, s), scale)
			},
			divergence: func(r, s float64) (div float64) {
				div = e.dconvDrDs() * (3.*r + 1.) / 2.
				return
			},
		},
		Edge1Vector: BaseVector{
			// Bottom edge
			eval: func(r, s float64) (v [2]float64) { return e1v(r, s) },
			dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e1v(r, s), f)
			},
			project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e1v(r, s), scale)
			},
			divergence: func(r, s float64) (div float64) {
				return e.dconvDrDs() * 2
			},
		},
		Edge2Vector: BaseVector{
			// Hypotenuse
			eval: func(r, s float64) (v [2]float64) { return e2v(r, s) },
			dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e2v(r, s), f)
			},
			project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e2v(r, s), scale)
			},
			divergence: func(r, s float64) (div float64) {
				return e.dconvDrDs() * 2 * sr2
			},
		},
		Edge3Vector: BaseVector{
			// Left edge
			eval: func(r, s float64) (v [2]float64) { return e3v(r, s) },
			dot: func(r, s float64, f [2]float64) (dot float64) {
				return Dot(e3v(r, s), f)
			},
			project: func(r, s float64, scale float64) (v [2]float64) {
				return Scale(e3v(r, s), scale)
			},
			divergence: func(r, s float64) (div float64) {
				return e.
					dconvDrDs() * 2
			},
		},
	}
	edgeLocation := 2 * e.NpInt
	switch e.P {
	case 1:
		g1 := R.DataP[edgeLocation]
		g2 := R.DataP[edgeLocation+1]
		e.Phi = e.ComposePhiRT1(g1, g2)
	case 2:
		g1 := R.DataP[edgeLocation]
		g2 := R.DataP[edgeLocation+1]
		g3 := R.DataP[edgeLocation+2]
		e.Phi = e.ComposePhiRT2(g1, g2, g3)
	default:
		panic("Ervin RT basis is not working with P>2")
		// e.InteriorPolyKBasis = NewJacobiBasis2D(e.P-1,
		// 	R.Copy().Subset(0, e.NpInt-1),
		// 	S.Copy().Subset(0, e.NpInt-1),
		// 	0, 0)
		// e.Phi = e.ComposePhiRTK(R.DataP[edgeLocation : edgeLocation+e.NpEdge])
	}
	return
}

func (e *ErvinRTBasis) conv(r float64) float64 { return (r + 1.) / 2. }

func (e *ErvinRTBasis) rconv(xi float64) float64 { return 2.*xi - 1. }

func (e *ErvinRTBasis) dconvDrDs() float64 { return 1. / 2. }

func (e *ErvinRTBasis) getFunctionType(j int) (param RTBasisFunctionType) {
	var (
		e4start    = 0
		e5start    = e.NpInt
		edge1Start = 2 * e.NpInt
		edge2Start = edge1Start + e.NpEdge
		edge3Start = edge2Start + e.NpEdge
		end        = edge3Start + e.NpEdge
	)
	switch {
	case j >= e4start && j < e5start:
		param = E4
	case j >= e5start && j < edge1Start:
		param = E5
	case j >= edge1Start && j < edge2Start:
		param = E1
	case j >= edge2Start && j < edge3Start:
		param = E2
	case j >= edge3Start && j < end:
		param = E3
	case j >= end:
		panic("j out of range")
	}
	return
}

func (e *ErvinRTBasis) getLpPolyTerm(j int, tBasis []float64) (lt VectorFunction) {
	var (
		param = e.getFunctionType(j)
		// sr2   = math.Sqrt2
	)
	jj := j - 2*e.NpInt
	bv := BaseVector{}
	switch param {
	case E1:
		bv = e.Edge1Vector
	case E2:
		bv = e.Edge2Vector
		jj -= e.NpEdge
	case E3:
		bv = e.Edge3Vector
		jj -= 2 * e.NpEdge
	default:
		panic("Lagrange polynomial is for edges only")
	}
	lagrange1D := func(r, s float64) (val float64) {
		var (
			t_rs float64
		)
		switch param {
		case E1:
			t_rs = r
		case E2:
			t_rs = s
		case E3:
			t_rs = -s
		}
		t := e.conv(t_rs)
		val = Lagrange1DPoly(t, tBasis, jj)
		return
	}

	lagrange1DDeriv := func(r, s float64) (deriv float64) {
		var (
			t_rs float64
		)
		switch param {
		case E1:
			t_rs = r
		case E2:
			t_rs = s
		case E3:
			t_rs = -s
		}
		t := e.conv(t_rs)
		deriv = e.dconvDrDs() * Lagrange1DPoly(t, tBasis, jj, Dr)
		return
	}

	Eval := func(r, s float64) (val float64) {
		val = lagrange1D(r, s)
		return
	}
	Gradient := func(r, s float64) (grad [2]float64) {
		switch param {
		case E1:
			grad[0] = lagrange1DDeriv(r, s)
		case E2:
			grad[1] = lagrange1DDeriv(r, s)
		case E3:
			grad[1] = -lagrange1DDeriv(r, s)
		}
		return
	}
	lt = VectorFunction{
		PolyMultiplier: PolynomialMultiplier{
			Eval:     Eval,
			Gradient: Gradient,
		},
		VectorBase: bv,
	}
	return
}

func (e *ErvinRTBasis) getBkPolyTerm(j int) (
	bk VectorFunction) {
	var (
		param = e.getFunctionType(j)
		jj    int
		bv    BaseVector
	)
	switch param {
	case E4:
		jj = j
		bv = e.E4Vector
	case E5:
		jj = j - e.NpInt
		bv = e.E5Vector
	default:
		panic("Bk polynomial is for interior only")
	}
	bk = VectorFunction{
		PolyMultiplier: PolynomialMultiplier{
			Eval: func(r, s float64) (val float64) {
				val = e.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj)
				val = e.conv(val)
				return
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				// grad = [2]float64{
				// 	e.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(R, S, jj, Dr),
				// 	e.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(R, S, jj, Ds),
				// }
				d1 := e.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Dr)
				d2 := e.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Ds)
				grad = [2]float64{d1, d2}
				return
			},
		},
		VectorBase: bv,
	}
	return
}

func (e *ErvinRTBasis) ComposePhiRTK(ePts []float64) (phi []VectorI) {
	phi = make([]VectorI, e.Np)
	ePtsConv := make([]float64, len(ePts))
	for j := range ePts {
		ePtsConv[j] = e.conv(ePts[j])
	}
	for j := 0; j < e.Np; j++ {
		switch e.getFunctionType(j) {
		case E1, E2, E3:
			phi[j] = e.getLpPolyTerm(j, ePtsConv)
		case E4, E5:
			phi[j] = e.getBkPolyTerm(j)
		default:
			panic("Bk polynomial wrong j")
		}
	}
	return
}

func (e *ErvinRTBasis) ComposePhiRT1(g1rs, g2rs float64) (phi []VectorI) {
	// Compose the basis for the RT1 element based on Ervin
	var (
		g1, g2   = e.conv(g1rs), e.conv(g2rs)
		constant = PolynomialMultiplier{
			Eval:     func(r, s float64) float64 { return 1. },
			Gradient: func(r, s float64) (grad [2]float64) { return },
		}

		l1Func = func(t_rs float64) float64 {
			return (e.conv(t_rs) - g2) / (g1 - g2)
		}
		l1Deriv = func() float64 {
			return e.dconvDrDs() / (g1 - g2)
		}

		l2Func = func(t_rs float64) float64 {
			return (e.conv(t_rs) - g1) / (g2 - g1)
		}
		l2Deriv = func() float64 { return -l1Deriv() }

		l1xiEdge1 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(r) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{l1Deriv(), 0}
			},
		}
		l2xiEdge1 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(r) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{l2Deriv(), 0}
			},
		}

		l1etaEdge2 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(s) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, l1Deriv()}
			},
		}
		l2etaEdge2 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(s) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, l2Deriv()}
			},
		}

		l1etaEdge3 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(s) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, l1Deriv()}
			},
		}
		l2etaEdge3 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(s) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, l2Deriv()}
			},
		}
	)
	phiVF := []VectorFunction{
		{PolyMultiplier: constant, VectorBase: e.E4Vector},
		{PolyMultiplier: constant, VectorBase: e.E5Vector},
		{PolyMultiplier: l1xiEdge1, VectorBase: e.Edge1Vector},
		{PolyMultiplier: l2xiEdge1, VectorBase: e.Edge1Vector},
		{PolyMultiplier: l1etaEdge2, VectorBase: e.Edge2Vector},
		{PolyMultiplier: l2etaEdge2, VectorBase: e.Edge2Vector},
		{PolyMultiplier: l2etaEdge3, VectorBase: e.Edge3Vector},
		{PolyMultiplier: l1etaEdge3, VectorBase: e.Edge3Vector},
	}
	phi = make([]VectorI, len(phiVF))
	for i := range phi {
		phi[i] = phiVF[i]
	}
	return
}

func (e *ErvinRTBasis) ComposePhiRT2(g1rs, g2rs, g3rs float64) (phi []VectorI) {
	// Compose the basis for the RT1 element based on Ervin
	var (
		g1, g2, g3 = e.conv(g1rs), e.conv(g2rs), e.conv(g3rs)
		// Get the two edge points for use in the lagrange terms for RT1
		div1 = 1. / ((g1 - g2) * (g1 - g3))
		div2 = 1. / ((g2 - g1) * (g2 - g3))
		div3 = 1. / ((g3 - g1) * (g3 - g2))

		l1func = func(t_rs float64) float64 {
			return (e.conv(t_rs) - g2) * (e.conv(t_rs) - g3) * div1
		}
		l2Func = func(t_rs float64) float64 {
			return (e.conv(t_rs) - g1) * (e.conv(t_rs) - g3) * div2
		}
		l3Func = func(t_rs float64) float64 {
			return (e.conv(t_rs) - g1) * (e.conv(t_rs) - g2) * div3
		}

		l1Deriv = func(t_rs float64) (deriv float64) {
			deriv = ((e.conv(t_rs) - g2) + (e.conv(t_rs) - g3)) * div1
			deriv *= e.dconvDrDs()
			return
		}
		l2Deriv = func(t_rs float64) (deriv float64) {
			deriv = ((e.conv(t_rs) - g1) + (e.conv(t_rs) - g3)) * div2
			deriv *= e.dconvDrDs()
			return
		}
		l3Deriv = func(t_rs float64) (deriv float64) {
			deriv = ((e.conv(t_rs) - g1) + (e.conv(t_rs) - g2)) * div3
			deriv *= e.dconvDrDs()
			return
		}

		l1xiEdge1 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l1Deriv(r), 0}
				return
			},
		}
		l2xiEdge1 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l2Deriv(r), 0}
				return
			},
		}
		l3xiEdge1 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l3Deriv(r), 0}
				return
			},
		}

		l1etaEdge2 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l1Deriv(s)}
				return
			},
		}
		l2etaEdge2 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l2Deriv(s)}
				return
			},
		}
		l3etaEdge2 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l3Deriv(s)}
				return
			},
		}

		l1etaEdge3 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l1Deriv(-s)}
				return
			},
		}
		l2etaEdge3 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l2Deriv(-s)}
				return
			},
		}
		l3etaEdge3 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l3Deriv(-s)}
				return
			},
		}

		e45mult1 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return 1. - e.conv(r) - e.conv(s)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{-0.5, -0.5}
			},
		}
		e45mult2 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return e.conv(r)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{0.5, 0}
			},
		}
		e45mult3 = PolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return e.conv(s)
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				return [2]float64{0, 0.5}
			},
		}
	)
	phiVF := []VectorFunction{
		{PolyMultiplier: e45mult1, VectorBase: e.E4Vector},
		{PolyMultiplier: e45mult2, VectorBase: e.E4Vector},
		{PolyMultiplier: e45mult3, VectorBase: e.E4Vector},
		{PolyMultiplier: e45mult1, VectorBase: e.E5Vector},
		{PolyMultiplier: e45mult2, VectorBase: e.E5Vector},
		{PolyMultiplier: e45mult3, VectorBase: e.E5Vector},
		{PolyMultiplier: l1xiEdge1, VectorBase: e.Edge1Vector},
		{PolyMultiplier: l2xiEdge1, VectorBase: e.Edge1Vector},
		{PolyMultiplier: l3xiEdge1, VectorBase: e.Edge1Vector},
		{PolyMultiplier: l1etaEdge2, VectorBase: e.Edge2Vector},
		{PolyMultiplier: l2etaEdge2, VectorBase: e.Edge2Vector},
		{PolyMultiplier: l3etaEdge2, VectorBase: e.Edge2Vector},
		{PolyMultiplier: l3etaEdge3, VectorBase: e.Edge3Vector},
		{PolyMultiplier: l2etaEdge3, VectorBase: e.Edge3Vector},
		{PolyMultiplier: l1etaEdge3, VectorBase: e.Edge3Vector},
	}
	phi = make([]VectorI, len(phiVF))
	for i := range phi {
		phi[i] = phiVF[i]
	}
	return
}
