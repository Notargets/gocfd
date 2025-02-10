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
		g1 := R.DataP[edgeLocation]
		g2 := R.DataP[edgeLocation+1]
		e.Phi = e.ComposePhiRT1(g1, g2)
	// case 2:
	// 	g1 := R.DataP[edgeLocation]
	// 	g2 := R.DataP[edgeLocation+1]
	// 	g3 := R.DataP[edgeLocation+2]
	// 	e.Phi = e.ComposePhiRT2(g1, g2, g3)
	default:
		e.InteriorPolyKBasis = NewJacobiBasis2D(e.P-1,
			R.Copy().Subset(0, e.NpInt-1),
			S.Copy().Subset(0, e.NpInt-1),
			0, 0)
		e.Phi = e.ComposePhiRTK(R.DataP[edgeLocation : edgeLocation+e.NpEdge])
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

func (e *ErvinRTBasis) getLpPolyTerm(j int, tBasis []float64) (lt BasisPolynomialTerm) {
	var (
		param = e.getFunctionType(j)
		bv    BasisVectorStruct
		jj    = j - 2*e.NpInt
		// sr2   = math.Sqrt2
	)
	switch param {
	case E1:
		bv = e.Edge1Vector
		// Multiply by 4 for Jacobian
		// bv.Divergence = func(r, s float64) (div float64) { return 1. }
	case E2:
		bv = e.Edge2Vector
		// bv.Divergence = func(r, s float64) (div float64) { return sr2 }
		jj -= e.NpEdge
	case E3:
		bv = e.Edge3Vector
		// bv.Divergence = func(r, s float64) (div float64) { return 1. }
		jj -= 2 * e.NpEdge
	default:
		panic("Lagrange polynomial is for edges only")
	}
	lagrange1D := func(r, s float64) (val float64) {
		var (
			div  = 1.
			mult = 1.
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
		tb_j := e.conv(tBasis[jj])
		for i := 0; i < e.NpEdge; i++ {
			if i != jj {
				tb_i := e.conv(tBasis[i])
				mult *= t - tb_i
				div *= tb_j - tb_i
			}
		}
		// fmt.Println("mult, div, t_rs, param, jj, r, s = ", mult, div, t_rs,
		// 	param, jj, r, s)
		return mult / div
	}
	lagrange1DDeriv := func(r, s float64) (deriv float64) {
		var (
			div  = 1.
			sum  float64
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
		tb_j := e.conv(tBasis[jj])
		for i := 0; i < e.NpEdge; i++ {
			if i != jj {
				tb_i := e.conv(tBasis[i])
				sum += t - tb_i
				div *= tb_j - tb_i
			}
		}
		return e.dconvDrDs() * sum / div
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
	lt = BasisPolynomialTerm{
		PolyMultiplier: BasisPolynomialMultiplier{
			Eval:     Eval,
			Gradient: Gradient,
		},
		BasisVector: bv,
	}
	return
}

func (e *ErvinRTBasis) getBkPolyTerm(j int) (
	bk BasisPolynomialTerm) {
	var (
		param = e.getFunctionType(j)
		jj    int
		bv    BasisVectorStruct
	)
	switch param {
	case E4:
		jj = j
		bv = e.E4Vector
		// bv.Divergence = func(r, s float64) (div float64) {
		// 	s = e.conv(s)
		// 	div = (3.*s + 1.) / 4. // Multiply by 4 for Jacobian
		// 	return
		// }
	case E5:
		jj = j - e.NpInt
		bv = e.E5Vector
		// bv.Divergence = func(r, s float64) (div float64) {
		// r = e.conv(r)
		// div = (3.*r + 1.) / 4. // Multiply by 4 for Jacobian
		// return
		// }
	default:
		panic("Bk polynomial is for interior only")
	}
	bk = BasisPolynomialTerm{
		PolyMultiplier: BasisPolynomialMultiplier{
			Eval: func(r, s float64) (val float64) {
				val = e.InteriorPolyKBasis.GetPolynomialAtJ(r, s, jj)
				return
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{
					e.InteriorPolyKBasis.GetPolynomialAtJ(r, s, jj, Dr),
					e.InteriorPolyKBasis.GetPolynomialAtJ(r, s, jj, Ds),
				}
				return
			},
		},
		BasisVector: bv,
	}
	return
}

func (e *ErvinRTBasis) ComposePhiRTK(ePts []float64) (phi []BasisPolynomialTerm) {
	// TODO: For the RTK to use the Divergence or derivatives from the E1-E5
	// TODO: vectors within the [0-1] coordinate system,
	// TODO: we need to multiply them by the Jacobian of the transform which is J=4
	// TODO: for the xi, eta = (r+1)/2, (s+1)/2 transform.
	// TODO: Just the derivatives in the [0,1] space need to be multiplied so
	//  that they are compatible with the derivatives in the [-1,1] space
	phi = make([]BasisPolynomialTerm, e.Np)
	for j := 0; j < e.Np; j++ {
		switch e.getFunctionType(j) {
		case E1, E2, E3:
			phi[j] = e.getLpPolyTerm(j, ePts)
		case E4, E5:
			phi[j] = e.getBkPolyTerm(j)
		default:
			panic("Bk polynomial wrong j")
		}
	}
	return
}

func (e *ErvinRTBasis) ComposePhiRT1(g1rs,
	g2rs float64) (phi []BasisPolynomialTerm) {
	// Compose the basis for the RT1 element based on Ervin
	var (
		g1, g2   = e.conv(g1rs), e.conv(g2rs)
		constant = BasisPolynomialMultiplier{
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

		l1xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(r) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{l1Deriv(), 0}
			},
		}
		l2xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(r) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{l2Deriv(), 0}
			},
		}

		l1etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(s) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, l1Deriv()}
			},
		}
		l2etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(s) },
			Gradient: func(r, s float64) [2]float64 {
				// return [2]float64{-l2Deriv(), l2Deriv()}
				return [2]float64{0, l2Deriv()}
			},
		}

		l1etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1Func(-s) },
			Gradient: func(r, s float64) [2]float64 {
				return [2]float64{0, -l1Deriv()}
			},
		}
		l2etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(-s) },
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

func (e *ErvinRTBasis) ComposePhiRT2(g1rs, g2rs, g3rs float64) (phi []BasisPolynomialTerm) {
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

		l1xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l1Deriv(r), 0}
				return
			},
		}
		l2xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l2Deriv(r), 0}
				return
			},
		}
		l3xiEdge1 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(r) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{l3Deriv(r), 0}
				return
			},
		}

		l1etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l1Deriv(s)}
				return
			},
		}
		l2etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l2Deriv(s)}
				return
			},
		}
		l3etaEdge2 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, l3Deriv(s)}
				return
			},
		}

		l1etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l1func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l1Deriv(-s)}
				return
			},
		}
		l2etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l2Func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l2Deriv(-s)}
				return
			},
		}
		l3etaEdge3 = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 { return l3Func(-s) },
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{0, -l3Deriv(-s)}
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
