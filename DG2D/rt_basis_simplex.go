package DG2D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type RTBasisSimplex struct {
	P                  int
	Np, NpInt, NpEdge  int
	E4Vector           BaseVector
	E5Vector           BaseVector
	Edge1Vector        BaseVector
	Edge2Vector        BaseVector
	Edge3Vector        BaseVector
	Phi                []VectorFunction
	InteriorPolyKBasis *JacobiBasis2D
}

func NewRTBasisSimplex(P int, R, S utils.Vector) (bs *RTBasisSimplex) {
	var (
		sr2   = math.Sqrt2
		Dot   = func(v1, v2 [2]float64) float64 { return v1[0]*v2[0] + v1[1]*v2[1] }
		Scale = func(v [2]float64, scale float64) [2]float64 {
			return [2]float64{v[0] * scale, v[1] * scale}
		}
		e1v = func(r, s float64) (v [2]float64) {
			// Bottom edge
			xi, eta := bs.conv(r), bs.conv(s)
			v = [2]float64{xi, eta - 1.}
			return
		}
		e2v = func(r, s float64) (v [2]float64) {
			// Hypotenuse
			xi, eta := bs.conv(r), bs.conv(s)
			v = [2]float64{sr2 * xi, sr2 * eta}
			return
		}
		e3v = func(r, s float64) (v [2]float64) {
			// Left edge
			xi, eta := bs.conv(r), bs.conv(s)
			v = [2]float64{xi - 1., eta}
			return
		}
		e4v = func(r, s float64) (v [2]float64) {
			xi, eta := bs.conv(r), bs.conv(s)
			v = [2]float64{eta * xi, eta * (eta - 1.)}
			return
		}
		e5v = func(r, s float64) (v [2]float64) {
			xi, eta := bs.conv(r), bs.conv(s)
			v = [2]float64{xi * (xi - 1.), eta * xi}
			return
		}
	)

	bs = &RTBasisSimplex{
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
				div = (3.*s + 1.) / 4.
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
				div = (3.*r + 1.) / 4.
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
				return bs.
					dconvDrDs() * 2
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
				return bs.
					dconvDrDs() * 2 * sr2
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
				return bs.
					dconvDrDs() * 2
			},
		},
	}
	edgeLocation := 2 * bs.NpInt
	bs.InteriorPolyKBasis = NewJacobiBasis2D(bs.P-1,
		R.Copy().Subset(0, bs.NpInt-1),
		S.Copy().Subset(0, bs.NpInt-1),
		0, 0)
	bs.Phi = bs.ComposePhiRTK(R.DataP[edgeLocation : edgeLocation+bs.NpEdge])
	return
}

func (bs *RTBasisSimplex) conv(r float64) float64 { return (r + 1.) / 2. }

func (bs *RTBasisSimplex) rconv(xi float64) float64 { return 2.*xi - 1. }

func (bs *RTBasisSimplex) dconvDrDs() float64 { return 1. / 2. }

func (bs *RTBasisSimplex) getFunctionType(j int) (param RTBasisFunctionType) {
	var (
		e4start    = 0
		e5start    = bs.NpInt
		edge1Start = 2 * bs.NpInt
		edge2Start = edge1Start + bs.NpEdge
		edge3Start = edge2Start + bs.NpEdge
		end        = edge3Start + bs.NpEdge
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

func (bs *RTBasisSimplex) getLpPolyTerm(j int, tBasis []float64) (lt VectorFunction) {
	var (
		param = bs.getFunctionType(j)
		bv    BaseVector
		jj    = j - 2*bs.NpInt
		// sr2   = math.Sqrt2
	)
	switch param {
	case E1:
		bv = bs.Edge1Vector
	case E2:
		bv = bs.Edge2Vector
		jj -= bs.NpEdge
	case E3:
		bv = bs.Edge3Vector
		jj -= 2 * bs.NpEdge
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
		t := bs.conv(t_rs)
		tb_j := bs.conv(tBasis[jj])
		for i := 0; i < bs.NpEdge; i++ {
			if i != jj {
				tb_i := bs.conv(tBasis[i])
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
		t := bs.conv(t_rs)
		tb_j := bs.conv(tBasis[jj])
		for i := 0; i < bs.NpEdge; i++ {
			if i != jj {
				tb_i := bs.conv(tBasis[i])
				sum += t - tb_i
				div *= tb_j - tb_i
			}
		}
		return bs.dconvDrDs() * sum / div
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

func (bs *RTBasisSimplex) getBkPolyTerm(j int) (
	bk VectorFunction) {
	var (
		param = bs.getFunctionType(j)
		jj    int
		bv    BaseVector
	)
	switch param {
	case E4:
		jj = j
		bv = bs.E4Vector
	case E5:
		jj = j - bs.NpInt
		bv = bs.E5Vector
	default:
		panic("Bk polynomial is for interior only")
	}
	bk = VectorFunction{
		PolyMultiplier: PolynomialMultiplier{
			Eval: func(r, s float64) (val float64) {
				val = bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj)
				return
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{
					bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Dr),
					bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Ds),
				}
				return
			},
		},
		VectorBase: bv,
	}
	return
}

func (bs *RTBasisSimplex) ComposePhiRTK(ePts []float64) (phi []VectorFunction) {
	phi = make([]VectorFunction, bs.Np)
	for j := 0; j < bs.Np; j++ {
		switch bs.getFunctionType(j) {
		case E1, E2, E3:
			phi[j] = bs.getLpPolyTerm(j, ePts)
		case E4, E5:
			phi[j] = bs.getBkPolyTerm(j)
		default:
			panic("Bk polynomial wrong j")
		}
	}
	return
}
