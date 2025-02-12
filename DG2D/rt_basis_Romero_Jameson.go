package DG2D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type RomeroJamesonRTBasis struct {
	P                  int
	Np, NpInt, NpEdge  int
	Phi                []BasisPolynomialTerm
	InteriorPolyKBasis *JacobiBasis2D
}

func NewRomeroJamesonRTBasis(P int, R, S utils.Vector) (rjb *RomeroJamesonRTBasis) {
	rjb = &RomeroJamesonRTBasis{
		P:      P,
		Np:     (P + 1) * (P + 3),
		NpInt:  P * (P + 1) / 2,
		NpEdge: P + 1,
	}
	edgeLocation := 2 * rjb.NpInt
	rjb.InteriorPolyKBasis = NewJacobiBasis2D(rjb.P-1,
		R.Copy().Subset(0, rjb.NpInt-1),
		S.Copy().Subset(0, rjb.NpInt-1),
		0, 0)
	rjb.Phi = rjb.ComposePhiRTK(R.DataP[edgeLocation : edgeLocation+rjb.NpEdge])
	return
}

func (rjb *RomeroJamesonRTBasis) getFunctionType(j int) (param RTBasisFunctionType) {
	var (
		e4start    = 0
		e5start    = rjb.NpInt
		edge1Start = 2 * rjb.NpInt
		edge2Start = edge1Start + rjb.NpEdge
		edge3Start = edge2Start + rjb.NpEdge
		end        = edge3Start + rjb.NpEdge
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

func (rjb *RomeroJamesonRTBasis) getLpPolyTermUnit(j int,
	tBasis []float64) (lt BasisPolynomialTerm) {
	var (
		oosr2       = 0.5 * math.Sqrt2
		Edge1Vector = BasisVectorStruct{
			// Bottom edge
			Eval: func(r, s float64) [2]float64 { return [2]float64{0, -1} },
			Dot:  func(r, s float64, f [2]float64) float64 { return -f[1] },
			Project: func(r, s float64, scale float64) [2]float64 {
				return [2]float64{0, -scale}
			},
			Divergence: func(r, s float64) (div float64) { return 0 },
		}
		Edge2Vector = BasisVectorStruct{
			// Hypotenuse
			Eval: func(r, s float64) [2]float64 {
				return [2]float64{oosr2, oosr2}
			},
			Dot: func(r, s float64, f [2]float64) float64 {
				return oosr2*f[0] + oosr2*f[1]
			},
			Project: func(r, s float64, scale float64) [2]float64 {
				return [2]float64{scale * oosr2, scale * oosr2}
			},
			Divergence: func(r, s float64) (div float64) { return 0 },
		}
		Edge3Vector = BasisVectorStruct{
			// Left edge
			Eval: func(r, s float64) [2]float64 { return [2]float64{-1, 0} },
			Dot:  func(r, s float64, f [2]float64) float64 { return -f[0] },
			Project: func(r, s float64, scale float64) [2]float64 {
				return [2]float64{-scale}
			},
			Divergence: func(r, s float64) (div float64) { return 0 },
		}
		param = rjb.getFunctionType(j)
		bv    BasisVectorStruct
		jj    = j - 2*rjb.NpInt
	)
	switch param {
	case E1:
		bv = Edge1Vector
	case E2:
		bv = Edge2Vector
		jj -= rjb.NpEdge
	case E3:
		bv = Edge3Vector
		jj -= 2 * rjb.NpEdge
	default:
		panic("Lagrange polynomial is for edges only")
	}
	lagrange1D := func(r, s float64) (val float64) {
		var (
			div  = 1.
			mult = 1.
			t    float64
		)
		switch param {
		case E1:
			t = r
		case E2:
			t = s
		case E3:
			t = -s
		}
		tb_j := tBasis[jj]
		for i := 0; i < rjb.NpEdge; i++ {
			if i != jj {
				tb_i := tBasis[i]
				mult *= t - tb_i
				div *= tb_j - tb_i
			}
		}
		return mult / div
	}
	lagrange1DDeriv := func(r, s float64) (deriv float64) {
		var (
			div = 1.
			sum float64
			t   float64
		)
		switch param {
		case E1:
			t = r
		case E2:
			t = s
		case E3:
			t = -s
		}
		tb_j := tBasis[jj]
		for i := 0; i < rjb.NpEdge; i++ {
			if i != jj {
				tb_i := tBasis[i]
				sum += t - tb_i
				div *= tb_j - tb_i
			}
		}
		return sum / div
	}
	edgeMultiplier := func(r, s float64) (val float64) {
		var (
			t float64
		)
		switch param {
		case E1:
			t = r
		case E2:
			t = s
		case E3:
			t = -s
		}
		val = (1 - t*t)
		return
	}
	edgeMultiplierDeriv := func(r, s float64) (deriv float64) {
		var (
			t float64
		)
		switch param {
		case E1:
			t = r
		case E2:
			t = s
		case E3:
			t = -s
		}
		deriv = -2 * t
		return
	}

	Eval := func(r, s float64) (val float64) {
		val = edgeMultiplier(r, s) * lagrange1D(r, s)
		return
	}
	Gradient := func(r, s float64) (grad [2]float64) {
		switch param {
		case E1:
			grad[0] = edgeMultiplierDeriv(r, s)*lagrange1D(r, s) + lagrange1DDeriv(r, s)*edgeMultiplier(r, s)
		case E2:
			grad[1] = edgeMultiplierDeriv(r, s)*lagrange1D(r, s) + lagrange1DDeriv(r, s)*edgeMultiplier(r, s)
		case E3:
			grad[1] = edgeMultiplierDeriv(r, s)*lagrange1D(r, s) - lagrange1DDeriv(r, s)*edgeMultiplier(r, s)
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

func (rjb *RomeroJamesonRTBasis) getBkPolyTermUnit(j int) (bk BasisPolynomialTerm) {
	var (
		E4Vector = BasisVectorStruct{
			// Interior vector E4
			Eval: func(r, s float64) [2]float64 { return [2]float64{1, 0} },
			Dot:  func(r, s float64, f [2]float64) float64 { return f[0] },
			Project: func(r, s float64, scale float64) [2]float64 {
				return [2]float64{scale, 0}
			},
			Divergence: func(r, s float64) (div float64) { return 0 },
		}
		E5Vector = BasisVectorStruct{
			// Interior vector E5
			Eval: func(r, s float64) [2]float64 { return [2]float64{0, 1} },
			Dot:  func(r, s float64, f [2]float64) float64 { return f[1] },
			Project: func(r, s float64, scale float64) [2]float64 {
				return [2]float64{0, scale}
			},
			Divergence: func(r, s float64) (div float64) { return 0 },
		}
		param = rjb.getFunctionType(j)
		jj    int
		bv    BasisVectorStruct
	)
	switch param {
	case E4:
		jj = j
		bv = E4Vector
	case E5:
		jj = j - rjb.NpInt
		bv = E5Vector
	default:
		panic("Bk polynomial is for interior only")
	}
	bk = BasisPolynomialTerm{
		PolyMultiplier: BasisPolynomialMultiplier{
			Eval: func(r, s float64) (val float64) {
				val = rjb.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj)
				return
			},
			Gradient: func(r, s float64) (grad [2]float64) {
				grad = [2]float64{
					rjb.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Dr),
					rjb.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Ds),
				}
				return
			},
		},
		BasisVector: bv,
	}
	return
}

func (rjb *RomeroJamesonRTBasis) ComposePhiRTK(ePts []float64) (phi []BasisPolynomialTerm) {
	phi = make([]BasisPolynomialTerm, rjb.Np)
	for j := 0; j < rjb.Np; j++ {
		switch rjb.getFunctionType(j) {
		case E1, E2, E3:
			// phi[j] = rjb.getLpPolyTerm(j, ePts)
			phi[j] = rjb.getLpPolyTermUnit(j, ePts)
		case E4, E5:
			// phi[j] = rjb.getBkPolyTerm(j)
			phi[j] = rjb.getBkPolyTermUnit(j)
		default:
			panic("Bk polynomial wrong j")
		}
	}
	return
}
