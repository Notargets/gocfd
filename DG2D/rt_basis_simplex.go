package DG2D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type RTBasisSimplex struct {
	P                 int
	Np, NpInt, NpEdge int
	Phi               []VectorI
	PKBasis           *JacobiBasis2D
}

func NewRTBasisSimplex(P int, R, S utils.Vector) (bs *RTBasisSimplex) {
	var (
		Np     = (P + 1) * (P + 3)
		NpInt  = P * (P + 1) / 2
		NpEdge = P + 1
	)
	bs = &RTBasisSimplex{
		P:      P,
		Np:     Np,
		NpInt:  NpInt,
		NpEdge: NpEdge,
		Phi:    make([]VectorI, Np),
	}
	bs.PKBasis = NewJacobiBasis2D(P-1,
		R.Copy().Subset(0, NpInt-1),
		S.Copy().Subset(0, NpInt-1),
		0, 0)
	edgeBegin := 2 * NpInt
	bs.ComposePhi(R.DataP[edgeBegin : edgeBegin+NpEdge])
	return
}

func (bs *RTBasisSimplex) ComposePhi(tBasis []float64) {
	var (
		// oosr2   = 0.5 * math.Sqrt2
		IJMap = make([][2]int, bs.NpInt)
	)

	var sk int
	for i := 0; i <= bs.P-1; i++ {
		for j := 0; j <= bs.P-1-i; j++ {
			IJMap[sk] = [2]int{i, j}
			sk++
		}
	}

	for j := 0; j < bs.Np; j++ {
		jj := j
		vectorEval := func(r, s float64) (v [2]float64) { return }
		divEval := func(r, s float64) (div float64) { return }
		evalPoly := func(r, s float64) (val float64) { return }
		evalPolyGradient := func(r, s float64) (grad [2]float64) { return }
		var polyMultiplier PolynomialMultiplier
		switch bs.getFunctionType(j) {
		case E4:
			jj = j
			vectorEval = func(r, s float64) (v [2]float64) {
				// Diminish at s boundaries
				return [2]float64{r * s, 1. - s*s}
			}
			divEval = func(r, s float64) (div float64) { return -s }
			evalPoly = func(r, s float64) float64 {
				return bs.PKBasis.GetOrthogonalPolynomialAtJ(r, s, jj)
			}
			evalPolyGradient = func(r, s float64) (grad [2]float64) {
				grad[0] = bs.PKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Dr)
				grad[1] = bs.PKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Ds)
				return
			}
			polyMultiplier = PolynomialMultiplier{
				Eval:     evalPoly,
				Gradient: evalPolyGradient,
			}
		case E5:
			jj = j - bs.NpInt
			vectorEval = func(r, s float64) (v [2]float64) {
				// Diminish at r boundaries
				return [2]float64{1. - r*r, r * s}
			}
			divEval = func(r, s float64) (div float64) { return -r }
			evalPoly = func(r, s float64) float64 {
				return bs.PKBasis.GetOrthogonalPolynomialAtJ(r, s, jj)
			}
			evalPolyGradient = func(r, s float64) (grad [2]float64) {
				grad[0] = bs.PKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Dr)
				grad[1] = bs.PKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Ds)
				return
			}
			polyMultiplier = PolynomialMultiplier{
				Eval:     evalPoly,
				Gradient: evalPolyGradient,
			}
		case E1:
			jj = j - 2*bs.NpInt
			vectorEval = func(r, s float64) (v [2]float64) {
				// Diminish away from bottom
				return [2]float64{r, (1. - s) / 2.}
			}
			divEval = func(r, s float64) (div float64) { return 1 / 2. }
			polyMultiplier = bs.getLpPolyTerm(j, tBasis)
		case E2:
			jj = j - 2*bs.NpInt - bs.NpEdge
			sr2t2 := 2. * math.Sqrt2
			vectorEval = func(r, s float64) (v [2]float64) {
				// Diminish away from Hypotenuse
				// 1 = r+s, (r+s+2)/3 is the contravariant
				return [2]float64{(r + s + 2.) / sr2t2, (r + s + 2.) / sr2t2}
			}
			divEval = func(r, s float64) (div float64) { return 2. / sr2t2 }
			polyMultiplier = bs.getLpPolyTerm(j, tBasis)
		case E3:
			jj = j - 2*bs.NpInt - 2*bs.NpEdge
			vectorEval = func(r, s float64) (v [2]float64) {
				// Diminish away from Left
				return [2]float64{(1. - r) / 2., s}
			}
			divEval = func(r, s float64) (div float64) { return 1 / 2. }
			polyMultiplier = bs.getLpPolyTerm(j, tBasis)
		default:
			panic("invalid function type")
		}

		bs.Phi[j] = VectorFunction{
			PolyMultiplier: polyMultiplier,
			VectorBase: BaseVector{
				eval: func(r, s float64) [2]float64 { return vectorEval(r, s) },
				dot: func(r, s float64, f [2]float64) (dot float64) {
					v := vectorEval(r, s)
					dot = v[0]*f[0] + v[1]*f[1]
					return
				},
				project: func(r, s, psi float64) (v [2]float64) {
					v = vectorEval(r, s)
					v[0] *= psi
					v[1] *= psi
					return
				},
				divergence: func(r, s float64) float64 { return divEval(r, s) },
			},
		}
	}
}

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

func (bs *RTBasisSimplex) getLpPolyTerm(j int, tBasis []float64) (pm PolynomialMultiplier) {
	var (
		param = bs.getFunctionType(j)
		sr2   = math.Sqrt2
	)
	jj := j - 2*bs.NpInt
	switch param {
	case E1:
		// Do nothing
	case E2:
		jj -= bs.NpEdge
	case E3:
		jj -= 2 * bs.NpEdge
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
		for i := 0; i < bs.NpEdge; i++ {
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
			div    = 1.
			sum, t float64
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
		for i := 0; i < bs.NpEdge; i++ {
			if i != jj {
				tb_i := tBasis[i]
				sum += t - tb_i
				div *= tb_j - tb_i
			}
		}
		return sum / div
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
			grad[1] = sr2 * lagrange1DDeriv(r, s)
		case E3:
			grad[1] = -lagrange1DDeriv(r, s)
		}
		return
	}
	pm = PolynomialMultiplier{
		Eval:     Eval,
		Gradient: Gradient,
	}
	return
}

// JacobiP computes the Jacobi polynomial of given order, alpha, and beta.
// Uses the standard three-term recurrence relation.
// JacobiP = func(coord float64, alpha float64, beta float64, order int) float64 {
// 	// Base cases
// 	if order == 0 {
// 		return 1.0
// 	}
// 	if order == 1 {
// 		return 0.5 * ((alpha - beta) + (alpha+beta+2.0)*coord)
// 	}
//
// 	// Recurrence relation
// 	var p0, p1 float64 = 1.0, 0.5 * ((alpha - beta) + (alpha+beta+2.0)*coord)
// 	var p2 float64
//
// 	for n := 1; n < order; n++ {
// 		a1 := 2.0 * float64(n) * (float64(n) + alpha + beta) * (2.0*float64(n) + alpha + beta - 2.0)
// 		a2 := (2.0*float64(n) + alpha + beta - 1.0) * ((alpha * alpha) - (beta * beta))
// 		a3 := (2.0*float64(n) + alpha + beta - 2.0) * (2.0*float64(n) + alpha + beta)
// 		a4 := 2.0 * (float64(n) + alpha) * (float64(n) + beta) * (2.0*float64(n) + alpha + beta)
//
// 		p2 = ((a2+a3*coord)*p1 - a4*p0) / a1
//
// 		// Shift for next iteration
// 		p0, p1 = p1, p2
// 	}
//
// 	return p2
// }
//
// // GradJacobiP computes the derivative of the Jacobi polynomial using a known identity.
// GradJacobiP = func(coord float64, alpha float64, beta float64, order int) float64 {
// 	if order == 0 {
// 		return 0.0
// 	}
// 	// Use the known relationship for Jacobi polynomial derivatives:
// 	// d/dx P_n^{(α,β)}(x) = 0.5 * (α + β + n + 1) * P_{n-1}^{(α+1,β+1)}(x)
// 	return 0.5 * (alpha + beta + float64(order) + 1.0) * JacobiP(coord, alpha+1, beta+1, order-1)
// }
