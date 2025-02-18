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
	PEdgeBasis        *JacobiBasis1D
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
	edge1Begin := 2 * NpInt
	edge2Begin := edge1Begin + NpEdge
	REdge := utils.NewVector(NpEdge, R.DataP[edge1Begin:edge2Begin])
	bs.PEdgeBasis = NewJacobiBasis1D(P, REdge, 0, 0)
	bs.PKBasis = NewJacobiBasis2D(P-1,
		R.Copy().Subset(0, NpInt-1),
		S.Copy().Subset(0, NpInt-1),
		0, 0)
	// bs.ComposePhi(R.DataP[edgeBegin : edgeBegin+NpEdge])
	bs.ComposePhi()
	return
}

// func (bs *RTBasisSimplex) ComposePhi(tBasis []float64) {
func (bs *RTBasisSimplex) ComposePhi() {
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

	// Constraints
	// ========= Left Edge (dot product = 0) constraint ================
	// Left normal is [-1,0], Edge functional: r = -1
	// For a vector function to vanish on the Left edge, we multiply an
	// r direction by [(r+1),]
	// ========= Left Edge (dim to (1,-1) constraint ===================
	// Multiply r direction by (r-1) to diminish approaching (1,s)
	// Multiply s direction by (s+1) to diminish approaching (r,-1)
	// [(r-1),(s+1)]
	// ========= Hypotenuse (dot product = 0) constraint ===============
	// Hypotenuse normal is [Sqrt2/2,Sqrt2/2], Edge functional: r+s = 1
	// For a vector function to vanish on the Left edge, we multiply an
	// [(1-r-s),(1-r-s)] term with the r and s vector components
	// ========= Hypotenuse (dim to (-1, -1) constraint ================
	// Multiply r direction by (r+1) to diminish approaching (-1,s)
	// Multiply s direction by (s+1) to diminish approaching (r,-1)
	// [(r+1),(s+1)]
	// ======== Bottom Edge (dot product = 0) constraint ===============
	// Bottom normal is [0,-1], Edge functional: s = -1
	// For a vector function to vanish on the Bottom edge, we multiply an
	// s direction by [,(s+1)]
	// ========== Bottom Edge (dim to (-1,1) constraint ===============
	// Multiply r direction by (r+1) to diminish approaching (-1,s)
	// Multiply s direction by (s-1) to diminish approaching (r,1)
	// [(r+1),(s-1)]
	for j := 0; j < bs.Np; j++ {
		jj := j
		vectorEval := func(r, s float64) (v [2]float64) { return }
		vectorDiv := func(r, s float64) (div float64) { return }
		evalPoly := func(r, s float64) (val float64) { return }
		evalPolyGradient := func(r, s float64) (grad [2]float64) { return }
		var polyMultiplier PolynomialMultiplier
		ftype := bs.getFunctionType(j)
		switch ftype {
		case E4, E5:
			switch ftype {
			case E4:
				jj = j
				vectorEval = func(r, s float64) (v [2]float64) {
					// E4 vector is tangent to all three edges:
					// E4 = [(s+1)(r+1)/4,(s+1)(s-1)/4]
					rp := r + 1.
					sp := s + 1.
					sm := s - 1.
					return [2]float64{sp * rp / 4., sp * sm / 4.}
				}
				vectorDiv = func(r, s float64) (div float64) { return (3.*s + 1) / 4. }
			case E5:
				jj = j - bs.NpInt
				vectorEval = func(r, s float64) (v [2]float64) {
					// E5 vector is tangent to all three edges:
					// E5 = [(r+1)(r-1)/4,(r+1)(s+1)/4]
					rp := r + 1.
					rm := r - 1.
					sp := s + 1.
					return [2]float64{rp * rm / 4., rp * sp / 4.}
				}
				vectorDiv = func(r, s float64) (div float64) { return (3.*r + 1) / 4. }
			}
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
		case E1, E2, E3:
			// polyMultiplier = bs.getLpPolyTerm(j, tBasis)
			polyMultiplier = bs.getEdgePolyTerm(j)
			switch ftype {
			case E1:
				jj = j - 2*bs.NpInt
				vectorEval = func(r, s float64) (v [2]float64) {
					// Bottom edge
					// dot zero with Left and Hypotenuse yields:
					// [(1-r-s)(r+1),(1-r-s)]
					// dim toward (-1,1) yields an additional: [(r+1),(s-1)]
					// v = [(1-r-s)(r+1),(1-r-s)(s-1)]
					h := 1. - r - s
					rp := r + 1.
					sm := s - 1.
					return [2]float64{h * rp, h * sm}
				}
				vectorDiv = func(r, s float64) (div float64) { return 2. - 3.*(r+s) }
			case E2:
				jj = j - 2*bs.NpInt - bs.NpEdge
				vectorEval = func(r, s float64) (v [2]float64) {
					// Hypotenuse
					// dot zero with Bottom and Left yields:
					// [(r+1),(s+1)]
					// dim toward (-1,-1) yields an additional: [(r+1),(s+1)]
					// [(r+1)(r+1),(s+1)(s+1)]
					rp := r + 1.
					sp := s + 1.
					return [2]float64{rp * rp, sp * sp}
				}
				vectorDiv = func(r, s float64) (div float64) { return 2.*(r+s) + 4. }
			case E3:
				jj = j - 2*bs.NpInt - 2*bs.NpEdge
				vectorEval = func(r, s float64) (v [2]float64) {
					// Left Edge
					// dot zero with Hypotenuse and Bottom yields:
					// [(1-r-s),(1-r-s)(s+1)]
					// dim toward (1,-1) yields an additional: [(r-1),(s+1)]
					// v = [(1-r-s)(r-1),(1-r-s)(s+1)]
					h := 1. - r - s
					rm := r - 1.
					sp := s + 1.
					return [2]float64{h * rm, h * sp}
				}
				vectorDiv = func(r, s float64) (div float64) { return 2. - 3.*(r+s) }
			default:
				panic("invalid function type")
			}
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
				divergence: func(r, s float64) float64 { return vectorDiv(r, s) },
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

func (bs *RTBasisSimplex) getEdgePolyTerm(j int) (pm PolynomialMultiplier) {
	var (
		param = bs.getFunctionType(j)
	)
	polyEval := func(r, s float64) (val float64) { return }
	polyDeriv := func(r, s float64) (val float64) { return }
	jj := j - 2*bs.NpInt
	switch param {
	case E1:
		polyEval = func(r, s float64) (val float64) {
			// Parameterization is r
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(r, jj)
		}
		polyDeriv = func(r, s float64) (val float64) {
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(r, jj, Dr)
		}
	case E2:
		jj -= bs.NpEdge
		polyEval = func(r, s float64) (val float64) {
			// Parameterization is s
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(s, jj)
		}
		polyDeriv = func(r, s float64) (val float64) {
			// Parameterization is (s+1)/2
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(s, jj, Dr)
		}
	case E3:
		jj -= 2 * bs.NpEdge
		polyEval = func(r, s float64) (val float64) {
			// Parameterization is -s
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(-s, jj)
		}
		polyDeriv = func(r, s float64) (val float64) {
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(-s, jj, Dr)
		}
	default:
		panic("Lagrange polynomial is for edges only")
	}
	pm = PolynomialMultiplier{
		Eval: func(r, s float64) float64 { return polyEval(r, s) },
		Gradient: func(r, s float64) (grad [2]float64) {
			switch param {
			case E1:
				grad[0] = polyDeriv(r, s)
			case E2:
				grad[1] = polyDeriv(r, s) / math.Sqrt2
			case E3:
				grad[1] = polyDeriv(r, s)
			}
			return
		}}
	return
}

func (bs *RTBasisSimplex) getLpPolyTerm(j int, tBasis []float64) (pm PolynomialMultiplier) {
	var (
		param = bs.getFunctionType(j)
	)
	jj := j - 2*bs.NpInt
	tB := []float64{}
	switch param {
	case E1:
		tB = tBasis
	case E2:
		jj -= bs.NpEdge
		tB = make([]float64, bs.NpEdge)
		for i := 0; i < bs.NpEdge; i++ {
			tB[i] = (tBasis[i] + 1.) / 2. // t = (s+1)/2
			// fmt.Printf("tB[%d,%f] = %f\n", i, tBasis[i], tB[i])
		}
	case E3:
		jj -= 2 * bs.NpEdge
		// tB = tBasis
		tB = make([]float64, bs.NpEdge)
		for i := 0; i < bs.NpEdge; i++ {
			tB[bs.NpEdge-1-i] = -tBasis[i] // t = -s
		}
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
			t = (s + 1.) / 2.
		case E3:
			t = -s
		}
		tb_j := tB[jj]
		for i := 0; i < bs.NpEdge; i++ {
			if i != jj {
				tb_i := tB[i]
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
			t = (s + 1.) / 2.
		case E3:
			t = -s
		}
		tb_j := tB[jj]
		for i := 0; i < bs.NpEdge; i++ {
			if i != jj {
				tb_i := tB[i]
				sum += t - tb_i
				div *= tb_j - tb_i
			}
		}
		return sum / div
	}

	pm = PolynomialMultiplier{
		Eval: func(r, s float64) float64 { return lagrange1D(r, s) },
		Gradient: func(r, s float64) (grad [2]float64) {
			switch param {
			case E1:
				grad[0] = lagrange1DDeriv(r, s)
			case E2:
				grad[1] = lagrange1DDeriv(r, s) / 2. // Parameterization is (s+1)/2
			case E3:
				grad[1] = -lagrange1DDeriv(r, s)
			}
			return
		}}
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
