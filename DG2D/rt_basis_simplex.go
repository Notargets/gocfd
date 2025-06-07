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
	bs.PEdgeBasis = NewJacobiBasis1D(P, REdge)
	bs.PKBasis = NewJacobiBasis2D(P-1, R.Copy().Subset(0, NpInt-1), S.Copy().Subset(0, NpInt-1))
	bs.ComposePhi()
	return
}

func (bs *RTBasisSimplex) ComposePhi() {
	// Constraints
	// ========= Left Edge (dot product = 0) constraint ================
	// Left normal is [-1,0], Edge functional: R = -1
	// For a vector function to vanish on the Left edge, we multiply an
	// R direction by [(R+1),]
	// ========= Left Edge (dim to (1,-1) constraint ===================
	// Multiply R direction by (R-1) to diminish approaching (1,S)
	// Multiply S direction by (S+1) to diminish approaching (R,-1)
	// [(R-1),(S+1)]
	// ========= Hypotenuse (dot product = 0) constraint ===============
	// Hypotenuse normal is [Sqrt2/2,Sqrt2/2], Edge functional: R+S = 1
	// For a vector function to vanish on the Left edge, we multiply an
	// [(1-R-S),(1-R-S)] term with the R and S vector components
	// ========= Hypotenuse (dim to (-1, -1) constraint ================
	// Multiply R direction by (R+1) to diminish approaching (-1,S)
	// Multiply S direction by (S+1) to diminish approaching (R,-1)
	// [(R+1),(S+1)]
	// ======== Bottom Edge (dot product = 0) constraint ===============
	// Bottom normal is [0,-1], Edge functional: S = -1
	// For a vector function to vanish on the Bottom edge, we multiply an
	// S direction by [,(S+1)]
	// ========== Bottom Edge (dim to (-1,1) constraint ===============
	// Multiply R direction by (R+1) to diminish approaching (-1,S)
	// Multiply S direction by (S-1) to diminish approaching (R,1)
	// [(R+1),(S-1)]
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
					// v = [ ηξ , η(η −1) ]
					// η = (R+1)/2, ξ = (S+1)/2
					// v = [ (R+1)(S+1)/4 , (R+1)((R+1)/4 -1/2) ]
					// v = [ (R+1)(S+1)/4 , (R-1)(R+1)/4 ]
					// E4 vector is tangent to all three edges:
					// E4 = [(S+1)(R+1)/4,(S+1)(S-1)/4]
					// E4' = [(S+1)(R+1),(S+1)(S-1)]
					rp := 0.5 * (r + 1.)
					sp := 0.5 * (s + 1.)
					sm := 0.5 * (s - 1.)
					// return [2]float64{sp * rp / 4., sp * sm / 4.}
					return [2]float64{sp * rp, sp * sm}
				}
				vectorDiv = func(r, s float64) (div float64) { return (3.*s + 1.) / 4. }
			case E5:
				jj = j - bs.NpInt
				vectorEval = func(r, s float64) (v [2]float64) {
					// E5 vector is tangent to all three edges:
					// E5 = [(R+1)(R-1)/4,(R+1)(S+1)/4]
					// E5' = [(R+1)(R-1),(R+1)(S+1)]
					rp := 0.5 * (r + 1.)
					rm := 0.5 * (r - 1.)
					sp := 0.5 * (s + 1.)
					// return [2]float64{rp * rm / 4., rp * sp / 4.}
					return [2]float64{rp * rm, rp * sp}
				}
				vectorDiv = func(r, s float64) (div float64) { return (3.*r + 1.) / 4. }
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
				vectorEval = func(r, s float64) (v [2]float64) {
					// Bottom edge
					// dot zero with Left and Hypotenuse yields:
					// [(1-R-S)(R+1),(1-R-S)]
					// dim toward (-1,1) yields an additional: [(R+1),(S-1)]
					// v = [(1-R-S)(R+1),(1-R-S)(S-1)]
					// h := 1. - R - S
					// v = [ ξ , η −1 ]
					// v' = [ (R+1), (S-1) ]
					rp := 0.5 * (r + 1.)
					sm := 0.5 * (s - 1.)
					return [2]float64{rp, sm}
				}
				vectorDiv = func(r, s float64) (div float64) { return 1. }
			case E2:
				sr2 := math.Sqrt2
				vectorEval = func(r, s float64) (v [2]float64) {
					// Hypotenuse
					// dot zero with Bottom and Left yields:
					// [(R+1),(S+1)]
					// dim toward (-1,-1) yields an additional: [(R+1),(S+1)]
					// [(R+1)(R+1),(S+1)(S+1)]
					// v = [√2 ξ , √2 η ]
					// v' = [√2 (R+1) , √2 (S+1) ]
					rp := 0.5 * (r + 1.)
					sp := 0.5 * (s + 1.)
					return [2]float64{sr2 * rp, sr2 * sp}
				}
				// vectorDiv = func(R, S float64) (div float64) { return 2.*(R+S) + 4. }
				vectorDiv = func(r, s float64) (div float64) { return 1. * sr2 }
			case E3:
				vectorEval = func(r, s float64) (v [2]float64) {
					// Left Edge
					// dot zero with Hypotenuse and Bottom yields:
					// [(1-R-S),(1-R-S)(S+1)]
					// dim toward (1,-1) yields an additional: [(R-1),(S+1)]
					// v = [(1-R-S)(R-1),(1-R-S)(S+1)]
					// h := 1. - R - S
					// v = [ ξ -1, η ]
					// v' = [ (R-1), (S+1) ]
					rm := 0.5 * (r - 1.)
					sp := 0.5 * (s + 1.)
					return [2]float64{rm, sp}
				}
				// vectorDiv = func(R, S float64) (div float64) { return 2.- 3.*(R+S) }
				vectorDiv = func(r, s float64) (div float64) { return 1. }
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
			// Parameterization is R
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(r, jj)
		}
		polyDeriv = func(r, s float64) (val float64) {
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(r, jj, Dr)
		}
	case E2:
		jj -= bs.NpEdge
		polyEval = func(r, s float64) (val float64) {
			// Parameterization is S
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(s, jj)
		}
		polyDeriv = func(r, s float64) (val float64) {
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(s, jj, Dr)
		}
	case E3:
		jj -= 2 * bs.NpEdge
		polyEval = func(r, s float64) (val float64) {
			// Parameterization is -S
			return bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(-s, jj)
		}
		polyDeriv = func(r, s float64) (val float64) {
			return -bs.PEdgeBasis.GetOrthogonalPolynomialAtJ(-s, jj, Dr)
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
				grad[1] = polyDeriv(r, s)
			case E3:
				grad[1] = polyDeriv(r, s)
			}
			return
		}}
	return
}
