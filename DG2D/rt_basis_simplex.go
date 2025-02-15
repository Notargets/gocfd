package DG2D

import (
	"math"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type RTBasisSimplex struct {
	P                  int
	Np, NpInt, NpEdge  int
	Phi                []VectorI
	InteriorPolyKBasis *JacobiBasis2D
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
	bs.InteriorPolyKBasis = NewJacobiBasis2D(P-1,
		R.Copy().Subset(0, NpInt-1),
		S.Copy().Subset(0, NpInt-1),
		0, 0)
	// bs.Phi = bs.ComposePhiRTK(R.DataP[edgeLocation : edgeLocation+bs.NpEdge])
	return
}

func (bs *RTBasisSimplex) ComposePhi() {
	var (
		evalFunc func(r, s float64) (v [2]float64)
		divFunc  func(a, b float64) float64
		oosr2    = 0.5 * math.Sqrt2
	)
	for j := 0; j < bs.Np; j++ {
		switch bs.getFunctionType(j) {
		case E4:
			jj := j
			evalFunc = func(r, s float64) (v [2]float64) {
				v[0] = bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(
					r, s, jj)
				return
			}
			divFunc = func(r, s float64) (div float64) {
				div = bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(
					r, s, jj, Dr)
				return
			}
		case E5:
			jj := j - bs.NpInt
			evalFunc = func(r, s float64) (v [2]float64) {
				v[1] = bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(
					r, s, jj)
				return
			}
			divFunc = func(r, s float64) (div float64) {
				div = bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(
					r, s, jj, Ds)
				return
			}
		case E1:
			jj := j - 2*bs.NpInt
			evalFunc = func(r, s float64) (v [2]float64) {
				v[1] = DG1D.JacobiP(
					utils.NewVector(1, []float64{r}), 0, 0, jj)[0]
				return
			}
			divFunc = func(r, s float64) (div float64) { return }
		case E2:
			jj := j - 2*bs.NpInt + bs.NpEdge
			evalFunc = func(r, s float64) (v [2]float64) {
				v[0] = oosr2 * DG1D.JacobiP(
					utils.NewVector(1, []float64{s - r}), 0, 0, jj)[0]
				v[1] = oosr2 * DG1D.JacobiP(
					utils.NewVector(1, []float64{s - r}), 0, 0, jj)[0]
				return
			}
			divFunc = func(r, s float64) (div float64) { return }
		case E3:
			jj := j - 2*bs.NpInt + 2*bs.NpEdge
			evalFunc = func(r, s float64) (v [2]float64) {
				v[0] = -DG1D.JacobiP(
					utils.NewVector(1, []float64{s}), 0, 0, jj)[0]
				return
			}
			divFunc = func(r, s float64) (div float64) { return }
		}
		bs.Phi[j] = BaseVector{
			eval: evalFunc,
			dot: func(r, s float64, f [2]float64) (dot float64) {
				v := evalFunc(r, s)
				dot = v[0]*f[0] + v[1]*f[1]
				return
			},
			project: func(r, s, psi float64) (v [2]float64) {
				v = evalFunc(r, s)
				v[0] *= psi
				v[1] *= psi
				return
			},
			divergence: divFunc,
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

// func (bs *RTBasisSimplex) getBkPolyTerm(j int) (
// 	bk VectorFunction) {
// 	var (
// 		param = bs.getFunctionType(j)
// 		jj    int
// 		bv    BaseVector
// 	)
// 	switch param {
// 	case E4:
// 		jj = j
// 		bv = bs.E4Vector
// 	case E5:
// 		jj = j - bs.NpInt
// 		bv = bs.E5Vector
// 	default:
// 		panic("Bk polynomial is for interior only")
// 	}
// 	bk = VectorFunction{
// 		PolyMultiplier: PolynomialMultiplier{
// 			Eval: func(r, s float64) (val float64) {
// 				val = bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj)
// 				return
// 			},
// 			Gradient: func(r, s float64) (grad [2]float64) {
// 				grad = [2]float64{
// 					bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Dr),
// 					bs.InteriorPolyKBasis.GetOrthogonalPolynomialAtJ(r, s, jj, Ds),
// 				}
// 				return
// 			},
// 		},
// 		VectorBase: bv,
// 	}
// 	return
// }

// func (bs *RTBasisSimplex) ComposePhiRTK(ePts []float64) (phi []VectorFunction) {
// 	phi = make([]VectorFunction, bs.Np)
// 	for j := 0; j < bs.Np; j++ {
// 		switch bs.getFunctionType(j) {
// 		case E1, E2, E3:
// 			phi[j] = bs.getLpPolyTerm(j, ePts)
// 		case E4, E5:
// 			phi[j] = bs.getBkPolyTerm(j)
// 		default:
// 			panic("Bk polynomial wrong j")
// 		}
// 	}
// 	return
// }
