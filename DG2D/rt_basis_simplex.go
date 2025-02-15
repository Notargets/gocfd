package DG2D

import (
	"math"

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
		IJMap    = make([][2]int, bs.Np)
		sk       int
		// JacobiP computes the Jacobi polynomial of given order, alpha, and beta.
		// Uses the standard three-term recurrence relation.
		JacobiP = func(coord float64, alpha float64, beta float64, order int) float64 {
			// Base cases
			if order == 0 {
				return 1.0
			}
			if order == 1 {
				return 0.5 * ((alpha - beta) + (alpha+beta+2.0)*coord)
			}

			// Recurrence relation
			var p0, p1 float64 = 1.0, 0.5 * ((alpha - beta) + (alpha+beta+2.0)*coord)
			var p2 float64

			for n := 1; n < order; n++ {
				a1 := 2.0 * float64(n) * (float64(n) + alpha + beta) * (2.0*float64(n) + alpha + beta - 2.0)
				a2 := (2.0*float64(n) + alpha + beta - 1.0) * ((alpha * alpha) - (beta * beta))
				a3 := (2.0*float64(n) + alpha + beta - 2.0) * (2.0*float64(n) + alpha + beta)
				a4 := 2.0 * (float64(n) + alpha) * (float64(n) + beta) * (2.0*float64(n) + alpha + beta)

				p2 = ((a2+a3*coord)*p1 - a4*p0) / a1

				// Shift for next iteration
				p0, p1 = p1, p2
			}

			return p2
		}
		// GradJacobiP computes the derivative of the Jacobi polynomial using a known identity.
		GradJacobiP = func(coord float64, alpha float64, beta float64, order int) float64 {
			if order == 0 {
				return 0.0
			}
			// Use the known relationship for Jacobi polynomial derivatives:
			// d/dx P_n^{(α,β)}(x) = 0.5 * (α + β + n + 1) * P_{n-1}^{(α+1,β+1)}(x)
			return 0.5 * (alpha + beta + float64(order) + 1.0) * JacobiP(coord, alpha+1, beta+1, order-1)
		}
	)

	for i := 0; i <= bs.P-1; i++ {
		for j := 0; j <= bs.P-1-i; j++ {
			IJMap[sk] = [2]int{i, j}
			sk++
		}
	}
	for j := 0; j < bs.Np; j++ {
		JP2D_j := func(j int) func(r, s float64) (val float64) {
			ii, jj := IJMap[j][0], IJMap[j][1]
			return func(r, s float64) (val float64) {
				val1 := JacobiP(r, 1, 0, ii)
				val2 := JacobiP(0.5*(s+1), 0, 1, jj)
				val = val1 * val2
				return
			}
		}
		JP2DDerivR_j := func(j int) func(r, s float64) (val float64) {
			ii, jj := IJMap[j][0], IJMap[j][1]
			return func(r, s float64) (df float64) {
				val2 := JacobiP(0.5*(s+1), 0, 1, jj)
				df1 := GradJacobiP(r, 1, 0, ii)
				df = val2 * df1
				return
			}
		}
		JP2DDerivS_j := func(j int) func(r, s float64) (val float64) {
			ii, jj := IJMap[j][0], IJMap[j][1]
			return func(r, s float64) (df float64) {
				val1 := JacobiP(r, 1, 0, ii)
				df2 := GradJacobiP(0.5*(s+1), 0, 1, jj)
				df = 0.5 * val1 * df2
				return
			}
		}
		switch bs.getFunctionType(j) {
		case E4:
			jj := j
			evalFunc = func(r, s float64) (v [2]float64) {
				v[0] = JP2D_j(jj)(r, s)
				return
			}
			divFunc = func(r, s float64) (div float64) {
				div = JP2DDerivR_j(jj)(r, s)
				return
			}
		case E5:
			jj := j - bs.NpInt
			evalFunc = func(r, s float64) (v [2]float64) {
				v[1] = JP2D_j(jj)(r, s)
				return
			}
			divFunc = func(r, s float64) (div float64) {
				div = JP2DDerivS_j(jj)(r, s)
				return
			}
		case E1:
			jj := j - 2*bs.NpInt
			evalFunc = func(r, s float64) (v [2]float64) {
				v[1] = r * JacobiP(r, 0, 0, jj) // Multiply by r
				return
			}
			divFunc = func(r, s float64) (div float64) {
				return JacobiP(r, 0, 0, jj) + r*GradJacobiP(r, 0, 0, jj)
			}
		case E2:
			jj := j - 2*bs.NpInt - bs.NpEdge
			evalFunc = func(r, s float64) (v [2]float64) {
				val := oosr2 * JacobiP(s-r, 0, 0, jj)
				v[0] = (s - r) * val // Multiply by (s-r)
				v[1] = (s - r) * val
				return
			}
			divFunc = func(r, s float64) (div float64) {
				val := oosr2 * JacobiP(s-r, 0, 0, jj)
				df := oosr2 * GradJacobiP(s-r, 0, 0, jj)
				return (s-r)*df + 2*val
			}
		case E3:
			jj := j - 2*bs.NpInt - 2*bs.NpEdge
			evalFunc = func(r, s float64) (v [2]float64) {
				v[0] = -s * JacobiP(s, 0, 0, jj) // Multiply by s
				return
			}
			divFunc = func(r, s float64) (div float64) {
				return -JacobiP(s, 0, 0, jj) - s*GradJacobiP(s, 0, 0, jj)
			}
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
