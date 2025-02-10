package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type JacobiBasis2D struct {
	P               int // Order
	Np              int // Dimension
	Alpha, Beta     float64
	V, Vinv, Vr, Vs utils.Matrix
}

func NewJacobiBasis2D(P int, R, S utils.Vector, Alpha, Beta float64) (jb2d *JacobiBasis2D) {
	jb2d = &JacobiBasis2D{
		P:     P,
		Np:    (P + 1) * (P + 2) / 2,
		Alpha: Alpha,
		Beta:  Beta,
	}
	if R.Len() != jb2d.Np {
		fmt.Printf("Length of R:%d, required length:%d\n", R.Len(), jb2d.Np)
		panic("Mismatch of length for basis coordinates")
	}
	jb2d.V = jb2d.Vandermonde2D(P, R, S)
	jb2d.Vinv = jb2d.V.InverseWithCheck()
	jb2d.Vr, jb2d.Vs = jb2d.GradVandermonde2D(P, R, S)
	return
}

func (jb2d *JacobiBasis2D) Vandermonde2D(N int, R, S utils.Vector) (V2D utils.Matrix) {
	V2D = utils.NewMatrix(R.Len(), jb2d.Np)
	var sk int
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			V2D.SetCol(sk, jb2d.Simplex2DP(R, S, i, j))
			sk++
		}
	}
	return
}

func (jb2d *JacobiBasis2D) GradVandermonde2D(N int, R, S utils.Vector) (V2Dr, V2Ds utils.Matrix) {
	var (
		Np = (N + 1) * (N + 2) / 2
		Nr = R.Len()
	)
	V2Dr, V2Ds = utils.NewMatrix(Nr, Np), utils.NewMatrix(Nr, Np)
	var sk int
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			ddr, dds := jb2d.GradSimplex2DP(R, S, i, j)
			V2Dr.M.SetCol(sk, ddr)
			V2Ds.M.SetCol(sk, dds)
			sk++
		}
	}
	return
}

func (jb2d *JacobiBasis2D) Simplex2DP(R, S utils.Vector, i, j int) (P []float64) {
	var (
		A, B = RStoAB(R, S)
		Np   = A.Len()
		bd   = B.DataP
	)
	h1 := DG1D.JacobiP(A, 0, jb2d.Beta, i)
	h2 := DG1D.JacobiP(B, float64(2*i+1), jb2d.Beta, j)
	P = make([]float64, Np)
	sq2 := math.Sqrt(2)
	for ii := range h1 {
		tv1 := sq2 * h1[ii] * h2[ii]
		tv2 := utils.POW(1-bd[ii], i)
		P[ii] = tv1 * tv2
	}
	return
}

func (jb2d *JacobiBasis2D) GradSimplex2DP(R, S utils.Vector, id, jd int) (ddr, dds []float64) {
	var (
		A, B   = RStoAB(R, S)
		ad, bd = A.DataP, B.DataP
	)
	fa := DG1D.JacobiP(A, 0, jb2d.Beta, id)
	dfa := DG1D.GradJacobiP(A, 0, jb2d.Beta, id)
	gb := DG1D.JacobiP(B, 2*float64(id)+1, jb2d.Beta, jd)
	dgb := DG1D.GradJacobiP(B, 2*float64(id)+1, jb2d.Beta, jd)
	// r-derivative
	// d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-B)) d/da
	ddr = make([]float64, len(gb))
	for i := range ddr {
		ddr[i] = dfa[i] * gb[i]
		if id > 0 {
			ddr[i] *= utils.POW(0.5*(1-bd[i]), id-1)
		}
		// Normalize
		ddr[i] *= math.Pow(2, float64(id)+0.5)
	}
	// s-derivative
	// d/ds = ((1+A)/2)/((1-B)/2) d/da + d/db
	dds = make([]float64, len(gb))
	for i := range dds {
		dds[i] = 0.5 * dfa[i] * gb[i] * (1 + ad[i])
		if id > 0 {
			dds[i] *= utils.POW(0.5*(1-bd[i]), id-1)
		}
		tmp := dgb[i] * utils.POW(0.5*(1-bd[i]), id)
		if id > 0 {
			tmp -= 0.5 * float64(id) * gb[i] * utils.POW(0.5*(1-bd[i]), id-1)
		}
		dds[i] += fa[i] * tmp
		// Normalize
		dds[i] *= math.Pow(2, float64(id)+0.5)
	}
	return
}

func (jb2d *JacobiBasis2D) PolynomialTerm(r, s float64, i, j int) (P float64) {
	P = jb2d.Simplex2DP(utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s}), i, j)[0]
	return
}

func (jb2d *JacobiBasis2D) PolynomialTermDr(r, s float64, i, j int) (dr float64) {
	ddrV, _ := jb2d.GradSimplex2DP(utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s}), i, j)
	return ddrV[0]
}

func (jb2d *JacobiBasis2D) PolynomialTermDs(r, s float64, i, j int) (ds float64) {
	_, ddsV := jb2d.GradSimplex2DP(utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s}), i, j)
	return ddsV[0]
}

func (jb2d *JacobiBasis2D) GetInterpMatrix(R, S utils.Vector) (Interp utils.Matrix) {
	/*
		Uses Jacobi polynomials as the basis function

		Compose a matrix of interpolating polynomials where each row represents one [r,s] location to be interpolated
		This matrix can then be multiplied by a single vector of function values at the polynomial nodes to produce a
		vector of interpolated values, one for each interpolation location
	*/
	var (
		N  = jb2d.P
		Np = jb2d.Np
	)
	// First compute polynomial terms, used by all polynomials
	polyTerms := make([]float64, R.Len()*Np)
	var sk int
	for ii, r := range R.DataP {
		s := S.DataP[ii]
		for i := 0; i <= N; i++ {
			for j := 0; j <= (N - i); j++ {
				polyTerms[sk] = jb2d.PolynomialTerm(r, s, i, j)
				sk++
			}
		}
	}
	ptV := utils.NewMatrix(R.Len(), Np, polyTerms).Transpose()
	Interp = jb2d.Vinv.Transpose().Mul(ptV).Transpose()
	return
}

func (jb2d *JacobiBasis2D) GetPolynomialEvaluation(r, s float64,
	derivO ...DerivativeDirection) (psi float64) {
	var (
		N     = jb2d.P
		deriv = None
	)
	if len(derivO) > 0 {
		deriv = derivO[0]
	}
	// Compute all polynomial terms and sum to form function value
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			switch deriv {
			case None:
				psi += jb2d.PolynomialTerm(r, s, i, j)
			case Dr:
				psi += jb2d.PolynomialTermDr(r, s, i, j)
			case Ds:
				psi += jb2d.PolynomialTermDs(r, s, i, j)
			}
		}
	}
	return
}

func (jb2d *JacobiBasis2D) GetPolynomialAtJ(r, s float64, j int,
	derivO ...DerivativeDirection) (phi float64) {
	// 2025: This produces the J-th polynomial of the 2D basis, where J is
	// numbered from 0 - Np, and Np = (P+1)(P+2)/2
	var (
		N     = jb2d.P
		deriv = None
		sk    int
	)
	if len(derivO) > 0 {
		deriv = derivO[0]
	}
	// Compute all polynomial terms and sum to form function value
	for i := 0; i <= N; i++ {
		for jj := 0; jj <= (N - i); jj++ {
			if sk == j {
				switch deriv {
				case None:
					phi = jb2d.PolynomialTerm(r, s, i, jj)
				case Dr:
					phi = jb2d.PolynomialTermDr(r, s, i, jj)
				case Ds:
					phi = jb2d.PolynomialTermDs(r, s, i, jj)
				}
				return
			}
			sk++
		}
	}
	return
}

func (jb2d *JacobiBasis2D) GetAllPolynomials(derivO ...DerivativeDirection) (
	PSI utils.Vector) {
	var (
		deriv = None
		m     utils.Matrix
	)
	if len(derivO) > 0 {
		deriv = derivO[0]
	}
	RowSum := func(m utils.Matrix, rowID int) (sum float64) {
		_, ns := m.Dims()
		for i := 0; i < ns; i++ {
			sum += m.At(rowID, i)
		}
		return
	}
	nr, _ := jb2d.V.Dims()
	PSI = utils.NewVector(nr)

	switch deriv {
	case None:
		m = jb2d.V
	case Dr:
		m = jb2d.Vr
	case Ds:
		m = jb2d.Vs
	}
	for i := 0; i < nr; i++ {
		PSI.DataP[i] = RowSum(m, i)
	}
	return
}

type LagrangePolynomial2D struct {
	P      int // Order
	Np     int // Dimension
	R, S   utils.Vector
	Coeffs utils.Matrix
	jb2d   *JacobiBasis2D
}
