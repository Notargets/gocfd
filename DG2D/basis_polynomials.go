package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type JacobiBasis1D struct {
	P           int // Order
	Np          int // Dimension
	Alpha, Beta float64
	R           utils.Vector
	V, Vinv, Vr utils.Matrix
}

func NewJacobiBasis1D(P int, R utils.Vector, Alpha, Beta float64) (jb1d *JacobiBasis1D) {
	jb1d = &JacobiBasis1D{
		P:     P,
		Np:    P + 1,
		Alpha: Alpha,
		Beta:  Beta,
		R:     R,
	}

	if R.Len() != jb1d.Np {
		fmt.Printf("Length of R:%d, required length:%d\n", R.Len(), jb1d.Np)
		panic("Mismatch of length for basis coordinates")
	}
	jb1d.V = jb1d.Vandermonde1D()
	jb1d.Vinv = jb1d.V.InverseWithCheck()
	jb1d.Vr = jb1d.GradVandermonde1D()
	return
}

func (jb1d *JacobiBasis1D) JacobiP(r float64, j int,
	derivO ...DerivativeDirection) float64 {
	if len(derivO) == 0 {
		return DG1D.JacobiP(utils.NewVector(1, []float64{r}),
			jb1d.Alpha, jb1d.Beta, j)[0]
	} else {
		return DG1D.GradJacobiP(utils.NewVector(1, []float64{r}),
			jb1d.Alpha, jb1d.Beta, j)[0]
	}
}

func (jb1d *JacobiBasis1D) Vandermonde1D() utils.Matrix {
	return DG1D.Vandermonde1D(jb1d.P, jb1d.R)
}

func (jb1d *JacobiBasis1D) GradVandermonde1D() utils.Matrix {
	return DG1D.GradVandermonde1D(jb1d.P, jb1d.R)
}

func (jb1d *JacobiBasis1D) GetOrthogonalPolynomialAtJ(r float64, j int,
	derivO ...DerivativeDirection) (phi float64) {
	// The coefficients of the J-th orthogonal polynomial are columns of Vinv
	// The value of the J-th polynomial evaluated at it's defining point is 1
	for i := 0; i < jb1d.Np; i++ {
		phi += jb1d.JacobiP(r, i, derivO...) * jb1d.Vinv.At(i, j)
	}
	return
}

type JacobiBasis2D struct {
	P               int // Order
	Np              int // Dimension
	Alpha, Beta     float64
	V, Vinv, Vr, Vs utils.Matrix
	IJ              [][2]int // I,J coordinates of basis polynomial indexed on j
}

func NewJacobiBasis2D(P int, R, S utils.Vector, Alpha, Beta float64) (jb2d *JacobiBasis2D) {
	jb2d = &JacobiBasis2D{
		P:     P,
		Np:    (P + 1) * (P + 2) / 2,
		Alpha: Alpha,
		Beta:  Beta,
	}
	jb2d.IJ = make([][2]int, jb2d.Np)
	var sk int
	for i := 0; i <= jb2d.P; i++ {
		for j := 0; j <= jb2d.P-i; j++ {
			jb2d.IJ[sk] = [2]int{i, j}
			sk++
		}
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
		deriv = None
	)
	if len(derivO) > 0 {
		deriv = derivO[0]
	}
	i, jj := jb2d.IJ[j][0], jb2d.IJ[j][1]
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

func (jb2d *JacobiBasis2D) GetOrthogonalPolynomialAtJ(r, s float64, j int,
	derivO ...DerivativeDirection) (phi float64) {
	// The coefficients of the J-th orthogonal polynomial are columns of Vinv
	// The value of the J-th polynomial evaluated at it's defining point is 1
	var (
		deriv = None
	)
	if len(derivO) > 0 {
		deriv = derivO[0]
	}
	for i := 0; i < jb2d.Np; i++ {
		ii, jj := jb2d.IJ[i][0], jb2d.IJ[i][1]
		switch deriv {
		case None:
			phi += jb2d.PolynomialTerm(r, s, ii, jj) * jb2d.Vinv.At(i, j)
		case Dr:
			phi += jb2d.PolynomialTermDr(r, s, ii, jj) * jb2d.Vinv.At(i, j)
		case Ds:
			phi += jb2d.PolynomialTermDs(r, s, ii, jj) * jb2d.Vinv.At(i, j)
		}
	}
	return
}
