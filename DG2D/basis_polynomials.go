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
	OrderAtJ        []int
	Order2DAtJ      [][2]int // Polynomial order in each direction at J modes
}

func NewJacobiBasis2D(P int, R, S utils.Vector, alpha, beta float64) *JacobiBasis2D {
	// This routine now calculates an Orthonormal basis using a Gramm-Schmidt
	// process. This assumes the use of the Williams-Shun-Jameson point
	// distribution for the element with its associated quadrature
	Np := (P + 1) * (P + 2) / 2
	jb := &JacobiBasis2D{
		P:          P,
		Np:         Np,
		Alpha:      alpha,
		Beta:       beta,
		OrderAtJ:   make([]int, Np),
		Order2DAtJ: make([][2]int, Np),
	}
	// fill OrderAtJ, Order2DAtJ
	sk := 0
	for i := 0; i <= P; i++ {
		for j := 0; j <= P-i; j++ {
			jb.OrderAtJ[sk] = i + j
			jb.Order2DAtJ[sk] = [2]int{i, j}
			sk++
		}
	}

	// 1) build raw Vandermonde and gradients
	jb.V = jb.Vandermonde2D(P, R, S)
	jb.Vr, jb.Vs = jb.GradVandermonde2D(P, R, S)
	jb.Vinv = jb.V.InverseWithCheck()
	return jb
}

// GetInterpMatrix builds interpolation from solution nodes -> given (R_e,S_e) edge nodes.
// It uses the same orthonormal basis Q so it's consistent with modal transforms.
func (jb *JacobiBasis2D) GetInterpMatrix(R_e, S_e utils.Vector) utils.Matrix {
	// raw Vandermonde at edge
	Vedge := jb.Vandermonde2D(jb.P, R_e, S_e)
	return Vedge.Mul(jb.Vinv)
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

func (jb2d *JacobiBasis2D) GradSimplex2DP(R, S utils.Vector, id, jd int) (ddr, dds []float64) {
	var (
		A, B   = RStoAB(R, S)
		ad, bd = A.DataP, B.DataP
	)
	fa := DG1D.JacobiP(A, 0, jb2d.Beta, id)
	dfa := DG1D.GradJacobiP(A, 0, jb2d.Beta, id)
	gb := DG1D.JacobiP(B, 2*float64(id)+1, jb2d.Beta, jd)
	dgb := DG1D.GradJacobiP(B, 2*float64(id)+1, jb2d.Beta, jd)
	// R-derivative
	// d/dr = da/dr d/da + db/dr d/db = (2/(1-S)) d/da = (2/(1-B)) d/da
	ddr = make([]float64, len(gb))
	for i := range ddr {
		ddr[i] = dfa[i] * gb[i]
		if id > 0 {
			ddr[i] *= utils.POW(0.5*(1-bd[i]), id-1)
		}
		// Normalize
		ddr[i] *= math.Pow(2, float64(id)+0.5)
	}
	// S-derivative
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
	i, jj := jb2d.Order2DAtJ[j][0], jb2d.Order2DAtJ[j][1]
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

func (jb2d *JacobiBasis2D) GetOrthogonalPolynomialAtJ(r, s float64, j int,
	derivO ...DerivativeDirection) (phi float64) {
	// The coefficients of the J-th orthogonal polynomial are columns of Vinv
	// The value of the J-th polynomial evaluated at it'S defining point is 1
	var (
		deriv = None
	)
	if len(derivO) > 0 {
		deriv = derivO[0]
	}
	for i := 0; i < jb2d.Np; i++ {
		ii, jj := jb2d.Order2DAtJ[i][0], jb2d.Order2DAtJ[i][1]
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

func Lagrange1DPoly(t float64, R []float64, j int, derivO ...DerivativeDirection) (p float64) {
	var (
		deriv = false
	)
	// This is a 1D polynomial, but the nodes are either in the RDir or SDir direction within 2D
	if len(derivO) != 0 {
		deriv = true
	}
	if !deriv {
		var (
			Np1D = len(R)
			XJ   = R[j]
			XI   = t
		)
		if j > Np1D-1 || j < 0 {
			panic("value of j larger than array or less than zero")
		}
		p = 1
		for m, XM := range R {
			if m != j {
				p *= (XI - XM) / (XJ - XM)
			}
		}
	} else {
		var (
			XJ = R[j]
			X  = t
		)
		for i, XI := range R {
			if i != j {
				pp := 1.
				for m, XM := range R {
					if m != i && m != j {
						pp *= (X - XM) / (XJ - XM)
					}
				}
				p += pp / (XJ - XI)
			}
		}
	}
	return
}

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
	// The value of the J-th polynomial evaluated at it'S defining point is 1
	for i := 0; i < jb1d.Np; i++ {
		phi += jb1d.JacobiP(r, i, derivO...) * jb1d.Vinv.At(i, j)
	}
	return
}
