package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type JacobiBasis2D struct {
	P                       int // Order
	Np                      int // Dimension
	Alpha, Beta             float64
	V, Vinv, Vr, Vs         utils.Matrix
	OrderAtJ                []int
	Order2DAtJ              [][2]int // Polynomial order in each direction at J modes
	VGS, VinvGS, VrGS, VsGS utils.Matrix
	GSTransfer              utils.Matrix
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

	// 2) compute raw modal mass: diag(Vraw^T W Vraw)
	w := WilliamsShunnJamesonWeights(P)
	W := utils.NewDiagMatrix(len(w), w)

	// 3) GS orthonormalize Vraw -> Q
	jb.VGS = OrthonormalizeBasis(jb.V, W)
	// rebuild Vinv = Q^T W
	jb.VinvGS = jb.VGS.Transpose().Mul(W)
	jb.GSTransfer = jb.Vinv.Mul(jb.VGS)

	// 4) scale VrRaw,VsRaw to match Q
	jb.VrGS = jb.Vr.Mul(jb.GSTransfer)
	jb.VsGS = jb.Vs.Mul(jb.GSTransfer)
	return jb
}

// GetInterpMatrixGS builds interpolation from solution nodes -> given (R_e,S_e) edge nodes.
// It uses the same orthonormal basis Q so it's consistent with modal transforms.
func (jb *JacobiBasis2D) GetInterpMatrixGS(R_e, S_e utils.Vector) utils.Matrix {
	// raw Vandermonde at edge
	Vedge := jb.Vandermonde2D(jb.P, R_e, S_e)
	return Vedge.Mul(jb.GSTransfer).Mul(jb.VinvGS)
}

// OrthonormalizeBasis as before (weighted GS)
func OrthonormalizeBasis(V, W utils.Matrix) utils.Matrix {
	Nq, Np := V.Dims()
	Q := utils.NewMatrix(Nq, Np)
	for j := 0; j < Np; j++ {
		// copy V[:,j]
		for i := 0; i < Nq; i++ {
			Q.Set(i, j, V.At(i, j))
		}
		// subtract projections
		for i := 0; i < j; i++ {
			dot := 0.0
			for k := 0; k < Nq; k++ {
				dot += Q.At(k, j) * W.At(k, k) * Q.At(k, i)
			}
			for k := 0; k < Nq; k++ {
				Q.Set(k, j, Q.At(k, j)-dot*Q.At(k, i))
			}
		}
		// normalize
		sum := 0.0
		for k := 0; k < Nq; k++ {
			sum += Q.At(k, j) * W.At(k, k) * Q.At(k, j)
		}
		inv := 1.0 / math.Sqrt(sum)
		for k := 0; k < Nq; k++ {
			Q.Set(k, j, Q.At(k, j)*inv)
		}
	}
	return Q
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

func (jb2d *JacobiBasis2D) GetModInterpMatrix(R, S utils.Vector,
	Nu, p float64, SmoothingIterations int) (Interp utils.Matrix) {
	/*
		Uses Jacobi polynomials as the basis function

		Compose a matrix of interpolating polynomials where each row represents one [R,S] location to be interpolated
		This matrix can then be multiplied by a single vector of function values at the polynomial nodes to produce a
		vector of interpolated values, one for each interpolation location
	*/
	var (
		N  = jb2d.P
		Np = jb2d.Np
	)
	if SmoothingIterations < 1 {
		SmoothingIterations = 1
	}
	CoefMods := jb2d.ModalFilter2D(Nu, p)
	// First compute polynomial terms, used by all polynomials
	polyTerms := make([]float64, R.Len()*Np)
	var sk int
	for ii, r := range R.DataP {
		s := S.DataP[ii]
		var sk2 int
		for i := 0; i <= N; i++ {
			for j := 0; j <= (N - i); j++ {
				CoefMod := math.Pow(CoefMods[sk2], float64(SmoothingIterations))
				polyTerms[sk] = CoefMod * jb2d.PolynomialTerm(r, s, i, j)
				sk++
				sk2++
			}
		}
	}
	ptV := utils.NewMatrix(R.Len(), Np, polyTerms).Transpose()
	Interp = jb2d.Vinv.Transpose().Mul(ptV).Transpose()
	return
}

func (jb2d *JacobiBasis2D) GetModulationMatrix(Nu, p float64) (ModMat utils.Matrix) {
	// This returns a ModMat matrix that will filter high modes out of an
	// existing field on Np nodal points in one matrix multiplication
	var (
		Np = jb2d.Np
		F  = utils.NewDiagMatrix(Np, jb2d.ModalFilter2D(Nu, p))
	)
	ModMat = jb2d.V.Mul(F).Mul(jb2d.Vinv)
	return
}

func (jb2d *JacobiBasis2D) ModalFilter2D(Nu, p float64) (CoeffModifier []float64) {
	// Tunables:
	// - Nu is a user-defined dissipation strength (tunable),
	// 0.1 to 0.5	Higher = stronger overall dissipation (aggressiveness)
	// - p controls sharpness (typically 1 to 4 — sharper if you only want to hit the highest modes).
	// 1 to 4	Higher = sharper cutoff (more selective to highest modes)
	var (
		P      = jb2d.P
		degree = jb2d.ModalDegree2D(P)
		Np     = len(degree)
	)
	CoeffModifier = make([]float64, Np)
	for i := 0; i < Np; i++ {
		CoeffModifier[i] = 1. - Nu*(math.Pow(degree[i]/float64(P), 2.*p))
	}
	return
}

func (jb2d *JacobiBasis2D) ModalDegree2D(P int) (degree []float64) {
	var (
		Np = (P + 1) * (P + 2) / 2
	)
	degree = make([]float64, Np)
	var sk int
	for i := 0; i <= P; i++ {
		for j := 0; j <= P-i; j++ {
			degree[sk] = math.Max(float64(i), float64(j))
			sk++
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
