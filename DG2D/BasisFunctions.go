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
	h1 := DG1D.JacobiP(A, jb2d.Alpha, jb2d.Beta, i)
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

	// Compute Jacobi polynomials and their derivatives
	fa := DG1D.JacobiP(A, jb2d.Alpha, jb2d.Beta, id)
	dfa := DG1D.GradJacobiP(A, jb2d.Alpha, jb2d.Beta, id) // Correct α use here
	gb := DG1D.JacobiP(B, 2*float64(id)+1, jb2d.Beta, jd)
	dgb := DG1D.GradJacobiP(B, 2*float64(id)+1, jb2d.Beta, jd)

	// Compute r-derivative: (2 / (1 - B)) d/da
	ddr = make([]float64, len(gb))
	for i := range ddr {
		ddr[i] = (2 / (1 - bd[i])) * dfa[i] * gb[i] // Include transformation factor
		if id > 0 {
			ddr[i] *= utils.POW(0.5*(1-bd[i]), id-1)
		}
		// Normalize
		ddr[i] *= math.Pow(2, float64(id)+0.5)
	}

	// Compute s-derivative: ((1 + A) / (1 - B)) d/da + d/db
	dds = make([]float64, len(gb))
	for i := range dds {
		// First term: ((1+A)/2)/((1-B)/2) * d/da
		dds[i] = ((1 + ad[i]) / (1 - bd[i])) * 0.5 * dfa[i] * gb[i]
		if id > 0 {
			dds[i] *= utils.POW(0.5*(1-bd[i]), id-1)
		}

		// Second term: d/db with explicit α adjustment
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

type LagrangeBasis2D struct {
	P, Np          int       // Order
	RNodes, SNodes []float64 // Nodes at which basis is defined
	JB2D           *JacobiBasis2D
}

func NewLagrangeBasis2D(P int, R, S utils.Vector) (lb2d *LagrangeBasis2D) {
	/*
		From Karniadakis and Sherwin's "Spectral/hp Element Methods for CFD" on page 124:
		"Since there is not a closed form expression for the Lagrange polynomial through an arbitrary set of points in the
		triangular region, it is necessary to express the Lagrange polynomial in terms of another polynomial which has a
		closed form definition..."

		Here we use the 2D simplex orthogonal polynomial from Hesthaven to produce a generalized Vandermonde matrix,
		from which we can use the inverse of the Vandermonde to multiply the Hestaven basis at locations of our
		choosing to obtain the Lagrange basis coefficients at that location, indirectly through the original basis.
	*/
	var (
		Np = (P + 1) * (P + 2) / 2
	)
	// Dimension: Np = (N+1)*(N+2)/2
	if Np != R.Len() || Np != S.Len() {
		err := fmt.Errorf("wrong length of arrays defining basis nodes, should be %d, is [Rlen,Slen] = [%d,%d]",
			Np, R.Len(), S.Len())
		panic(err)
	}
	lb2d = &LagrangeBasis2D{
		P:      P,
		Np:     Np,
		RNodes: R.DataP,
		SNodes: S.DataP,
		JB2D:   NewJacobiBasis2D(P, R, S, 0, 0),
	}
	return
}

func (lb2d *LagrangeBasis2D) GetInterpMatrix(R, S utils.Vector) (Interp utils.Matrix) {
	var (
		Np = lb2d.Np
		N  = lb2d.P
	)
	Interp = utils.NewMatrix(R.Len(), Np) // Will transpose after fill
	var sk int
	for I := 0; I <= N; I++ {
		for J := 0; J <= N-I; J++ {
			Interp.SetCol(sk, lb2d.JB2D.Simplex2DP(R, S, I, J))
			sk++
		}
	}
	Interp = lb2d.JB2D.Vinv.Transpose().Mul(Interp.Transpose()).Transpose()
	return
}

func (lb2d *LagrangeBasis2D) GetGradInterpMatrices(R, S utils.Vector) (InterpDR, InterpDS utils.Matrix) {
	var (
		Np = lb2d.Np
		N  = lb2d.P
	)
	InterpDR = utils.NewMatrix(R.Len(), Np) // Will transpose after fill
	InterpDS = utils.NewMatrix(R.Len(), Np) // Will transpose after fill
	var sk int
	for I := 0; I <= N; I++ {
		for J := 0; J <= N-I; J++ {
			dR, dS := lb2d.JB2D.GradSimplex2DP(R, S, I, J)
			InterpDR.SetCol(sk, dR)
			InterpDS.SetCol(sk, dS)
			sk++
		}
	}
	InterpDR = lb2d.JB2D.Vinv.Transpose().Mul(InterpDR.Transpose()).Transpose()
	InterpDS = lb2d.JB2D.Vinv.Transpose().Mul(InterpDS.Transpose()).Transpose()
	return
}

func (lb2d *LagrangeBasis2D) BasisPolynomial(R, S utils.Vector, i, j int) (P []float64) {
	/*
		[i,j] are the coordinates of the basis polynomial term
		R and S are the location to get values from the polynomial
	*/
	Interp := lb2d.GetInterpMatrix(R, S)
	// fmt.Printf("P = %d, term number[%d,%d] = %d\n", lb2d.P, i, j, lb2d.getTermNumber(i, j))
	P = Interp.Col(lb2d.getTermNumber(i, j)).DataP
	return
}

func (lb2d *LagrangeBasis2D) PolynomialTerm(r, s float64, i, j int) (p float64) {
	return lb2d.BasisPolynomial(utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s}), i, j)[0]
}

func (lb2d *LagrangeBasis2D) PolynomialTermDr(r, s float64, i, j int) (dr float64) {
	DR, _ := lb2d.GradBasisPolynomial(utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s}), i, j)
	return DR[0]
}

func (lb2d *LagrangeBasis2D) PolynomialTermDs(r, s float64, i, j int) (ds float64) {
	_, DS := lb2d.GradBasisPolynomial(utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s}), i, j)
	return DS[0]
}

func (lb2d *LagrangeBasis2D) GradBasisPolynomial(R, S utils.Vector, i, j int) (DrTerms, DsTerms []float64) {
	/*
		[i,j] are the coordinates of the basis polynomial term
		R and S are the location to get values from the polynomial
	*/
	InterpDR, InterpDS := lb2d.GetGradInterpMatrices(R, S)
	InterpDR = InterpDR
	InterpDS = InterpDS
	DrTerms = InterpDR.Col(lb2d.getTermNumber(i, j)).DataP
	DsTerms = InterpDS.Col(lb2d.getTermNumber(i, j)).DataP
	return
}

func (lb2d *LagrangeBasis2D) getTermNumber(i, j int) (sk int) {
	var (
		N = lb2d.P
	)
	for I := 0; I <= N; I++ {
		for J := 0; J <= N-I; J++ {
			if I == i && J == j {
				return
			}
			sk++
		}
	}
	sk = -1
	return
}
