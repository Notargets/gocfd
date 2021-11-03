package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

func Vandermonde2D(N int, R, S utils.Vector) (V2D utils.Matrix) {
	V2D = utils.NewMatrix(R.Len(), (N+1)*(N+2)/2)
	var sk int
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			V2D.SetCol(sk, Simplex2DP(R, S, i, j))
			sk++
		}
	}
	return
}

// Purpose  : Compute (x,y) nodes in equilateral triangle for
//            polynomial of order N
func Nodes2D(N int) (x, y utils.Vector) {
	var (
		alpha                                                               float64
		Np                                                                  = (N + 1) * (N + 2) / 2
		L1, L2, L3                                                          utils.Vector
		blend1, blend2, blend3, warp1, warp2, warp3, warpf1, warpf2, warpf3 []float64
	)
	L1, L2, L3, x, y =
		utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np)
	l1d, l2d, l3d, xd, yd := L1.DataP, L2.DataP, L3.DataP, x.DataP, y.DataP
	blend1, blend2, blend3, warp1, warp2, warp3 =
		make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np)

	alpopt := []float64{
		0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
		0.9800, 1.0999, 1.2832, 1.3648, 1.4773,
		1.4959, 1.5743, 1.5770, 1.6223, 1.6258,
	}
	if N < 16 {
		alpha = alpopt[N-1]
	} else {
		alpha = 5. / 3.
	}
	// Create equidistributed nodes on equilateral triangle
	fn := 1. / float64(N)
	var sk int
	for n := 0; n < N+1; n++ {
		for m := 0; m < (N + 1 - n); m++ {
			l1d[sk] = float64(n) * fn
			l3d[sk] = float64(m) * fn
			sk++
		}
	}
	for i := range xd {
		l2d[i] = 1 - l1d[i] - l3d[i]
		xd[i] = l3d[i] - l2d[i]
		yd[i] = (2*l1d[i] - l3d[i] - l2d[i]) / math.Sqrt(3)
		// Compute blending function at each node for each edge
		blend1[i] = 4 * l2d[i] * l3d[i]
		blend2[i] = 4 * l1d[i] * l3d[i]
		blend3[i] = 4 * l1d[i] * l2d[i]
	}
	// Amount of warp for each node, for each edge
	warpf1 = Warpfactor(N, L3.Copy().Subtract(L2))
	warpf2 = Warpfactor(N, L1.Copy().Subtract(L3))
	warpf3 = Warpfactor(N, L2.Copy().Subtract(L1))
	// Combine blend & warp
	for i := range warpf1 {
		warp1[i] = blend1[i] * warpf1[i] * (1 + utils.POW(alpha*l1d[i], 2))
		warp2[i] = blend2[i] * warpf2[i] * (1 + utils.POW(alpha*l2d[i], 2))
		warp3[i] = blend3[i] * warpf3[i] * (1 + utils.POW(alpha*l3d[i], 2))
	}
	// Accumulate deformations associated with each edge
	for i := range xd {
		xd[i] += warp1[i] + math.Cos(2*math.Pi/3)*warp2[i] + math.Cos(4*math.Pi/3)*warp3[i]
		yd[i] += math.Sin(2*math.Pi/3)*warp2[i] + math.Sin(4*math.Pi/3)*warp3[i]
	}
	return
}

func Warpfactor(N int, rout utils.Vector) (warpF []float64) {
	var (
		Nr   = rout.Len()
		Pmat = utils.NewMatrix(N+1, Nr)
	)
	// Compute LGL and equidistant node distribution
	LGLr := DG1D.JacobiGL(0, 0, N)
	req := utils.NewVector(N+1).Linspace(-1, 1)
	Veq := DG1D.Vandermonde1D(N, req)
	// Evaluate Lagrange polynomial at rout
	for i := 0; i < (N + 1); i++ {
		Pmat.M.SetRow(i, DG1D.JacobiP(rout, 0, 0, i))
	}
	Lmat := Veq.Transpose().LUSolve(Pmat)
	// Compute warp factor
	warp := Lmat.Transpose().Mul(LGLr.Subtract(req).ToMatrix())
	// Scale factor
	zerof := rout.Copy().Apply(func(val float64) (res float64) {
		if math.Abs(val) < (1.0 - (1e-10)) {
			res = 1.
		}
		return
	})
	sf := zerof.Copy().ElMul(rout).Apply(func(val float64) (res float64) {
		res = 1 - val*val
		return
	})
	w2 := warp.Copy()
	warp.ElDiv(sf.ToMatrix()).Add(w2.ElMul(zerof.AddScalar(-1).ToMatrix()))
	warpF = warp.DataP
	return
}

func GradVandermonde2D(N int, R, S utils.Vector) (V2Dr, V2Ds utils.Matrix) {
	var (
		Np = (N + 1) * (N + 2) / 2
		Nr = R.Len()
	)
	V2Dr, V2Ds = utils.NewMatrix(Nr, Np), utils.NewMatrix(Nr, Np)
	var sk int
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			ddr, dds := GradSimplex2DP(R, S, i, j)
			V2Dr.M.SetCol(sk, ddr)
			V2Ds.M.SetCol(sk, dds)
			sk++
		}
	}
	return
}

func Simplex2DP(R, S utils.Vector, i, j int) (P []float64) {
	var (
		A, B = RStoAB(R, S)
		Np   = A.Len()
		bd   = B.DataP
	)
	h1 := DG1D.JacobiP(A, 0, 0, i)
	h2 := DG1D.JacobiP(B, float64(2*i+1), 0, j)
	P = make([]float64, Np)
	sq2 := math.Sqrt(2)
	for ii := range h1 {
		tv1 := sq2 * h1[ii] * h2[ii]
		tv2 := utils.POW(1-bd[ii], i)
		P[ii] = tv1 * tv2
	}
	return
}

func GradSimplex2DP(R, S utils.Vector, id, jd int) (ddr, dds []float64) {
	var (
		A, B   = RStoAB(R, S)
		ad, bd = A.DataP, B.DataP
	)
	fa := DG1D.JacobiP(A, 0, 0, id)
	dfa := DG1D.GradJacobiP(A, 0, 0, id)
	gb := DG1D.JacobiP(B, 2*float64(id)+1, 0, jd)
	dgb := DG1D.GradJacobiP(B, 2*float64(id)+1, 0, jd)
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

func Simplex2DPTerm(r, s float64, i, j int) (P float64) {
	P = Simplex2DP(utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s}), i, j)[0]
	return
}

func GradSimplex2DPTerm(r, s float64, i, j int) (ddr, dds float64) {
	ddrV, ddsV := GradSimplex2DP(utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s}), i, j)
	ddr, dds = ddrV[0], ddsV[0]
	return
}

func RStoAB(R, S utils.Vector) (a, b utils.Vector) {
	var (
		Np     = R.Len()
		rd, sd = R.DataP, S.DataP
	)
	ad, bd := make([]float64, Np), make([]float64, Np)
	for n, sval := range sd {
		/*
			if sval != 1 {
				ad[n] = 2*(1+rd[n])/(1-sval) - 1
			} else {
				ad[n] = -1
			}
			bd[n] = sval
		*/
		ad[n], bd[n] = rsToab(rd[n], sval)
	}
	a, b = utils.NewVector(Np, ad), utils.NewVector(Np, bd)
	return
}

func rsToab(r, s float64) (a, b float64) {
	if s != 1 {
		a = 2*(1+r)/(1-s) - 1
	} else {
		a = -1
	}
	b = s
	return
}

// function [r,s] = xytors(x,y)
// Purpose : Transfer from (x,y) in equilateral triangle
//           to (r,s) coordinates in standard triangle
func XYtoRS(x, y utils.Vector) (r, s utils.Vector) {
	r, s = utils.NewVector(x.Len()), utils.NewVector(x.Len())
	var (
		xd, yd = x.DataP, y.DataP
		rd, sd = r.DataP, s.DataP
	)
	sr3 := math.Sqrt(3)
	for i := range xd {
		l1 := (sr3*yd[i] + 1) / 3
		l2 := (-3*xd[i] - sr3*yd[i] + 2) / 6
		l3 := (3*xd[i] - sr3*yd[i] + 2) / 6
		rd[i] = -l2 + l3 - l1
		sd[i] = -l2 - l3 + l1
	}
	return
}

func CalculateElementLocalGeometry(EToV utils.Matrix, VX, VY, R, S utils.Vector) (X, Y utils.Matrix) {
	/*
		For input values of vector field [R,S], transform them into element local [X,Y]
	*/
	va, vb, vc := EToV.Col(0), EToV.Col(1), EToV.Col(2)
	X = R.Copy().Add(S).Scale(-1).Outer(VX.SubsetIndex(va.ToIndex())).Add(
		R.Copy().AddScalar(1).Outer(VX.SubsetIndex(vb.ToIndex()))).Add(
		S.Copy().AddScalar(1).Outer(VX.SubsetIndex(vc.ToIndex()))).Scale(0.5)
	Y = R.Copy().Add(S).Scale(-1).Outer(VY.SubsetIndex(va.ToIndex())).Add(
		R.Copy().AddScalar(1).Outer(VY.SubsetIndex(vb.ToIndex()))).Add(
		S.Copy().AddScalar(1).Outer(VY.SubsetIndex(vc.ToIndex()))).Scale(0.5)
	return
}

type LagrangeBasis1D struct {
	P       int       // Order
	Np      int       // Dimension of basis = N+1
	Weights []float64 // Barycentric weights, one per basis polynomial
	Nodes   []float64 // Nodes at which basis is defined
}

func NewLagrangeBasis1D(R []float64) (lb *LagrangeBasis1D) {
	/*
		At a given order P, there are (P+1) basis polynomials representing that order
		To recover a basis polynomial we need to specifiy:
		`	P = The order of the basis
			j = The basis polynomial number within the basis
			R = The points used to define the basis, (P+1) dimension
	*/
	lb = &LagrangeBasis1D{
		P:       len(R) - 1,
		Np:      len(R),
		Weights: make([]float64, len(R)),
		Nodes:   R,
	}
	// Calculate the weight for each basis function j
	for j := 0; j < lb.Np; j++ {
		lb.Weights[j] = 1.
	}
	for j := 0; j < lb.Np; j++ {
		for i := 0; i < lb.Np; i++ {
			if i != j {
				lb.Weights[j] /= R[j] - R[i]
			}
		}
	}
	return
}

func (lb *LagrangeBasis1D) GetInterpolationMatrix(R []float64) (im utils.Matrix) {
	/*
			Provided function values at each of the P+1 nodes, interpolate a new function value at location r
			Note that the points in R are not necessarily the defining points of the basis, and are not necessarily at the
		    same points within F, the provided set of function values at the nodes of the basis
	*/
	var (
		fj = make([]float64, len(R)) // temporary storage for each basis function evaluation
	)
	im = utils.NewMatrix(len(R), lb.Np) // Rows are for evaluation points, columns for basis
	for j := 0; j < lb.Np; j++ {        // For each basis function
		fj = lb.EvaluateBasisPolynomial(R, j)
		for i, val := range fj {
			im.Set(i, j, val)
		}
	}
	return
}

func (lb *LagrangeBasis1D) Interpolate(R []float64, F []float64) (f []float64) {
	/*
			Provided function values at each of the P+1 nodes, interpolate a new function value at location r
			Note that the points in R are not necessarily the defining points of the basis, and are not necessarily at the
		    same points within F, the provided set of function values at the nodes of the basis
	*/
	var (
		fj = make([]float64, len(R)) // temporary storage for each basis function evaluation
	)
	for j := 0; j < lb.Np; j++ { // For each basis function
		fj = lb.EvaluateBasisPolynomial(R, j)
		for i := range R {
			f[i] += fj[i] * F[j]
		}
	}
	return
}

func (lb *LagrangeBasis1D) EvaluateBasisPolynomial(R []float64, j int) (f []float64) {
	/*
		This evaluates a single basis polynomial (the jth) within the basis for order P at all points in R
		Note that the points in R are not necessarily the defining points of the basis
	*/
	f = make([]float64, len(R))
	for i, r := range R {
		f[i] = lb.evaluateL(r) * lb.Weights[j]
		if math.Abs(r-lb.Nodes[j]) < 0.0000000001 {
			f[i] = 1.
		} else {
			f[i] /= (r - lb.Nodes[j])
		}
	}
	return
}

func (lb *LagrangeBasis1D) evaluateL(r float64) (f float64) {
	/*
		This is the polynomial term in the Barycentric version of the Lagrange polynomial basis
		It is not specific to the jth polynomial, but applies to all the individual basis polynomials
	*/
	f = 1.
	for _, rr := range lb.Nodes {
		f *= (r - rr)
	}
	return
}

type LagrangeBasis2D struct {
	P, Np                      int          // Order
	RNodes, SNodes             []float64    // Nodes at which basis is defined
	V, Vinv, Vr, Vs            utils.Matrix // Inverse of Vandermonde matrix and derivatives for orthogonal non-Lagrange 2D polynomial basis
	Interp, InterpDR, InterpDS utils.Matrix
}

func NewLagrangeBasis2D(P int, R, S utils.Vector) (lb2d *LagrangeBasis2D) {
	/*
			From Karniadakis and Sherwin's "Spectral/hp Element Methods for CFD" on page 124:
			"Since there is not a closed form expression for the Lagrange polynomial through an arbitrary set of points in the
			triangular region, it is necessary to express the Lagrange polynomial in terms of another polynomial which has a
			closed form definition..."

			Here we use the 2D simplex orthogonal polynomial from Hesthaven to produce a generalized Vandermonde matrix, from
		    which we can use the inverse of the Vandermonde to multiply the other basis at locations of our choosing to obtain
			the Lagrange basis coefficients.
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
	}
	lb2d.V = Vandermonde2D(P, R, S)
	lb2d.Vinv = lb2d.V.InverseWithCheck()
	lb2d.Vr, lb2d.Vs = GradVandermonde2D(P, R, S)
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
			Interp.SetCol(sk, Simplex2DP(R, S, I, J))
			sk++
		}
	}
	Interp = lb2d.Vinv.Transpose().Mul(Interp.Transpose())
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
			dR, dS := GradSimplex2DP(R, S, I, J)
			InterpDR.SetCol(sk, dR)
			InterpDS.SetCol(sk, dS)
			sk++
		}
	}
	InterpDR = lb2d.Vinv.Transpose().Mul(InterpDR.Transpose())
	InterpDS = lb2d.Vinv.Transpose().Mul(InterpDS.Transpose())
	return
}

func (lb2d *LagrangeBasis2D) EvaluateBasisPolynomialTerm(R, S utils.Vector, i, j int) (P []float64) {
	/*
		[i,j] are the coordinates of the basis polynomial term
		R and S are the location to get values from the polynomial
	*/
	Interp := lb2d.GetInterpMatrix(R, S)
	P = Interp.Row(lb2d.getTermNumber(i, j)).DataP
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
	return
}
