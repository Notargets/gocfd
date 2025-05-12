package DG1D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
)

type NODE_TYPE uint

const (
	GAUSS NODE_TYPE = iota
	GAUSS_LOBATO
)

type Elements1D struct {
	K, Np, Nfp, NFaces                int
	R, VX, FMask                      utils.Vector
	EToV, EToE, EToF                  utils.Matrix
	X, Dr, Rx, FScale, NX, LIFT       utils.Matrix
	V, Vinv                           utils.Matrix
	VmapM, VmapP, VmapB, VmapI, VmapO utils.Index
	MapB, MapI, MapO                  utils.Index
}

func NewElements1D(N int, VX utils.Vector, EToV utils.Matrix, ntA ...NODE_TYPE) (el *Elements1D) {
	var (
		K, NFaces = EToV.Dims()
		Nfp       = 1 // One point per face in 1D
	)
	// N is the polynomial degree, Np is the number of interpolant points = N+1
	el = &Elements1D{
		K:      K,
		Np:     N + 1,
		Nfp:    Nfp,
		NFaces: NFaces,
		VX:     VX,
		EToV:   EToV,
	}
	var nt NODE_TYPE
	if len(ntA) == 0 {
		nt = GAUSS_LOBATO
	} else {
		nt = ntA[0]
	}
	el.Startup1D(nt)
	el.V.SetReadOnly("V")
	el.Vinv.SetReadOnly("Vinv")
	el.EToV.SetReadOnly("EToV")
	el.EToE.SetReadOnly("EToE")
	el.EToF.SetReadOnly("EToF")
	el.X.SetReadOnly("R")
	el.Dr.SetReadOnly("Dr")
	el.Rx.SetReadOnly("Rx")
	el.FScale.SetReadOnly("FScale")
	el.NX.SetReadOnly("NX")
	el.LIFT.SetReadOnly("LIFT")
	return
}

func JacobiGL(alpha, beta float64, N int) (R utils.Vector) {
	var (
		x    = make([]float64, N+1)
		xint utils.Vector
	)
	if N == 1 {
		x[0] = -1
		x[1] = 1
		R = utils.NewVector(N+1, x)
		return
	}
	xint, _ = JacobiGQ(alpha+1, beta+1, N-2)
	x[0] = -1
	x[N] = 1
	var iter int
	dataXint := xint.V.RawVector().Data
	for i := 1; i < N; i++ {
		// x[i] = xint.AtVec(iter)
		x[i] = dataXint[iter]
		iter++
	}
	R = utils.NewVector(len(x), x)
	return
}

func LegendreZeros(N int, reverseO ...bool) (X []float64) {
	/*
		For a given order N, there are N locations within [-1,1] that are the zeros of the Legendre polynomial
		These zeros were tabulated within "TABLE OF THE ZEROS OF THE LEGENDRE POLYNOMIALS OF ORDER 1-16 AND THE
		WEIGHT COEFFICIENTS FOR GAUSS' MECHANICAL QUADRATURE FORMULA1" by Lowan, et. al.
	*/
	switch N {
	case 0:
		X = []float64{0}
	case 1:
		X = []float64{-0.577350269189626, 0.577350269189626}
	case 2:
		X = []float64{-0.774596669241483, 0, 0.774596669241483}
	case 3:
		X = []float64{-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053}
	case 4:
		X = []float64{-0.906179845938664, -0.538469310105683, 0.000000000000000, 0.538469310105683, 0.906179845938664}
	case 5:
		X = []float64{
			-0.932469514203152, -0.661209386466265, -0.238619186083197,
			0.238619186083197, 0.661209386466265, 0.932469514203152,
		}
	case 6:
		X = []float64{
			-0.949107912342759, -0.741531185599394, -0.405845151377397,
			0.0,
			0.405845151377397, 0.741531185599394, 0.949107912342759,
		}
	case 7:
		X = []float64{
			-0.960289856497536, -0.796666477413627, -0.525532409916329, -0.183434642495650,
			0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536,
		}
	case 8:
		X = []float64{
			-0.968160239507626, -0.836031107326636, -0.613371432700590, -0.324253423403809,
			0.000000000000000,
			0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626,
		}
	case 9:
		X = []float64{
			-0.973906528517172, -0.865063366688985, -0.679409568299024, -0.433395394129247, -0.148874338981631,
			0.148874338981631, 0.433395394129247, 0.679409568299024, 0.865063366688985, 0.973906528517172,
		}
	case 10:
		X = []float64{
			-0.978228658146057, -0.887062599768095, -0.730152005574049, -0.519096129110681, -0.269543155952345,
			0.000000000000000,
			0.269543155952345, 0.519096129110681, 0.730152005574049, 0.887062599768095, 0.978228658146057,
		}
	case 11:
		X = []float64{
			-0.981560634246719, -0.904117256370475, -0.769902674194305, -0.587317954286617, -0.367831498918180, -0.125333408511469,
			0.125333408511469, 0.367831498918180, 0.587317954286617, 0.769902674194305, 0.904117256370475, 0.981560634246719,
		}
	case 12:
		X = []float64{
			-0.984183054718588, -0.917598399222978, -0.801578090733310, -0.642349339440340, -0.448492751036447, -0.230458315955135,
			0.000000000000000,
			0.230458315955135, 0.448492751036447, 0.642349339440340, 0.801578090733310, 0.917598399222978, 0.984183054718588,
		}
	case 13:
		X = []float64{
			-0.986283808696812, -0.928434883663574, -0.827201315069765, -0.687292904811685, -0.515248636358154,
			-0.319112368927890, -0.108054948707344, 0.108054948707344, 0.319112368927890, 0.515248636358154,
			0.687292904811685, 0.827201315069765, 0.928434883663574, 0.986283808696812,
		}
	case 14:
		X = []float64{
			-0.987992518020485, -0.937273392400706, -0.848206583410427, -0.724417731360170, -0.570972172608539,
			-0.394151347077563, -0.201194093997435,
			0.000000000000000,
			0.201194093997435, 0.394151347077563, 0.570972172608539, 0.724417731360170, 0.848206583410427,
			0.937273392400706, 0.987992518020485,
		}
	case 15:
		X = []float64{
			-0.989400934991650, -0.944575023073233, -0.865631202387832, -0.755404408355003, -0.617876244402644,
			-0.458016777657227, -0.281603550779259, -0.095012509837637,
			0.095012509837637, 0.281603550779259, 0.458016777657227, 0.617876244402644, 0.755404408355003,
			0.865631202387832, 0.944575023073233, 0.989400934991650,
		}
	default:
		err := fmt.Errorf("distribution undefined for order %d", N)
		panic(err)
	}
	var reverse bool
	if len(reverseO) != 0 {
		reverse = reverseO[0]
	}
	if reverse {
		length := len(X)
		X2 := make([]float64, length)
		for i, x := range X {
			X2[length-1-i] = x
		}
		X = X2
	}
	return
}

func JacobiGQ(alpha, beta float64, N int) (X, W utils.Vector) {
	var (
		x, w       []float64
		fac        float64
		h1, d0, d1 []float64
		VVr        *mat.Dense
	)
	if N == 0 {
		x = []float64{-(alpha - beta) / (alpha + beta + 2.)}
		w = []float64{2.}
		return utils.NewVector(len(x), x), utils.NewVector(len(w), w)
	}

	h1 = make([]float64, N+1)
	for i := 0; i < N+1; i++ {
		h1[i] = 2*float64(i) + alpha + beta
	}

	// main diagonal: diag(-1/2*(alpha^2-beta^2)./(h1+2)./h1)
	d0 = make([]float64, N+1)
	fac = -.5 * (alpha*alpha - beta*beta)
	for i := 0; i < N+1; i++ {
		val := h1[i]
		d0[i] = fac / (val * (val + 2.))
	}
	// Handle division by zero
	eps := 1.e-16
	if alpha+beta < 10*eps {
		d0[0] = 0.
	}

	// 1st upper diagonal: diag(2./(h1(1:N)+2).*sqrt((1:N).*((1:N)+alpha+beta) .* ((1:N)+alpha).*((1:N)+beta)./(h1(1:N)+1)./(h1(1:N)+3)),1);
	// for (i=1; i<=N; ++i) { d1(i)=2.0/(h1(i)+2.0)*sqrt(i*(i+alpha+beta)*(i+alpha)*(i+beta)/(h1(i)+1)/(h1(i)+3.0)); }
	var ip1 float64
	d1 = make([]float64, N)
	for i := 0; i < N; i++ {
		ip1 = float64(i + 1)
		val := h1[i]
		d1[i] = 2. / (val + 2.)
		d1[i] *= math.Sqrt(ip1 * (ip1 + alpha + beta) * (ip1 + alpha) * (ip1 + beta) / ((val + 1.) * (val + 3.)))
	}

	JJ := utils.NewSymTriDiagonal(d0, d1)

	var eig mat.EigenSym
	ok := eig.Factorize(JJ, true)
	if !ok {
		panic("eigenvalue decomposition failed")
	}
	x = eig.Values(x)
	X = utils.NewVector(N+1, x)

	VVr = mat.NewDense(len(x), len(x), nil)
	eig.VectorsTo(VVr)
	W = utils.NewVector(len(x), VVr.RawRowView(0)).POW(2).Scale(Gamma0(alpha, beta))
	return X, W
}

func Vandermonde1D(N int, R utils.Vector) (V utils.Matrix) {
	V = utils.NewMatrix(R.Len(), N+1)
	for j := 0; j < N+1; j++ {
		V.SetCol(j, JacobiP(R, 0, 0, j))
	}
	return
}

func JacobiP(r utils.Vector, alpha, beta float64, N int) (p []float64) {
	// This outputs the Jacobi polynomial of degree N, output at the points R
	// Each point input generates an output value that is the complete value
	// of the polynomial,
	// e.g. the summation of all polynomial terms for that location
	var (
		Nc = r.Len()
	)
	rg := 1. / math.Sqrt(Gamma0(alpha, beta))
	if N == 0 {
		p = utils.ConstArray(Nc, rg)
		return
	}
	Np1 := N + 1
	pl := make([]float64, Np1*Nc)
	var iter int
	for i := 0; i < Nc; i++ {
		pl[i+iter] = rg
	}

	iter += Nc // Increment to next row
	ab := alpha + beta
	rg1 := 1. / math.Sqrt(Gamma1(alpha, beta))
	for i := 0; i < Nc; i++ {
		pl[i+iter] = rg1 * ((ab+2.0)*r.AtVec(i)/2.0 + (alpha-beta)/2.0)
	}

	if N == 1 {
		p = pl[iter : iter+Nc]
		return
	}

	a1 := alpha + 1.
	b1 := beta + 1.
	ab1 := ab + 1.
	aold := 2.0 * math.Sqrt(a1*b1/(ab+3.0)) / (ab + 2.0)
	PL := mat.NewDense(Np1, Nc, pl)
	var xrow []float64
	for i := 0; i < N-1; i++ {
		ip1 := float64(i + 1)
		ip2 := float64(ip1 + 1)
		h1 := 2.0*ip1 + ab
		anew := 2.0 / (h1 + 2.0) * math.Sqrt(ip2*(ip1+ab1)*(ip1+a1)*(ip1+b1)/(h1+1.0)/(h1+3.0))
		bnew := -(alpha*alpha - beta*beta) / h1 / (h1 + 2.0)
		x_bnew := utils.NewVector(r.Len()).SetScalar(-bnew)
		x_bnew.Add(r)
		xi := PL.RawRowView(i)
		xip1 := PL.RawRowView(i + 1)
		xrow = make([]float64, len(xi))
		for j := range xi {
			xrow[j] = (-aold*xi[j] + x_bnew.AtVec(j)*xip1[j]) / anew
		}
		PL.SetRow(i+2, xrow)
		aold = anew
	}
	p = PL.RawRowView(N)
	return
}

func GradJacobiP(r utils.Vector, alpha, beta float64, N int) (p []float64) {
	if N == 0 {
		p = make([]float64, r.Len())
		return
	}
	p = JacobiP(r, alpha+1, beta+1, N-1)
	fN := float64(N)
	fac := math.Sqrt(fN * (fN + alpha + beta + 1))
	for i, val := range p {
		p[i] = val * fac
	}
	return
}

func GradVandermonde1D(N int, r utils.Vector) (Vr utils.Matrix) {
	Vr = utils.NewMatrix(r.Len(), N+1)
	for i := 0; i < N+1; i++ {
		Vr.SetCol(i, GradJacobiP(r, 0, 0, i))
	}
	return
}

func Lift1D(V utils.Matrix, Np, Nfaces, Nfp int) (LIFT utils.Matrix) {
	Emat := utils.NewMatrix(Np, Nfaces*Nfp)
	Emat.Set(0, 0, 1)
	Emat.Set(Np-1, 1, 1)
	LIFT = V.Mul(V.Transpose()).Mul(Emat)
	return
}

func Normals1D(Nfaces, Nfp, K int) (NX utils.Matrix) {
	nx := make([]float64, Nfaces*Nfp*K)
	for i := 0; i < K; i++ {
		nx[i] = -1
		nx[i+K] = 1
	}
	NX = utils.NewMatrix(Nfp*Nfaces, K, nx)
	return
}

func GeometricFactors1D(Dr, X utils.Matrix) (J, Rx utils.Matrix) {
	J = Dr.Mul(X)
	Rx = J.Copy().POW(-1)
	return
}

func (el *Elements1D) LagrangeInterpolant(r float64) (Li utils.Matrix) {
	Li = utils.NewMatrix(1, el.R.Len()).AddScalar(1)
	LiData := Li.RawMatrix().Data
	for j := range el.R.RawVector().Data {
		LiData[j] *= LagrangePolyAtJ(r, el.R.DataP, j)
	}
	return
}

func LagrangePolyAtJ(r float64, R []float64, j int) (s float64) {
	/*
		From https://en.wikipedia.org/wiki/Lagrange_polynomial
		This evaluates the Lagrange polynomial at term J for location R[j]

		The equivalent Newton polynomial is more efficient for repetitive usage
	*/
	var (
		k = len(R)
	)
	if j > k || j < 0 {
		panic("value of j larger than array or less than zero")
	}
	s = 1
	for i, rBasis := range R {
		if i == j {
			continue
		}
		metric := (r - rBasis) / (R[j] - rBasis)
		s *= metric
	}
	return
}

// func Lagrange1DPoly(r float64, R []float64, j int, derivO ...int) (p float64) {
// 	/*
// 		Evaluates the Lagrange polynomial or its derivative, defined on points R,
// 		centered on the j-th point R[j], and evaluated at point r
// 		The order of the polynomial is P = len(R)-1
// 	*/
// 	var (
// 		deriv = 0
// 	)
// 	if len(derivO) != 0 {
// 		deriv = derivO[0]
// 	}
// 	switch deriv {
// 	case 0:
// 		var (
// 			Np1D = len(R)
// 			XJ   = R[j]
// 			XI   = r
// 		)
// 		if j > Np1D-1 || j < 0 {
// 			panic("value of j larger than array or less than zero")
// 		}
// 		p = 1
// 		for m, XM := range R {
// 			if m != j {
// 				p *= (XI - XM) / (XJ - XM)
// 			}
// 		}
// 	case 1:
// 		/*
// 			From https://math.stackexchange.com/questions/1105160/evaluate-derivative-of-lagrange-polynomials-at-construction-points
// 		*/
// 		var (
// 			XJ = R[j]
// 			x  = r
// 		)
// 		for l, XL := range R {
// 			if l != j {
// 				pp := 1.
// 				for m, XM := range R {
// 					if m != j && m != l {
// 						pp *= (x - XM) / (XJ - XM)
// 					}
// 				}
// 				p += pp / (XJ - XL)
// 			}
// 		}
// 	default:
// 		panic("invalid derivo")
// 	}
// 	return
// }
