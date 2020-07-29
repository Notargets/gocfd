package DG1D

import (
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
	el.X.SetReadOnly("X")
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
		//x[i] = xint.AtVec(iter)
		x[i] = dataXint[iter]
		iter++
	}
	R = utils.NewVector(len(x), x)
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
		x_bnew := utils.NewVector(r.Len()).Set(-bnew)
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

func GradVandermonde1D(r utils.Vector, N int) (Vr utils.Matrix) {
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
		LiData[j] *= LagrangePolyAtJ(r, el.R.Data(), j)
	}
	return
}

func LagrangePolyAtJ(t float64, R []float64, j int) (s float64) {
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
	valj := R[j]
	s = 1
	for i, val := range R {
		if i == j {
			continue
		}
		metric := (t - val) / (valj - val)
		s *= metric
	}
	return
}

func JacobiP2D(u, v float64, i, j int) (p float64) {
	/*
		This is a 2D normalized polynomial basis for the unit triangle with vertices:
					(0,0), (0,1), (1, 0)
		Arguments:
					(u,v) is the position within the triangle where the basis is evaluated
		    		(i,j) is the index within the basis, which has (K+1)(K+2)/2 terms
	*/
	a := utils.NewVector(1, []float64{2*u/(1-v) - 1})
	b := utils.NewVector(1, []float64{2*v - 1})
	h1 := JacobiP(a, 0, 0, i)
	h2 := JacobiP(b, float64(2*i+1), 0, i)
	p = math.Sqrt(2) * utils.POW(2*(1-v), j) * h1[0] * h2[0] / 4 // The divide by 4 is the Jacobian that converts from a different triangle
	return
}
