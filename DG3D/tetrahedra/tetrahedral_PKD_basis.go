package tetrahedra

import (
	"math"
	"sort"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
)

// TetBasis represents the Proriol-Koornwinder-Dubiner basis on reference tetrahedron
// Note: This basis is orthogonal with respect to the weighted inner product on the tetrahedron,
// but NOT orthonormal with respect to the standard L2 inner product
type TetBasis struct {
	N       int
	Np      int // Number of basis functions = (P+1)(P+2)(P+3)/6
	R, S, T utils.Vector
	V, VInv utils.Matrix // Vandermonde matrix and inverse
	// Mass and differentiation matrices
	M, MInv    utils.Matrix // Mass matrix, Inverse mass matrix
	Dr, Ds, Dt utils.Matrix // Differentiation matrices
	LIFT       utils.Matrix // Lift matrix for surface integrals
	Fmask      [][]int      // Face node indices
	Nfp        int          // Number of face points per face
}

func NewTetBasis(N int) (tb *TetBasis) {
	tb = &TetBasis{
		N:  N,
		Np: (N + 1) * (N + 2) * (N + 3) / 6,
	}
	tb.R, tb.S, tb.T = Nodes3D(N)
	tb.V = tb.ComputeVandermonde(tb.R, tb.S, tb.T)
	tb.VInv = tb.V.InverseWithCheck()
	tb.LIFT = tb.Lift3D()
	Vr, Vs, Vt := tb.ComputeGradVandermonde(tb.R, tb.S, tb.T)
	tb.Dr = Vr.Mul(tb.VInv)
	tb.Ds = Vs.Mul(tb.VInv)
	tb.Dt = Vt.Mul(tb.VInv)
	tb.M = tb.VInv.Transpose().Mul(tb.VInv)
	inverse, err := tb.M.Inverse()
	if err != nil {
		return nil
	}
	tb.MInv = inverse
	return
}

func JacobiGL(alpha, beta float64, N int) []float64 {
	if N == 0 {
		return []float64{0.0}
	}

	x := make([]float64, N+1)
	if N == 1 {
		x[0] = -1.0
		x[1] = 1.0
		return x
	}

	// Interior Gauss-Jacobi points
	xint := JacobiGQ(alpha+1, beta+1, N-1)

	x[0] = -1.0
	for i := 0; i < N-1; i++ {
		x[i+1] = xint[i]
	}
	x[N] = 1.0

	return x
}

func JacobiGQ(alpha, beta float64, N int) []float64 {
	if N == 0 {
		return []float64{-(alpha - beta) / (alpha + beta + 2)}
	}

	// Form symmetric matrix from recurrence
	J := utils.NewMatrix(N+1, N+1)
	h1 := make([]float64, N+1)

	for i := 0; i <= N; i++ {
		h1[i] = 2*float64(i) + alpha + beta
	}

	// Diagonal
	for i := 0; i <= N; i++ {
		J.Set(i, i, -(alpha*alpha-beta*beta)/(h1[i]*(h1[i]+2)))
	}

	// Super/sub diagonal
	for i := 0; i < N; i++ {
		fi := float64(i + 1)
		v := 2 * fi * (fi + alpha + beta) * (fi + alpha) * (fi + beta) /
			(h1[i] * (h1[i] + 1) * (h1[i] + 2))
		v = math.Sqrt(v)
		J.Set(i, i+1, v)
		J.Set(i+1, i, v)
	}

	// Compute eigenvalues using gonum's eigenvalue solver
	var eig mat.Eigen
	ok := eig.Factorize(J.M, mat.EigenLeft)
	if !ok {
		panic("eigenvalue decomposition failed")
	}

	// Extract eigenvalues (real parts)
	values := eig.Values(nil)
	x := make([]float64, N+1)
	for i := 0; i <= N; i++ {
		x[i] = real(values[i])
	}

	// Sort eigenvalues in ascending order
	sort.Float64s(x)

	return x
}

func Simplex3DP(r, s, t utils.Vector, i, j, k int) []float64 {
	// Convert to collapsed coordinates
	a, b, c := RSTtoABC(r, s, t)

	n := r.Len()

	// Compute Jacobi polynomials
	h1 := DG1D.JacobiP(a, 0, 0, i)
	h2 := DG1D.JacobiP(b, float64(2*i+1), 0, j)
	h3 := DG1D.JacobiP(c, float64(2*(i+j)+2), 0, k)

	// Compute the PKD polynomial
	P := make([]float64, n)

	for idx := 0; idx < n; idx++ {
		bi := b.At(idx)
		ci := c.At(idx)

		// Match C++ exactly: tv1 = 2.0*sqrt(2.0)*h1.dm(h2)
		tv1 := 2.0 * math.Sqrt(2.0) * h1[idx] * h2[idx]

		// Match C++ exactly: tv2 = pow(1.0-b,(double)i)
		tv2 := math.Pow(1.0-bi, float64(i))

		// Match C++ exactly: tv3 = h3.dm(pow(1.0-c,(double)(i+j)))
		tv3 := h3[idx] * math.Pow(1.0-ci, float64(i+j))

		// Match C++ exactly: (*P) = tv1.dm( tv2.dm(tv3) )
		P[idx] = tv1 * tv2 * tv3
	}

	return P
}

func (tb *TetBasis) BuildFmask3D() (fmask [][]int) {
	var (
		Np      = tb.Np
		r, s, t = tb.R, tb.S, tb.T
	)
	fmask = make([][]int, 4)

	// Face tolerances
	NODETOL := 1e-10

	// Temporary storage for face nodes
	for f := 0; f < 4; f++ {
		fmask[f] = make([]int, 0)
	}

	// Face 1: t = -1 (matching C++: abs(1+t) < NODETOL)
	for i := 0; i < Np; i++ {
		if math.Abs(1.0+t.At(i)) < NODETOL {
			fmask[0] = append(fmask[0], i)
		}
	}

	// Face 2: s = -1 (matching C++: abs(1+s) < NODETOL)
	for i := 0; i < Np; i++ {
		if math.Abs(1.0+s.At(i)) < NODETOL {
			fmask[1] = append(fmask[1], i)
		}
	}

	// Face 3: r+s+t = -1 (matching C++: abs(1+r+s+t) < NODETOL)
	for i := 0; i < Np; i++ {
		if math.Abs(1.0+r.At(i)+s.At(i)+t.At(i)) < NODETOL {
			fmask[2] = append(fmask[2], i)
		}
	}

	// Face 4: r = -1 (matching C++: abs(1+r) < NODETOL)
	for i := 0; i < Np; i++ {
		if math.Abs(1.0+r.At(i)) < NODETOL {
			fmask[3] = append(fmask[3], i)
		}
	}

	return
}

func (tb *TetBasis) Lift3D() (LIFT utils.Matrix) {
	tb.Fmask = tb.BuildFmask3D()
	tb.Nfp = (tb.N + 1) * (tb.N + 2) / 2 // Number of face points
	var (
		N       = tb.N
		Np      = tb.Np
		V       = tb.V
		R, S, T = tb.R, tb.S, tb.T
		fmask   = tb.Fmask
	)
	Nfaces := 4

	// Create face mass matrix
	Emat := utils.NewMatrix(Np, Nfaces*tb.Nfp)

	// Build the Emat matrix
	for face := 0; face < Nfaces; face++ {
		// Extract face coordinates according to C++ mapping
		// Note: C++ uses 1-based face indexing, we use 0-based
		var faceR, faceS utils.Vector

		switch face {
		case 0: // C++ face 1: use R,S coordinates
			faceR = R.SubsetIndex(fmask[face])
			faceS = S.SubsetIndex(fmask[face])
		case 1: // C++ face 2: use R,T coordinates
			faceR = R.SubsetIndex(fmask[face])
			faceS = T.SubsetIndex(fmask[face])
		case 2: // C++ face 3: use S,T coordinates
			faceR = S.SubsetIndex(fmask[face])
			faceS = T.SubsetIndex(fmask[face])
		case 3: // C++ face 4: use S,T coordinates
			faceR = S.SubsetIndex(fmask[face])
			faceS = T.SubsetIndex(fmask[face])
		}

		// Compute face Vandermonde and mass matrix
		faceV := DG2D.Vandermonde2D(N, faceR, faceS)
		massFace := faceV.Mul(faceV.Transpose()).InverseWithCheck()

		// Fill Emat
		// C++: JJ.reset((face-1)*Nfp+1, face*Nfp) creates 1-based indices
		// Go: we need 0-based indices: face*Nfp to (face+1)*Nfp-1
		for i := 0; i < tb.Nfp; i++ {
			for j := 0; j < tb.Nfp; j++ {
				// C++ assigns massFace to Emat(idr, JJ)
				// idr = Fmask(All,face) = Fmask[face]
				// JJ columns are face*Nfp+j (0-based)
				Emat.Set(fmask[face][i], face*tb.Nfp+j, massFace.At(i, j))
			}
		}
	}

	// Compute LIFT = V*V'*Emat
	LIFT = V.Mul(V.Transpose()).Mul(Emat)

	return
}

func (tb *TetBasis) EvaluateBasis(r, s, t float64) (phi []float64) {
	phi = make([]float64, tb.Np)

	// Create vectors for single evaluation point
	rv := utils.NewVector(1, []float64{r})
	sv := utils.NewVector(1, []float64{s})
	tv := utils.NewVector(1, []float64{t})

	// Loop over all basis functions
	idx := 0
	for i := 0; i <= tb.N; i++ {
		for j := 0; j <= tb.N-i; j++ {
			for k := 0; k <= tb.N-i-j; k++ {
				P := Simplex3DP(rv, sv, tv, i, j, k)
				phi[idx] = P[0]
				idx++
			}
		}
	}
	return
}

func (tb *TetBasis) ComputeVandermonde(r, s, t utils.Vector) (V utils.Matrix) {
	n := r.Len()
	V = utils.NewMatrix(n, tb.Np)

	// Loop over all basis functions
	sk := 0
	for i := 0; i <= tb.N; i++ {
		for j := 0; j <= tb.N-i; j++ {
			for k := 0; k <= tb.N-i-j; k++ {
				P := Simplex3DP(r, s, t, i, j, k)
				V.SetCol(sk, P)
				sk++
			}
		}
	}

	return V
}

func (tb *TetBasis) GradSimplex3DP(r, s, t utils.Vector, id, jd,
	kd int) (dmodedr, dmodeds, dmodedt []float64) {
	n := r.Len()
	dmodedr = make([]float64, n)
	dmodeds = make([]float64, n)
	dmodedt = make([]float64, n)

	// Convert to collapsed coordinates
	a, b, c := RSTtoABC(r, s, t)

	// Compute Jacobi polynomials
	fa := DG1D.JacobiP(a, 0, 0, id)
	gb := DG1D.JacobiP(b, float64(2*id+1), 0, jd)
	hc := DG1D.JacobiP(c, float64(2*id+2*jd+2), 0, kd)

	// Compute derivatives
	dfa := DG1D.GradJacobiP(a, 0, 0, id)
	dgb := DG1D.GradJacobiP(b, float64(2*id+1), 0, jd)
	dhc := DG1D.GradJacobiP(c, float64(2*id+2*jd+2), 0, kd)

	// Compute each derivative component
	for i := 0; i < n; i++ {
		ai := a.At(i)
		bi := b.At(i)
		ci := c.At(i)

		// r-derivative
		V3Dr := dfa[i] * (gb[i] * hc[i])
		if id > 0 {
			V3Dr = V3Dr * math.Pow(0.5*(1-bi), float64(id-1))
		}
		if id+jd > 0 {
			V3Dr = V3Dr * math.Pow(0.5*(1-ci), float64(id+jd-1))
		}
		dmodedr[i] = V3Dr

		// s-derivative
		V3Ds := 0.5 * (1 + ai) * V3Dr
		tmp := dgb[i] * math.Pow(0.5*(1-bi), float64(id))
		if id > 0 {
			tmp = tmp + (-0.5*float64(id))*(gb[i]*math.Pow(0.5*(1-bi), float64(id-1)))
		}
		if id+jd > 0 {
			tmp = tmp * math.Pow(0.5*(1-ci), float64(id+jd-1))
		}
		tmp = fa[i] * (tmp * hc[i])
		V3Ds = V3Ds + tmp
		dmodeds[i] = V3Ds

		// t-derivative
		V3Dt := 0.5*(1+ai)*V3Dr + 0.5*(1+bi)*tmp
		tmp2 := dhc[i] * math.Pow(0.5*(1-ci), float64(id+jd))
		if id+jd > 0 {
			tmp2 = tmp2 - 0.5*float64(id+jd)*(hc[i]*math.Pow(0.5*(1-ci), float64(id+jd-1)))
		}
		tmp2 = fa[i] * (gb[i] * tmp2)
		tmp2 = tmp2 * math.Pow(0.5*(1-bi), float64(id))
		V3Dt = V3Dt + tmp2
		dmodedt[i] = V3Dt

		// normalize
		normFactor := math.Pow(2, float64(2*id+jd)+1.5)
		dmodedr[i] = dmodedr[i] * normFactor
		dmodeds[i] = dmodeds[i] * normFactor
		dmodedt[i] = dmodedt[i] * normFactor
	}

	return
}

func (tb *TetBasis) EvaluateBasisGrad(r, s, t float64) (dphidr, dphids, dphidt []float64) {
	dphidr = make([]float64, tb.Np)
	dphids = make([]float64, tb.Np)
	dphidt = make([]float64, tb.Np)

	// Create vectors for single point
	rv := utils.NewVector(1, []float64{r})
	sv := utils.NewVector(1, []float64{s})
	tv := utils.NewVector(1, []float64{t})

	// Loop over all basis functions
	idx := 0
	for i := 0; i <= tb.N; i++ {
		for j := 0; j <= tb.N-i; j++ {
			for k := 0; k <= tb.N-i-j; k++ {
				dr, ds, dt := tb.GradSimplex3DP(rv, sv, tv, i, j, k)
				dphidr[idx] = dr[0]
				dphids[idx] = ds[0]
				dphidt[idx] = dt[0]
				idx++
			}
		}
	}

	return dphidr, dphids, dphidt
}

func (tb *TetBasis) ComputeGradVandermonde(r, s, t utils.Vector) (Vr, Vs, Vt utils.Matrix) {
	n := r.Len()
	Vr = utils.NewMatrix(n, tb.Np)
	Vs = utils.NewMatrix(n, tb.Np)
	Vt = utils.NewMatrix(n, tb.Np)

	// Loop over all basis functions
	sk := 0
	for i := 0; i <= tb.N; i++ {
		for j := 0; j <= tb.N-i; j++ {
			for k := 0; k <= tb.N-i-j; k++ {
				dr, ds, dt := tb.GradSimplex3DP(r, s, t, i, j, k)
				Vr.SetCol(sk, dr)
				Vs.SetCol(sk, ds)
				Vt.SetCol(sk, dt)
				sk++
			}
		}
	}

	return Vr, Vs, Vt
}

func evalshift(p int, pval float64, L1, L2, L3 []float64) (dx, dy []float64) {
	n := len(L1)
	dx = make([]float64, n)
	dy = make([]float64, n)

	// 1) compute Gauss-Lobatto-Legendre node distribution
	gaussX := JacobiGL(0, 0, p)
	// Negate the values (C++: gaussX = -JacobiGL(0,0,p))
	for i := range gaussX {
		gaussX[i] = -gaussX[i]
	}

	// 3) compute blending function at each node for each edge
	blend1 := make([]float64, n)
	blend2 := make([]float64, n)
	blend3 := make([]float64, n)

	for i := 0; i < n; i++ {
		blend1[i] = L2[i] * L3[i]
		blend2[i] = L1[i] * L3[i]
		blend3[i] = L1[i] * L2[i]
	}

	// 4) amount of warp for each node, for each edge
	tv1 := make([]float64, n)
	tv2 := make([]float64, n)
	tv3 := make([]float64, n)

	for i := 0; i < n; i++ {
		tv1[i] = L3[i] - L2[i]
		tv2[i] = L1[i] - L3[i]
		tv3[i] = L2[i] - L1[i]
	}

	warpfactor1 := evalwarp(p, gaussX, tv1)
	warpfactor2 := evalwarp(p, gaussX, tv2)
	warpfactor3 := evalwarp(p, gaussX, tv3)

	// Scale by 4.0
	for i := 0; i < n; i++ {
		warpfactor1[i] *= 4.0
		warpfactor2[i] *= 4.0
		warpfactor3[i] *= 4.0
	}

	// 5) combine blend & warp
	warp1 := make([]float64, n)
	warp2 := make([]float64, n)
	warp3 := make([]float64, n)

	for i := 0; i < n; i++ {
		warp1[i] = blend1[i] * warpfactor1[i] * (1.0 + pval*pval*L1[i]*L1[i])
		warp2[i] = blend2[i] * warpfactor2[i] * (1.0 + pval*pval*L2[i]*L2[i])
		warp3[i] = blend3[i] * warpfactor3[i] * (1.0 + pval*pval*L3[i]*L3[i])
	}

	// 6) evaluate shift in equilateral triangle
	TWOPI := 2.0 * math.Pi
	FOURPI := 4.0 * math.Pi

	for i := 0; i < n; i++ {
		dx[i] = 1.0*warp1[i] + math.Cos(TWOPI/3.0)*warp2[i] + math.Cos(FOURPI/3.0)*warp3[i]
		dy[i] = 0.0*warp1[i] + math.Sin(TWOPI/3.0)*warp2[i] + math.Sin(FOURPI/3.0)*warp3[i]
	}

	return dx, dy
}

func evalwarp(p int, gaussX []float64, xnodes []float64) []float64 {
	n := len(xnodes)
	warp := make([]float64, n)

	// Create equidistant nodes
	xeq := make([]float64, p+1)
	for i := 0; i <= p; i++ {
		xeq[i] = -1.0 + 2.0*float64(i)/float64(p)
	}

	// For each evaluation point
	for i := 0; i < n; i++ {
		x := xnodes[i]

		// Compute Lagrange interpolation from equidistant to GLL nodes
		warpval := 0.0
		for j := 0; j <= p; j++ {
			// Lagrange basis at xeq[j] evaluated at x
			lagrange := 1.0
			for k := 0; k <= p; k++ {
				if k != j {
					lagrange *= (x - xeq[k]) / (xeq[j] - xeq[k])
				}
			}
			warpval += lagrange * (gaussX[j] - xeq[j])
		}

		// Scale factor
		sf := 1.0 - x*x
		warp[i] = warpval * sf
	}

	return warp
}
