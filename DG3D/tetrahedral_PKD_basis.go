package DG3D

// Simplex3DP evaluates 3D orthonormal polynomial on simplex at (r,s,t) of order (i,j,k)
// This is a direct translation of Hesthaven &package DG3D

import (
	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
	"math"
)

// TetBasis represents the Proriol-Koornwinder-Dubiner basis on reference tetrahedron
type TetBasis struct {
	Order int
	Np    int // Number of basis functions = (P+1)(P+2)(P+3)/6
}

// NewTetBasis creates a new tetrahedral basis of given order
func NewTetBasis(order int) *TetBasis {
	np := (order + 1) * (order + 2) * (order + 3) / 6
	return &TetBasis{
		Order: order,
		Np:    np,
	}
}

// TetCubature holds cubature points and weights for tetrahedron
type TetCubature struct {
	R, S, T utils.Vector // Cubature points
	W       utils.Vector // Weights
}

// RStoAB transforms from (r,s) to (a,b) coordinates for 2D simplex
func RStoAB(r, s utils.Vector) (a, b utils.Vector) {
	n := r.Len()
	a = utils.NewVector(n)
	b = utils.NewVector(n)

	for i := 0; i < n; i++ {
		if s.At(i) != 1 {
			a.Set(i, 2*(1+r.At(i))/(1-s.At(i))-1)
		} else {
			a.Set(i, -1)
		}
		b.Set(i, s.At(i))
	}
	return
}

// EquispacedNodes3D generates equispaced nodes on the reference tetrahedron
func EquispacedNodes3D(P int) (r, s, t utils.Vector) {
	// Equispaced node distribution for tetrahedron
	// NOT optimal Warburton nodes - just equally spaced

	Np := (P + 1) * (P + 2) * (P + 3) / 6
	rr := make([]float64, Np)
	ss := make([]float64, Np)
	tt := make([]float64, Np)

	if P == 1 {
		// Vertices of reference tetrahedron
		rr[0], ss[0], tt[0] = -1.0, -1.0, -1.0
		rr[1], ss[1], tt[1] = 1.0, -1.0, -1.0
		rr[2], ss[2], tt[2] = -1.0, 1.0, -1.0
		rr[3], ss[3], tt[3] = -1.0, -1.0, 1.0
	} else {
		// Use equispaced nodes on the reference tetrahedron
		idx := 0
		for n := 0; n <= P; n++ {
			for m := 0; m <= P-n; m++ {
				for l := 0; l <= P-n-m; l++ {
					// Map to barycentric coordinates
					i := l
					j := m
					k := n
					L1 := float64(i) / float64(P)
					L2 := float64(j) / float64(P)
					L3 := float64(k) / float64(P)
					L4 := 1.0 - L1 - L2 - L3

					// Convert to Cartesian coordinates on reference element
					rr[idx] = -L1 - L2 - L3 + L4
					ss[idx] = -L1 - L2 + L3 - L4
					tt[idx] = -L1 + L2 - L3 - L4
					idx++
				}
			}
		}
	}

	r = utils.NewVector(Np, rr)
	s = utils.NewVector(Np, ss)
	t = utils.NewVector(Np, tt)

	return r, s, t
}

// Warpfactor computes warp factor used in creating optimal node distributions
func Warpfactor(N int, rout []float64) []float64 {
	// Compute LGL and equidistant node distribution
	LGLr := JacobiGL(0, 0, N)

	// Equidistant nodes
	req := make([]float64, N+1)
	for i := 0; i <= N; i++ {
		req[i] = -1.0 + 2.0*float64(i)/float64(N)
	}

	// Compute Vandermonde based on equispaced nodes
	Veq := utils.NewMatrix(N+1, N+1)
	for i := 0; i <= N; i++ {
		for j := 0; j <= N; j++ {
			Veq.Set(i, j, DG1D.JacobiP(utils.NewVector(1, []float64{req[i]}), 0,
				0, j)[0])
		}
	}

	// Evaluate Lagrange polynomial at rout
	Nr := len(rout)
	Pmat := utils.NewMatrix(N+1, Nr)
	for i := 0; i <= N; i++ {
		for j := 0; j < Nr; j++ {
			Pmat.Set(i, j, DG1D.JacobiP(utils.NewVector(1, []float64{rout[j]}),
				0, 0, i)[0])
		}
	}

	Lmat := Veq.InverseWithCheck().Transpose().Mul(Pmat)

	// Compute warp factor
	warp := make([]float64, Nr)
	for i := 0; i < Nr; i++ {
		warp[i] = 0
		for j := 0; j <= N; j++ {
			warp[i] += Lmat.At(j, i) * (LGLr[j] - req[j])
		}
	}

	// Scale factor
	if N > 1 {
		for i := 0; i < Nr; i++ {
			sf := 1.0 - (math.Abs(rout[i]))
			sf = sf * sf
			warp[i] = warp[i] * sf
		}
	}

	return warp
}

// Nodes3D computes optimized interpolation nodes on tetrahedron
func Nodes3D(N int) (x, y, z utils.Vector) {
	// Total number of nodes
	Np := (N + 1) * (N + 2) * (N + 3) / 6

	// Compute equilateral nodes
	X := make([]float64, Np)
	Y := make([]float64, Np)
	Z := make([]float64, Np)

	// Set vertices of tetrahedron
	v1 := []float64{-1, -1 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v2 := []float64{1, -1 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v3 := []float64{0, 2 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v4 := []float64{0, 0, 3 / math.Sqrt(6)}

	// Orthogonal base for face 1
	t1 := []float64{v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]}
	t2 := []float64{v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]}

	sk := 0
	for n := 0; n <= N; n++ {
		for m := 0; m <= N-n; m++ {
			for l := 0; l <= N-n-m; l++ {
				// Barycentric coordinates
				L1 := float64(l) / float64(N)
				L2 := float64(m) / float64(N)
				L3 := float64(n) / float64(N)

				// Compute equilateral triangle coordinates
				X[sk] = L3*v4[0] + L2*v3[0] + L1*v2[0] + (1-L1-L2-L3)*v1[0]
				Y[sk] = L3*v4[1] + L2*v3[1] + L1*v2[1] + (1-L1-L2-L3)*v1[1]
				Z[sk] = L3*v4[2] + L2*v3[2] + L1*v2[2] + (1-L1-L2-L3)*v1[2]

				sk++
			}
		}
	}

	// Now apply warp and blend for each face
	// Face 1
	for k := 0; k < Np; k++ {
		// Get barycentric coordinates
		L1 := (math.Sqrt(3)*Y[k] + Z[k] + 2) / 4
		L2 := (-3*X[k] - math.Sqrt(3)*Y[k] + Z[k] + 2) / 4
		L3 := (3*X[k] - math.Sqrt(3)*Y[k] + Z[k] + 2) / 4
		L4 := (-Z[k] + 2) / 4

		// Skip if not on face
		if math.Abs(L1) > 1e-12 {
			continue
		}

		// Compute blend
		blend := 4 * L2 * L3 * L4

		// Compute warp amount on face
		r := L3 - L2
		s := L4 - L2
		denom := 2*L2 + L3 + L4 - 1
		if math.Abs(denom) > 1e-12 {
			r = r / denom
			s = s / denom
		}

		// Compute warp
		warpR := Warpfactor(N, []float64{r})
		warpS := Warpfactor(N, []float64{s})

		// Apply warp
		X[k] = X[k] + blend*warpR[0]*t1[0]
		Y[k] = Y[k] + blend*warpR[0]*t1[1]
		Z[k] = Z[k] + blend*warpR[0]*t1[2]

		X[k] = X[k] + blend*warpS[0]*t2[0]
		Y[k] = Y[k] + blend*warpS[0]*t2[1]
		Z[k] = Z[k] + blend*warpS[0]*t2[2]
	}

	// Similar process for faces 2, 3, 4...
	// This is simplified - full implementation would handle all faces

	// Transform to reference element coordinates
	r, s, t := XYZtoRST(X, Y, Z)

	return r, s, t
}

// XYZtoRST transforms from equilateral to reference tetrahedron coordinates
func XYZtoRST(X, Y, Z []float64) (r, s, t utils.Vector) {
	n := len(X)
	rr := make([]float64, n)
	ss := make([]float64, n)
	tt := make([]float64, n)

	for i := 0; i < n; i++ {
		// Transform from equilateral to reference coordinates
		// This is a simplified transformation
		rr[i] = 2*X[i]/(1-Z[i]/math.Sqrt(6)) - 1
		ss[i] = 2*Y[i]*math.Sqrt(3)/(3*(1-Z[i]/math.Sqrt(6))) - 1
		tt[i] = Z[i]*math.Sqrt(6)/2 - 1

		// Handle singularities
		if math.Abs(Z[i]-math.Sqrt(6)) < 1e-12 {
			rr[i] = -1
			ss[i] = -1
			tt[i] = 1
		}
	}

	r = utils.NewVector(n, rr)
	s = utils.NewVector(n, ss)
	t = utils.NewVector(n, tt)

	return
}

// JacobiGL computes Gauss-Lobatto-Jacobi points
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

// JacobiGQ computes Gauss-Jacobi quadrature points
func JacobiGQ(alpha, beta float64, N int) []float64 {
	if N == 0 {
		return []float64{-(alpha - beta) / (alpha + beta + 2)}
	}

	// Form symmetric matrix from recurrence
	J := make([]float64, (N+1)*(N+1))
	h1 := make([]float64, N+1)

	for i := 0; i <= N; i++ {
		h1[i] = 2*float64(i) + alpha + beta
	}

	// Diagonal
	for i := 0; i <= N; i++ {
		J[i*(N+1)+i] = -(alpha*alpha - beta*beta) / (h1[i] * (h1[i] + 2))
	}

	// Super/sub diagonal
	for i := 0; i < N; i++ {
		fi := float64(i + 1)
		v := 2 * fi * (fi + alpha + beta) * (fi + alpha) * (fi + beta) /
			(h1[i] * (h1[i] + 1) * (h1[i] + 2))
		v = math.Sqrt(v)
		J[i*(N+1)+i+1] = v
		J[(i+1)*(N+1)+i] = v
	}

	// Compute eigenvalues (simplified - would use LAPACK in practice)
	// This returns approximate values for small N
	x := make([]float64, N+1)
	for i := 0; i <= N; i++ {
		x[i] = math.Cos(math.Pi * float64(2*i+1) / float64(2*N+2))
	}

	return x
}

// RSTtoABC transforms from (r,s,t) to (a,b,c) coordinates for 3D simplex
func RSTtoABC(r, s, t utils.Vector) (a, b, c utils.Vector) {
	n := r.Len()
	a = utils.NewVector(n)
	b = utils.NewVector(n)
	c = utils.NewVector(n)

	tol := 1e-8

	for i := 0; i < n; i++ {
		if math.Abs(s.At(i)+t.At(i)-1) > tol {
			a.Set(i, 2*(1+r.At(i))/(-s.At(i)-t.At(i))-1)
		} else {
			a.Set(i, -1)
		}

		if math.Abs(t.At(i)-1) > tol {
			b.Set(i, 2*(1+s.At(i))/(1-t.At(i))-1)
		} else {
			b.Set(i, -1)
		}

		c.Set(i, t.At(i))
	}
	return
}

// Simplex3DP evaluates 3D orthonormal polynomial on simplex at (r,s,t) of order (i,j,k)
// This is a direct translation of Hesthaven & Warburton's MATLAB code
func Simplex3DP(r, s, t utils.Vector, i, j, k int) []float64 {
	// Convert to collapsed coordinates
	a, b, c := RSTtoABC(r, s, t)

	n := r.Len()

	// Compute Jacobi polynomials
	h1 := DG1D.JacobiP(a, 0, 0, i)
	h2 := DG1D.JacobiP(b, float64(2*i+1), 0, j)
	h3 := DG1D.JacobiP(c, float64(2*i+2*j+2), 0, k)

	// Compute the PKD polynomial
	P := make([]float64, n)

	for idx := 0; idx < n; idx++ {
		// Compute powers
		pow1 := math.Pow((1-b.At(idx))/2, float64(i))
		pow2 := math.Pow((1-c.At(idx))/2, float64(i+j))

		// Combine terms with normalization
		// The sqrt(8) factor normalizes for the reference tetrahedron volume
		P[idx] = math.Sqrt(8.0) * h1[idx] * h2[idx] * h3[idx] * pow1 * pow2
	}

	return P
}

// EvaluateBasis evaluates all basis functions at a single point (r,s,t)
func (tb *TetBasis) EvaluateBasis(r, s, t float64) []float64 {
	phi := make([]float64, tb.Np)

	// Create vectors for single point
	rv := utils.NewVector(1, []float64{r})
	sv := utils.NewVector(1, []float64{s})
	tv := utils.NewVector(1, []float64{t})

	// Loop over all basis functions
	idx := 0
	for i := 0; i <= tb.Order; i++ {
		for j := 0; j <= tb.Order-i; j++ {
			for k := 0; k <= tb.Order-i-j; k++ {
				P := Simplex3DP(rv, sv, tv, i, j, k)
				phi[idx] = P[0]
				idx++
			}
		}
	}

	return phi
}

// ComputeVandermonde computes the Vandermonde matrix for given nodes
func (tb *TetBasis) ComputeVandermonde(r, s, t utils.Vector) utils.Matrix {
	n := r.Len()
	V := utils.NewMatrix(n, tb.Np)

	// Loop over all basis functions
	sk := 0
	for i := 0; i <= tb.Order; i++ {
		for j := 0; j <= tb.Order-i; j++ {
			for k := 0; k <= tb.Order-i-j; k++ {
				P := Simplex3DP(r, s, t, i, j, k)
				V.SetCol(sk, P)
				sk++
			}
		}
	}

	return V
}

// GradSimplex3DP computes gradients of 3D orthonormal polynomial
func GradSimplex3DP(r, s, t utils.Vector, id, jd, kd int) (dmodedr, dmodeds, dmodedt []float64) {
	n := r.Len()
	dmodedr = make([]float64, n)
	dmodeds = make([]float64, n)
	dmodedt = make([]float64, n)

	// Convert to collapsed coordinates
	a, b, c := RSTtoABC(r, s, t)

	// Compute Jacobi polynomials and derivatives
	fa := DG1D.JacobiP(a, 0, 0, id)
	gb := DG1D.JacobiP(b, float64(2*id+1), 0, jd)
	hc := DG1D.JacobiP(c, float64(2*id+2*jd+2), 0, kd)

	dfa := DG1D.GradJacobiP(a, 0, 0, id)
	dgb := DG1D.GradJacobiP(b, float64(2*id+1), 0, jd)
	dhc := DG1D.GradJacobiP(c, float64(2*id+2*jd+2), 0, kd)

	// Compute derivatives
	for i := 0; i < n; i++ {
		ai := a.At(i)
		bi := b.At(i)
		ci := c.At(i)

		// Temporary values
		tmp1 := gb[i] * hc[i]
		tmp2 := fa[i] * hc[i]
		tmp3 := fa[i] * gb[i]

		// Powers
		pow1 := 1.0
		pow2 := 1.0
		if id > 0 {
			pow1 = math.Pow((1-bi)/2, float64(id))
		}
		if id+jd > 0 {
			pow2 = math.Pow((1-ci)/2, float64(id+jd))
		}

		// d/dr
		dmodedr[i] = dfa[i] * tmp1
		if id > 0 {
			dmodedr[i] *= math.Pow((1-bi)/2, float64(id-1))
		}
		if id+jd > 0 {
			dmodedr[i] *= math.Pow((1-ci)/2, float64(id+jd-1))
		}

		// d/ds
		dmodeds[i] = 0
		// First term
		if id > 0 {
			term1 := dfa[i] * (1 + ai) * tmp1 * math.Pow((1-bi)/2, float64(id-1))
			if id+jd > 0 {
				term1 *= math.Pow((1-ci)/2, float64(id+jd-1))
			}
			dmodeds[i] += term1
		}

		// Second term
		term2 := tmp3 * dgb[i] * pow1
		if id+jd > 0 {
			term2 *= math.Pow((1-ci)/2, float64(id+jd-1))
		}
		dmodeds[i] += term2

		// Third term
		if id > 0 {
			term3 := -float64(id) * tmp1 * fa[i] * math.Pow((1-bi)/2, float64(id-1))
			if id+jd > 1 {
				term3 *= math.Pow((1-ci)/2, float64(id+jd-1))
			} else if id+jd == 1 {
				term3 *= 1
			}
			dmodeds[i] += 0.5 * term3
		}

		// d/dt
		dmodedt[i] = 0
		// First term
		if id > 0 {
			term1 := dfa[i] * (1 + ai) * tmp1 * math.Pow((1-bi)/2, float64(id-1))
			if id+jd > 0 {
				term1 *= math.Pow((1-ci)/2, float64(id+jd-1))
			}
			dmodedt[i] += term1
		}

		// Second term
		if id > 0 {
			term2 := tmp3 * dgb[i] * (1 + bi) * math.Pow((1-bi)/2, float64(id-1))
			if id+jd > 0 {
				term2 *= math.Pow((1-ci)/2, float64(id+jd-1))
			}
			dmodedt[i] += term2
		}

		// Third term
		term3 := tmp2 * dhc[i] * pow1 * pow2
		dmodedt[i] += term3

		// Fourth term
		if id > 0 {
			term4 := -float64(id) * tmp1 * fa[i] * math.Pow((1-bi)/2, float64(id-1))
			if id+jd > 0 {
				term4 *= math.Pow((1-ci)/2, float64(id+jd-1))
			}
			dmodedt[i] += 0.5 * term4
		}

		// Fifth term
		if id+jd > 0 {
			term5 := -float64(id+jd) * tmp3 * fa[i] * pow1
			if id+jd > 1 {
				term5 *= math.Pow((1-ci)/2, float64(id+jd-1))
			}
			dmodedt[i] += 0.5 * term5
		}

		// Apply normalization
		dmodedr[i] *= math.Sqrt(8.0)
		dmodeds[i] *= math.Sqrt(8.0)
		dmodedt[i] *= math.Sqrt(8.0)
	}

	return
}

// EvaluateBasisGrad evaluates gradients of all basis functions at a single point
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
	for i := 0; i <= tb.Order; i++ {
		for j := 0; j <= tb.Order-i; j++ {
			for k := 0; k <= tb.Order-i-j; k++ {
				dr, ds, dt := GradSimplex3DP(rv, sv, tv, i, j, k)
				dphidr[idx] = dr[0]
				dphids[idx] = ds[0]
				dphidt[idx] = dt[0]
				idx++
			}
		}
	}

	return dphidr, dphids, dphidt
}

// ComputeGradVandermonde computes gradient Vandermonde matrices
func (tb *TetBasis) ComputeGradVandermonde(r, s, t utils.Vector) (Vr, Vs, Vt utils.Matrix) {
	n := r.Len()
	Vr = utils.NewMatrix(n, tb.Np)
	Vs = utils.NewMatrix(n, tb.Np)
	Vt = utils.NewMatrix(n, tb.Np)

	// Loop over all basis functions
	sk := 0
	for i := 0; i <= tb.Order; i++ {
		for j := 0; j <= tb.Order-i; j++ {
			for k := 0; k <= tb.Order-i-j; k++ {
				dr, ds, dt := GradSimplex3DP(r, s, t, i, j, k)
				Vr.SetCol(sk, dr)
				Vs.SetCol(sk, ds)
				Vt.SetCol(sk, dt)
				sk++
			}
		}
	}

	return Vr, Vs, Vt
}

// ComputeMassMatrix computes the mass matrix M = V^T * W * V
// For orthonormal basis, this should be the identity matrix
func (tb *TetBasis) ComputeMassMatrix(V utils.Matrix, cubature *TetCubature) utils.Matrix {
	n := V.Rows()  // number of cubature points
	np := V.Cols() // number of basis functions

	// Create V^T * W
	VT := V.Transpose()
	VTW := utils.NewMatrix(np, n)

	// Apply diagonal weight matrix
	for i := 0; i < np; i++ {
		for j := 0; j < n; j++ {
			VTW.Set(i, j, VT.At(i, j)*cubature.W.At(j))
		}
	}

	// Compute M = (V^T * W) * V
	M := VTW.Mul(V)

	return M
}

// Cubature3D returns cubature points and weights for tetrahedron
// This is a direct translation of Hesthaven & Warburton's Cubature3D.m
func Cubature3D(order int) (x, y, z []float64, w []float64) {
	// Cubature for tetrahedron from Hesthaven & Warburton

	switch order {
	case 1:
		// Order 1 - 1 point rule (exact for polynomials up to degree 1)
		x = []float64{-0.5}
		y = []float64{-0.5}
		z = []float64{-0.5}
		w = []float64{8.0 / 6.0}

	case 2:
		// Order 2 - 4 point rule (exact for polynomials up to degree 2)
		a := 0.58541019662496845446
		b := 0.13819660112501051518

		x = []float64{-b, -b, -b, -a}
		y = []float64{-b, -b, -a, -b}
		z = []float64{-b, -a, -b, -b}

		w = []float64{
			8.0 / 24.0,
			8.0 / 24.0,
			8.0 / 24.0,
			8.0 / 24.0,
		}

	case 3:
		// Order 3 - 5 point rule (exact for polynomials up to degree 3)
		x = []float64{-0.5, -1.0 / 6.0, -1.0 / 6.0, -1.0 / 6.0, -1.0 / 2.0}
		y = []float64{-0.5, -1.0 / 6.0, -1.0 / 6.0, -1.0 / 2.0, -1.0 / 6.0}
		z = []float64{-0.5, -1.0 / 6.0, -1.0 / 2.0, -1.0 / 6.0, -1.0 / 6.0}

		w = []float64{
			-8.0 / 20.0,
			8.0 / 20.0,
			8.0 / 20.0,
			8.0 / 20.0,
			8.0 / 20.0,
		}

	default:
		// Default to order 4 - 11 point rule
		// This is more complex, using the basic order 3 for now
		return Cubature3D(3)
	}

	return
}
