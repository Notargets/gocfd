package DG3D

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

// FaceCubature holds cubature points and weights for triangular faces
type FaceCubature struct {
	R, S utils.Vector // Cubature points on 2D triangle
	W    utils.Vector // Weights
}

// GeometricFactors stores metric terms for physical elements
type GeometricFactors struct {
	X, Y, Z    utils.Matrix // Physical coordinates at each node
	Rx, Ry, Rz utils.Matrix // dr/dx, dr/dy, dr/dz
	Sx, Sy, Sz utils.Matrix // ds/dx, ds/dy, ds/dz
	Tx, Ty, Tz utils.Matrix // dt/dx, dt/dy, dt/dz
	J          utils.Matrix // Jacobian determinant
}

// FaceGeometricFactors stores metric terms for element faces
type FaceGeometricFactors struct {
	Nx, Ny, Nz utils.Matrix // Outward normal components
	SJ         utils.Matrix // Surface Jacobian (area scaling)
	Fscale     utils.Matrix // Face integration scaling = sJ/J(face)
}

// Element3D represents a tetrahedral element
type Element3D struct {
	*TetBasis

	// Node information
	R, S, T utils.Vector // Reference coordinates
	X, Y, Z utils.Vector // Physical coordinates

	// Mass and differentiation matrices
	V          utils.Matrix // Vandermonde matrix
	Vr, Vs, Vt utils.Matrix // Gradient Vandermonde matrices
	Dr, Ds, Dt utils.Matrix // Differentiation matrices
	M          utils.Matrix // Mass matrix
	invM       utils.Matrix // Inverse mass matrix

	// Surface integration
	LIFT       utils.Matrix   // Lift matrix for surface integrals
	Fmask      [][]int        // Face node indices
	Fx, Fy, Fz []utils.Matrix // Physical coordinates on each face

	// Connectivity
	VmapM  []int // Maps face nodes to volume nodes
	VmapP  []int // Maps face nodes to neighbor volume nodes
	MapM   []int // Maps face nodes to unique face nodes
	MapP   []int // Maps to neighbor face nodes
	BCType []int // Boundary condition types per face

	// Geometric factors
	GeomFactors     *GeometricFactors
	FaceGeomFactors *FaceGeometricFactors
}

// ConnectivityArrays stores mesh connectivity information
type ConnectivityArrays struct {
	EToE [][]int // Element to element connectivity
	EToF [][]int // Element to face connectivity
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

func RSTtoABC(r, s, t utils.Vector) (a, b, c utils.Vector) {
	n := r.Len()
	a = utils.NewVector(n)
	b = utils.NewVector(n)
	c = utils.NewVector(n)

	tol := 1e-8

	for i := 0; i < n; i++ {
		if math.Abs(s.At(i)+t.At(i)) > tol { // Changed from s+t-1 to s+t
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

// EquispacedNodes3D generates equispaced nodes on the reference tetrahedron
func EquispacedNodes3D(P int) (r, s, t utils.Vector) {
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
			Veq.Set(i, j, DG1D.JacobiP(utils.NewVector(1, []float64{req[i]}), 0, 0, j)[0])
		}
	}

	// Evaluate Lagrange polynomial at rout
	Nr := len(rout)
	Pmat := utils.NewMatrix(N+1, Nr)
	for i := 0; i <= N; i++ {
		for j := 0; j < Nr; j++ {
			Pmat.Set(i, j, DG1D.JacobiP(utils.NewVector(1, []float64{rout[j]}), 0, 0, i)[0])
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
	// Choose optimized blending parameter
	alphastore := []float64{
		0.0000, 0.0000, 0.0000, 0.1002, 1.1332,
		1.5608, 1.3413, 1.2577, 1.1603, 1.10153,
		0.6080, 0.4523, 0.8856, 0.8717, 0.9655,
	}

	alpha := 1.0
	if N <= 14 { // 0-based indexing, so N=14 corresponds to p=15 in C++
		alpha = alphastore[N]
	}

	// Total number of nodes and tolerance
	Np := (N + 1) * (N + 2) * (N + 3) / 6
	tol := 1e-10
	sqrt3 := math.Sqrt(3.0)
	sqrt6 := math.Sqrt(6.0)

	// Create equidistributed nodes
	r, s, t := EquiNodes3D(N)

	// Compute barycentric coordinates from (r,s,t)
	L1 := make([]float64, Np)
	L2 := make([]float64, Np)
	L3 := make([]float64, Np)
	L4 := make([]float64, Np)

	for i := 0; i < Np; i++ {
		L1[i] = (1.0 + t.At(i)) / 2.0
		L2[i] = (1.0 + s.At(i)) / 2.0
		L3[i] = -(1.0 + r.At(i) + s.At(i) + t.At(i)) / 2.0
		L4[i] = (1.0 + r.At(i)) / 2.0
	}

	// Set vertices of tetrahedron
	v1 := []float64{-1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v2 := []float64{1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v3 := []float64{0.0, 2.0 / sqrt3, -1.0 / sqrt6}
	v4 := []float64{0.0, 0.0, 3.0 / sqrt6}

	// Form undeformed coordinates
	X := make([]float64, Np)
	Y := make([]float64, Np)
	Z := make([]float64, Np)

	for i := 0; i < Np; i++ {
		X[i] = L3[i]*v1[0] + L4[i]*v2[0] + L2[i]*v3[0] + L1[i]*v4[0]
		Y[i] = L3[i]*v1[1] + L4[i]*v2[1] + L2[i]*v3[1] + L1[i]*v4[1]
		Z[i] = L3[i]*v1[2] + L4[i]*v2[2] + L2[i]*v3[2] + L1[i]*v4[2]
	}

	// Shift coordinates for warping
	shiftX := make([]float64, Np)
	shiftY := make([]float64, Np)
	shiftZ := make([]float64, Np)

	// Orthogonal axis tangents on faces 1-4
	t1 := make([][]float64, 4)
	t2 := make([][]float64, 4)

	// Face 1
	t1[0] = []float64{v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]}
	t2[0] = []float64{v3[0] - 0.5*(v1[0]+v2[0]), v3[1] - 0.5*(v1[1]+v2[1]), v3[2] - 0.5*(v1[2]+v2[2])}

	// Face 2
	t1[1] = []float64{v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]}
	t2[1] = []float64{v4[0] - 0.5*(v1[0]+v2[0]), v4[1] - 0.5*(v1[1]+v2[1]), v4[2] - 0.5*(v1[2]+v2[2])}

	// Face 3
	t1[2] = []float64{v3[0] - v2[0], v3[1] - v2[1], v3[2] - v2[2]}
	t2[2] = []float64{v4[0] - 0.5*(v2[0]+v3[0]), v4[1] - 0.5*(v2[1]+v3[1]), v4[2] - 0.5*(v2[2]+v3[2])}

	// Face 4
	t1[3] = []float64{v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]}
	t2[3] = []float64{v4[0] - 0.5*(v1[0]+v3[0]), v4[1] - 0.5*(v1[1]+v3[1]), v4[2] - 0.5*(v1[2]+v3[2])}

	// Normalize tangents
	for n := 0; n < 4; n++ {
		norm1 := math.Sqrt(t1[n][0]*t1[n][0] + t1[n][1]*t1[n][1] + t1[n][2]*t1[n][2])
		norm2 := math.Sqrt(t2[n][0]*t2[n][0] + t2[n][1]*t2[n][1] + t2[n][2]*t2[n][2])
		for j := 0; j < 3; j++ {
			t1[n][j] /= norm1
			t2[n][j] /= norm2
		}
	}

	// Warp and blend for each face
	for face := 0; face < 4; face++ {
		var La, Lb, Lc, Ld []float64

		switch face {
		case 0: // Face 1
			La, Lb, Lc, Ld = L1, L2, L3, L4
		case 1: // Face 2
			La, Lb, Lc, Ld = L2, L1, L3, L4
		case 2: // Face 3
			La, Lb, Lc, Ld = L3, L1, L4, L2
		case 3: // Face 4
			La, Lb, Lc, Ld = L4, L1, L3, L2
		}

		// Compute warp tangential to face
		warp1, warp2 := WarpShiftFace3D(N, alpha, alpha, La, Lb, Lc, Ld)

		// Compute volume blending
		blend := make([]float64, Np)
		denom := make([]float64, Np)

		for i := 0; i < Np; i++ {
			blend[i] = Lb[i] * Lc[i] * Ld[i]
			denom[i] = (Lb[i] + 0.5*La[i]) * (Lc[i] + 0.5*La[i]) * (Ld[i] + 0.5*La[i])

			if denom[i] > tol {
				blend[i] = (1.0 + alpha*alpha*La[i]*La[i]) * blend[i] / denom[i]
			}
		}

		// Compute warp & blend
		for i := 0; i < Np; i++ {
			shiftX[i] += blend[i]*warp1[i]*t1[face][0] + blend[i]*warp2[i]*t2[face][0]
			shiftY[i] += blend[i]*warp1[i]*t1[face][1] + blend[i]*warp2[i]*t2[face][1]
			shiftZ[i] += blend[i]*warp1[i]*t1[face][2] + blend[i]*warp2[i]*t2[face][2]
		}

		// Fix face warp
		for i := 0; i < Np; i++ {
			if La[i] < tol {
				count := 0
				if Lb[i] > tol {
					count++
				}
				if Lc[i] > tol {
					count++
				}
				if Ld[i] > tol {
					count++
				}

				if count < 3 {
					shiftX[i] = warp1[i]*t1[face][0] + warp2[i]*t2[face][0]
					shiftY[i] = warp1[i]*t1[face][1] + warp2[i]*t2[face][1]
					shiftZ[i] = warp1[i]*t1[face][2] + warp2[i]*t2[face][2]
				}
			}
		}
	}

	// Apply shift
	for i := 0; i < Np; i++ {
		X[i] += shiftX[i]
		Y[i] += shiftY[i]
		Z[i] += shiftZ[i]
	}

	x = utils.NewVector(Np, X)
	y = utils.NewVector(Np, Y)
	z = utils.NewVector(Np, Z)

	// Apply shift
	for i := 0; i < Np; i++ {
		X[i] += shiftX[i]
		Y[i] += shiftY[i]
		Z[i] += shiftZ[i]
	}

	// Transform to reference element coordinates
	r, s, t = XYZtoRST(X, Y, Z)

	return r, s, t
}

// EquiNodes3D generates equispaced nodes on the reference tetrahedron
func EquiNodes3D(N int) (r, s, t utils.Vector) {
	Np := (N + 1) * (N + 2) * (N + 3) / 6
	rr := make([]float64, Np)
	ss := make([]float64, Np)
	tt := make([]float64, Np)

	sk := 0
	for n := 0; n <= N; n++ {
		for m := 0; m <= N-n; m++ {
			for l := 0; l <= N-n-m; l++ {
				rr[sk] = -1.0 + 2.0*float64(l)/float64(N)
				ss[sk] = -1.0 + 2.0*float64(m)/float64(N)
				tt[sk] = -1.0 + 2.0*float64(n)/float64(N)
				sk++
			}
		}
	}

	r = utils.NewVector(Np, rr)
	s = utils.NewVector(Np, ss)
	t = utils.NewVector(Np, tt)

	return r, s, t
}

// WarpShiftFace3D computes warp factor used in creating 3D Warp & Blend nodes
func WarpShiftFace3D(p int, pval, pval2 float64, L1, L2, L3, L4 []float64) (warpx, warpy []float64) {
	// Compute warp factors using evalshift
	warpx, warpy = evalshift(p, pval, L2, L3, L4)
	return warpx, warpy
}

// XYZtoRST transforms from equilateral to reference tetrahedron coordinates
func XYZtoRST(X, Y, Z []float64) (r, s, t utils.Vector) {
	n := len(X)
	rr := make([]float64, n)
	ss := make([]float64, n)
	tt := make([]float64, n)

	v1 := []float64{-1, -1 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v2 := []float64{1, -1 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v3 := []float64{0, 2 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v4 := []float64{0, 0, 3 / math.Sqrt(6)}

	for i := 0; i < n; i++ {
		// Solve for barycentric coordinates
		A := utils.NewMatrix(3, 3)
		A.Set(0, 0, v2[0]-v1[0])
		A.Set(0, 1, v3[0]-v1[0])
		A.Set(0, 2, v4[0]-v1[0])
		A.Set(1, 0, v2[1]-v1[1])
		A.Set(1, 1, v3[1]-v1[1])
		A.Set(1, 2, v4[1]-v1[1])
		A.Set(2, 0, v2[2]-v1[2])
		A.Set(2, 1, v3[2]-v1[2])
		A.Set(2, 2, v4[2]-v1[2])

		b := utils.NewMatrix(3, 1)
		b.Set(0, 0, X[i]-v1[0])
		b.Set(1, 0, Y[i]-v1[1])
		b.Set(2, 0, Z[i]-v1[2])

		lambda := A.LUSolve(b)

		L1 := lambda.At(0, 0)
		L2 := lambda.At(1, 0)
		L3 := lambda.At(2, 0)
		L4 := 1 - L1 - L2 - L3

		// Convert to reference coordinates
		rr[i] = -L4 - L2 - L3 + L1
		ss[i] = -L4 - L1 + L2 - L3
		tt[i] = -L4 - L1 - L2 + L3
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

// Simplex3DP evaluates 3D polynomial on simplex at (r,s,t) of order (i,j,k)
// Exactly matching C++ implementation
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

// BuildFmask3D builds face node masks for all 4 faces of tetrahedron
// BuildFmask3D builds face node masks for all 4 faces of tetrahedron
// BuildFmask3D builds face node masks for all 4 faces of tetrahedron
func BuildFmask3D(r, s, t utils.Vector, N int) [][]int {
	Np := r.Len()
	fmask := make([][]int, 4)

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

	return fmask
}

// Lift3D computes lift matrix for surface integrals
func Lift3D(N int, fmask [][]int, V utils.Matrix, R, S, T utils.Vector) utils.Matrix {
	Np := V.Rows()
	Nfaces := 4
	Nfp := (N + 1) * (N + 2) / 2 // Number of face points

	// Create face mass matrix
	Emat := utils.NewMatrix(Np, Nfaces*Nfp)

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
		for i := 0; i < Nfp; i++ {
			for j := 0; j < Nfp; j++ {
				// C++ assigns massFace to Emat(idr, JJ)
				// idr = Fmask(All,face) = fmask[face]
				// JJ columns are face*Nfp+j (0-based)
				Emat.Set(fmask[face][i], face*Nfp+j, massFace.At(i, j))
			}
		}
	}

	// Compute LIFT = V*V'*Emat
	LIFT := V.Mul(V.Transpose()).Mul(Emat)

	return LIFT
}

// GeometricFactors3D computes metric terms for physical elements
func GeometricFactors3D(x, y, z, Dr, Ds, Dt utils.Matrix) *GeometricFactors {
	xr := Dr.Mul(x)
	xs := Ds.Mul(x)
	xt := Dt.Mul(x)
	yr := Dr.Mul(y)
	ys := Ds.Mul(y)
	yt := Dt.Mul(y)
	zr := Dr.Mul(z)
	zs := Ds.Mul(z)
	zt := Dt.Mul(z)

	// Compute Jacobian
	J := utils.NewMatrix(xr.Rows(), xr.Cols())
	for i := 0; i < xr.Rows(); i++ {
		for j := 0; j < xr.Cols(); j++ {
			J.Set(i, j, xr.At(i, j)*(ys.At(i, j)*zt.At(i, j)-zs.At(i, j)*yt.At(i, j))-
				yr.At(i, j)*(xs.At(i, j)*zt.At(i, j)-zs.At(i, j)*xt.At(i, j))+
				zr.At(i, j)*(xs.At(i, j)*yt.At(i, j)-ys.At(i, j)*xt.At(i, j)))
		}
	}

	// Compute inverse metric terms
	rx := utils.NewMatrix(xr.Rows(), xr.Cols())
	ry := utils.NewMatrix(xr.Rows(), xr.Cols())
	rz := utils.NewMatrix(xr.Rows(), xr.Cols())
	sx := utils.NewMatrix(xr.Rows(), xr.Cols())
	sy := utils.NewMatrix(xr.Rows(), xr.Cols())
	sz := utils.NewMatrix(xr.Rows(), xr.Cols())
	tx := utils.NewMatrix(xr.Rows(), xr.Cols())
	ty := utils.NewMatrix(xr.Rows(), xr.Cols())
	tz := utils.NewMatrix(xr.Rows(), xr.Cols())

	for i := 0; i < xr.Rows(); i++ {
		for j := 0; j < xr.Cols(); j++ {
			rx.Set(i, j, (ys.At(i, j)*zt.At(i, j)-zs.At(i, j)*yt.At(i, j))/J.At(i, j))
			ry.Set(i, j, -(xs.At(i, j)*zt.At(i, j)-zs.At(i, j)*xt.At(i, j))/J.At(i, j))
			rz.Set(i, j, (xs.At(i, j)*yt.At(i, j)-ys.At(i, j)*xt.At(i, j))/J.At(i, j))

			sx.Set(i, j, -(yr.At(i, j)*zt.At(i, j)-zr.At(i, j)*yt.At(i, j))/J.At(i, j))
			sy.Set(i, j, (xr.At(i, j)*zt.At(i, j)-zr.At(i, j)*xt.At(i, j))/J.At(i, j))
			sz.Set(i, j, -(xr.At(i, j)*yt.At(i, j)-yr.At(i, j)*xt.At(i, j))/J.At(i, j))

			tx.Set(i, j, (yr.At(i, j)*zs.At(i, j)-zr.At(i, j)*ys.At(i, j))/J.At(i, j))
			ty.Set(i, j, -(xr.At(i, j)*zs.At(i, j)-zr.At(i, j)*xs.At(i, j))/J.At(i, j))
			tz.Set(i, j, (xr.At(i, j)*ys.At(i, j)-yr.At(i, j)*xs.At(i, j))/J.At(i, j))
		}
	}

	return &GeometricFactors{
		X: x, Y: y, Z: z,
		Rx: rx, Ry: ry, Rz: rz,
		Sx: sx, Sy: sy, Sz: sz,
		Tx: tx, Ty: ty, Tz: tz,
		J: J,
	}
}

// Normals3D computes outward facing normals at element faces
func Normals3D(geom *GeometricFactors, fmask [][]int, K int) *FaceGeometricFactors {
	Nfaces := 4
	Nfp := len(fmask[0])

	nx := utils.NewMatrix(Nfp*Nfaces, K)
	ny := utils.NewMatrix(Nfp*Nfaces, K)
	nz := utils.NewMatrix(Nfp*Nfaces, K)
	sJ := utils.NewMatrix(Nfp*Nfaces, K)

	// Face 1: t = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[0][i]
			nx.Set(i, k, -geom.Tx.At(vid, k))
			ny.Set(i, k, -geom.Ty.At(vid, k))
			nz.Set(i, k, -geom.Tz.At(vid, k))
			sJ.Set(i, k, geom.J.At(vid, k)*math.Sqrt(
				nx.At(i, k)*nx.At(i, k)+ny.At(i, k)*ny.At(i, k)+nz.At(i, k)*nz.At(i, k)))

			// Normalize
			nx.Set(i, k, nx.At(i, k)/sJ.At(i, k)*geom.J.At(vid, k))
			ny.Set(i, k, ny.At(i, k)/sJ.At(i, k)*geom.J.At(vid, k))
			nz.Set(i, k, nz.At(i, k)/sJ.At(i, k)*geom.J.At(vid, k))
		}
	}

	// Face 2: s = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[1][i]
			nx.Set(Nfp+i, k, -geom.Sx.At(vid, k))
			ny.Set(Nfp+i, k, -geom.Sy.At(vid, k))
			nz.Set(Nfp+i, k, -geom.Sz.At(vid, k))
			sJ.Set(Nfp+i, k, geom.J.At(vid, k)*math.Sqrt(
				nx.At(Nfp+i, k)*nx.At(Nfp+i, k)+ny.At(Nfp+i, k)*ny.At(Nfp+i, k)+nz.At(Nfp+i, k)*nz.At(Nfp+i, k)))

			// Normalize
			nx.Set(Nfp+i, k, nx.At(Nfp+i, k)/sJ.At(Nfp+i, k)*geom.J.At(vid, k))
			ny.Set(Nfp+i, k, ny.At(Nfp+i, k)/sJ.At(Nfp+i, k)*geom.J.At(vid, k))
			nz.Set(Nfp+i, k, nz.At(Nfp+i, k)/sJ.At(Nfp+i, k)*geom.J.At(vid, k))
		}
	}

	// Face 3: r+s+t = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[2][i]
			nx.Set(2*Nfp+i, k, geom.Rx.At(vid, k)+geom.Sx.At(vid, k)+geom.Tx.At(vid, k))
			ny.Set(2*Nfp+i, k, geom.Ry.At(vid, k)+geom.Sy.At(vid, k)+geom.Ty.At(vid, k))
			nz.Set(2*Nfp+i, k, geom.Rz.At(vid, k)+geom.Sz.At(vid, k)+geom.Tz.At(vid, k))
			sJ.Set(2*Nfp+i, k, geom.J.At(vid, k)*math.Sqrt(
				nx.At(2*Nfp+i, k)*nx.At(2*Nfp+i, k)+ny.At(2*Nfp+i, k)*ny.At(2*Nfp+i, k)+nz.At(2*Nfp+i, k)*nz.At(2*Nfp+i, k)))

			// Normalize
			nx.Set(2*Nfp+i, k, nx.At(2*Nfp+i, k)/sJ.At(2*Nfp+i, k)*geom.J.At(vid, k))
			ny.Set(2*Nfp+i, k, ny.At(2*Nfp+i, k)/sJ.At(2*Nfp+i, k)*geom.J.At(vid, k))
			nz.Set(2*Nfp+i, k, nz.At(2*Nfp+i, k)/sJ.At(2*Nfp+i, k)*geom.J.At(vid, k))
		}
	}

	// Face 4: r = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[3][i]
			nx.Set(3*Nfp+i, k, -geom.Rx.At(vid, k))
			ny.Set(3*Nfp+i, k, -geom.Ry.At(vid, k))
			nz.Set(3*Nfp+i, k, -geom.Rz.At(vid, k))
			sJ.Set(3*Nfp+i, k, geom.J.At(vid, k)*math.Sqrt(
				nx.At(3*Nfp+i, k)*nx.At(3*Nfp+i, k)+ny.At(3*Nfp+i, k)*ny.At(3*Nfp+i, k)+nz.At(3*Nfp+i, k)*nz.At(3*Nfp+i, k)))

			// Normalize
			nx.Set(3*Nfp+i, k, nx.At(3*Nfp+i, k)/sJ.At(3*Nfp+i, k)*geom.J.At(vid, k))
			ny.Set(3*Nfp+i, k, ny.At(3*Nfp+i, k)/sJ.At(3*Nfp+i, k)*geom.J.At(vid, k))
			nz.Set(3*Nfp+i, k, nz.At(3*Nfp+i, k)/sJ.At(3*Nfp+i, k)*geom.J.At(vid, k))
		}
	}

	// Compute Fscale
	Fscale := utils.NewMatrix(Nfp*Nfaces, K)
	for f := 0; f < Nfaces; f++ {
		for i := 0; i < Nfp; i++ {
			for k := 0; k < K; k++ {
				vid := fmask[f%4][i]
				Fscale.Set(f*Nfp+i, k, sJ.At(f*Nfp+i, k)/geom.J.At(vid, k))
			}
		}
	}

	return &FaceGeometricFactors{
		Nx: nx, Ny: ny, Nz: nz,
		SJ:     sJ,
		Fscale: Fscale,
	}
}

// BuildMaps3D creates connectivity maps for DG
func BuildMaps3D(K int, Np int, Nfp int, fmask [][]int, x, y, z utils.Matrix,
	conn *ConnectivityArrays) (VmapM, VmapP, MapM, MapP []int) {

	Nfaces := 4

	// Build global face to node mapping
	nodeIds := make(map[[3]float64]int)
	NODETOL := 1e-7

	// Unique node identification
	Nnodes := 0
	for k := 0; k < K; k++ {
		for n := 0; n < Np; n++ {
			key := [3]float64{
				math.Round(x.At(n, k)/NODETOL) * NODETOL,
				math.Round(y.At(n, k)/NODETOL) * NODETOL,
				math.Round(z.At(n, k)/NODETOL) * NODETOL,
			}
			if _, exists := nodeIds[key]; !exists {
				nodeIds[key] = Nnodes
				Nnodes++
			}
		}
	}

	// Volume to global node mapping
	vmapM := make([]int, Nfp*Nfaces*K)
	vmapP := make([]int, Nfp*Nfaces*K)
	mapM := make([]int, Nfp*Nfaces*K)
	mapP := make([]int, Nfp*Nfaces*K)

	// Initialize with identity mapping
	for i := range vmapM {
		vmapM[i] = i
		vmapP[i] = i
		mapM[i] = i
		mapP[i] = i
	}

	// Create face to face mapping
	for k1 := 0; k1 < K; k1++ {
		for f1 := 0; f1 < Nfaces; f1++ {
			k2 := conn.EToE[k1][f1]
			f2 := conn.EToF[k1][f1]

			if k2 < k1 || (k2 == k1 && f2 < f1) {
				continue // Only process each face pair once
			}

			// Find matching nodes
			for i := 0; i < Nfp; i++ {
				n1 := fmask[f1][i]
				key1 := [3]float64{
					math.Round(x.At(n1, k1)/NODETOL) * NODETOL,
					math.Round(y.At(n1, k1)/NODETOL) * NODETOL,
					math.Round(z.At(n1, k1)/NODETOL) * NODETOL,
				}

				for j := 0; j < Nfp; j++ {
					n2 := fmask[f2][j]
					key2 := [3]float64{
						math.Round(x.At(n2, k2)/NODETOL) * NODETOL,
						math.Round(y.At(n2, k2)/NODETOL) * NODETOL,
						math.Round(z.At(n2, k2)/NODETOL) * NODETOL,
					}

					if key1 == key2 {
						id1 := k1*Nfaces*Nfp + f1*Nfp + i
						id2 := k2*Nfaces*Nfp + f2*Nfp + j

						vmapM[id1] = k1*Np + n1
						vmapM[id2] = k2*Np + n2
						vmapP[id1] = k2*Np + n2
						vmapP[id2] = k1*Np + n1

						mapP[id1] = id2
						mapP[id2] = id1
						break
					}
				}
			}
		}
	}

	return vmapM, vmapP, mapM, mapP
}

// Connect3D builds element to element connectivity
func Connect3D(EToV [][]int) *ConnectivityArrays {
	K := len(EToV)
	Nfaces := 4

	// Build face to vertex mapping for tetrahedron
	vn := [][]int{
		{0, 1, 2}, // Face 1
		{0, 1, 3}, // Face 2
		{1, 2, 3}, // Face 3
		{0, 2, 3}, // Face 4
	}

	// Create face to element+face mapping
	faces := make(map[[3]int]struct{ elem, face int })

	for k := 0; k < K; k++ {
		for f := 0; f < Nfaces; f++ {
			// Get vertices of this face
			v := make([]int, 3)
			for i := 0; i < 3; i++ {
				v[i] = EToV[k][vn[f][i]]
			}

			// Sort vertices to create unique key
			sort.Ints(v)
			key := [3]int{v[0], v[1], v[2]}

			if match, exists := faces[key]; exists {
				// Found matching face
				// Update both elements
				if k != match.elem {
					// Different elements share this face
					// This is an internal face
				}
			} else {
				faces[key] = struct{ elem, face int }{k, f}
			}
		}
	}

	// Build connectivity arrays
	EToE := make([][]int, K)
	EToF := make([][]int, K)

	for k := 0; k < K; k++ {
		EToE[k] = make([]int, Nfaces)
		EToF[k] = make([]int, Nfaces)

		// Initialize with self-connectivity
		for f := 0; f < Nfaces; f++ {
			EToE[k][f] = k
			EToF[k][f] = f
		}

		// Find actual connections
		for f := 0; f < Nfaces; f++ {
			v := make([]int, 3)
			for i := 0; i < 3; i++ {
				v[i] = EToV[k][vn[f][i]]
			}
			sort.Ints(v)
			key := [3]int{v[0], v[1], v[2]}

			// Search all faces for match
			for fkey, fdata := range faces {
				if fkey == key && fdata.elem != k {
					EToE[k][f] = fdata.elem
					EToF[k][f] = fdata.face
					break
				}
			}
		}
	}

	return &ConnectivityArrays{
		EToE: EToE,
		EToF: EToF,
	}
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

// evalshift computes two-dimensional Warp & Blend transform
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

// evalwarp evaluates the warp function
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
