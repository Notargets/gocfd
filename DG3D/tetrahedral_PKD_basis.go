package DG3D

import (
	"math"
	"sort"

	"github.com/notargets/gocfd/DG1D"
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
			Veq.Set(i, j, DG1D.JacobiP(utils.NewVector(1, []float64{req[i]}),
				0, 0, j)[0])
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

	// Create equidistant nodes
	L1, L2, L3, L4 := make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np)

	sk := 0
	for n := 0; n <= N; n++ {
		for m := 0; m <= N-n; m++ {
			for l := 0; l <= N-n-m; l++ {
				L1[sk] = float64(l) / float64(N)
				L2[sk] = float64(m) / float64(N)
				L3[sk] = float64(n) / float64(N)
				L4[sk] = 1.0 - L1[sk] - L2[sk] - L3[sk]
				sk++
			}
		}
	}

	// Set vertices of tetrahedron
	v1 := []float64{-1, -1 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v2 := []float64{1, -1 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v3 := []float64{0, 2 / math.Sqrt(3), -1 / math.Sqrt(6)}
	v4 := []float64{0, 0, 3 / math.Sqrt(6)}

	// Create initial equilateral coordinates
	X := make([]float64, Np)
	Y := make([]float64, Np)
	Z := make([]float64, Np)

	for i := 0; i < Np; i++ {
		X[i] = L1[i]*v2[0] + L2[i]*v3[0] + L3[i]*v4[0] + L4[i]*v1[0]
		Y[i] = L1[i]*v2[1] + L2[i]*v3[1] + L3[i]*v4[1] + L4[i]*v1[1]
		Z[i] = L1[i]*v2[2] + L2[i]*v3[2] + L3[i]*v4[2] + L4[i]*v1[2]
	}

	// Face 1: L1 = 0
	faceMask := make([]bool, Np)
	for i := 0; i < Np; i++ {
		faceMask[i] = math.Abs(L1[i]) < 1e-10
	}

	// Create orthogonal tangent vectors for face 1
	t1 := []float64{v3[0] - v4[0], v3[1] - v4[1], v3[2] - v4[2]}
	t2 := []float64{v2[0] - v4[0], v2[1] - v4[1], v2[2] - v4[2]}

	WarpShiftFace3D(N, faceMask, L2, L3, L4, t1, t2, X, Y, Z)

	// Face 2: L2 = 0
	for i := 0; i < Np; i++ {
		faceMask[i] = math.Abs(L2[i]) < 1e-10
	}

	t1 = []float64{v1[0] - v3[0], v1[1] - v3[1], v1[2] - v3[2]}
	t2 = []float64{v4[0] - v3[0], v4[1] - v3[1], v4[2] - v3[2]}

	WarpShiftFace3D(N, faceMask, L1, L3, L4, t1, t2, X, Y, Z)

	// Face 3: L3 = 0
	for i := 0; i < Np; i++ {
		faceMask[i] = math.Abs(L3[i]) < 1e-10
	}

	t1 = []float64{v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]}
	t2 = []float64{v4[0] - v1[0], v4[1] - v1[1], v4[2] - v1[2]}

	WarpShiftFace3D(N, faceMask, L1, L2, L4, t1, t2, X, Y, Z)

	// Face 4: L4 = 0
	for i := 0; i < Np; i++ {
		faceMask[i] = math.Abs(L4[i]) < 1e-10
	}

	t1 = []float64{v3[0] - v2[0], v3[1] - v2[1], v3[2] - v2[2]}
	t2 = []float64{v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]}

	WarpShiftFace3D(N, faceMask, L1, L2, L3, t1, t2, X, Y, Z)

	// Transform to reference element coordinates
	r, s, t := XYZtoRST(X, Y, Z)

	return r, s, t
}

// WarpShiftFace3D applies warp and blend to a face
func WarpShiftFace3D(N int, faceMask []bool, L1, L2, L3, t1, t2 []float64, X, Y, Z []float64) {
	Np := len(X)

	// Compute blend factor for each node
	blend := make([]float64, Np)
	for i := 0; i < Np; i++ {
		blend[i] = 4 * L1[i] * L2[i] * L3[i]
	}

	// For each face node
	for i := 0; i < Np; i++ {
		if !faceMask[i] {
			continue
		}

		// Compute r,s coordinates on face
		denom := L1[i] + L2[i] + L3[i]
		if math.Abs(denom) < 1e-10 {
			continue
		}

		r := (L2[i] - L1[i]) / denom
		s := (L3[i] - L1[i]) / denom

		// Evaluate warp
		warpR := Warpfactor(N, []float64{r})
		warpS := Warpfactor(N, []float64{s})

		// Apply warp and blend
		X[i] += blend[i] * warpR[0] * t1[0]
		Y[i] += blend[i] * warpR[0] * t1[1]
		Z[i] += blend[i] * warpR[0] * t1[2]

		X[i] += blend[i] * warpS[0] * t2[0]
		Y[i] += blend[i] * warpS[0] * t2[1]
		Z[i] += blend[i] * warpS[0] * t2[2]
	}
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

// Simplex3DP evaluates 3D polynomial on simplex at (r,s,t) of order (i,j,k)
// The PKD basis is orthogonal with respect to the weighted inner product
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

	// Face 1: t = -1
	for i := 0; i < Np; i++ {
		if math.Abs(t.At(i)+1) < NODETOL {
			fmask[0] = append(fmask[0], i)
		}
	}

	// Face 2: s = -1
	for i := 0; i < Np; i++ {
		if math.Abs(s.At(i)+1) < NODETOL {
			fmask[1] = append(fmask[1], i)
		}
	}

	// Face 3: r+s+t = -1
	for i := 0; i < Np; i++ {
		if math.Abs(r.At(i)+s.At(i)+t.At(i)+1) < NODETOL {
			fmask[2] = append(fmask[2], i)
		}
	}

	// Face 4: r = -1
	for i := 0; i < Np; i++ {
		if math.Abs(r.At(i)+1) < NODETOL {
			fmask[3] = append(fmask[3], i)
		}
	}

	return fmask
}

// Lift3D computes lift matrix for surface integrals
func Lift3D(N int, fmask [][]int, V utils.Matrix) utils.Matrix {
	Np := V.Rows()
	Nfaces := 4
	Nfp := (N + 1) * (N + 2) / 2 // Number of face points

	// Create face mass matrix
	Emat := utils.NewMatrix(Np, Nfaces*Nfp)

	// Reference triangle nodes for faces
	faceR, faceS := Nodes2D(N)
	faceV := Vandermonde2D(N, faceR, faceS)
	massFace := faceV.Mul(faceV.Transpose()).InverseWithCheck()

	// Build the Emat matrix
	for face := 0; face < Nfaces; face++ {
		if len(fmask[face]) != Nfp {
			panic("face mask size mismatch")
		}

		// Extract face nodes
		for i := 0; i < Nfp; i++ {
			for j := 0; j < Nfp; j++ {
				Emat.Set(fmask[face][i], face*Nfp+j, massFace.At(i, j))
			}
		}
	}

	// Compute LIFT = V*V'*Emat
	LIFT := V.Mul(V.Transpose()).Mul(Emat)

	return LIFT
}

// Nodes2D computes nodes on reference triangle
func Nodes2D(N int) (r, s utils.Vector) {
	Np := (N + 1) * (N + 2) / 2
	rr := make([]float64, Np)
	ss := make([]float64, Np)

	// Create equidistant nodes
	sk := 0
	for n := 0; n <= N; n++ {
		for m := 0; m <= N-n; m++ {
			L1 := float64(m) / float64(N)
			L2 := float64(n) / float64(N)
			L3 := 1.0 - L1 - L2
			_ = L3

			// Convert to xy coordinates
			rr[sk] = 2*L1 - 1
			ss[sk] = 2*L2 - 1
			sk++
		}
	}

	r = utils.NewVector(Np, rr)
	s = utils.NewVector(Np, ss)

	// Apply warping for optimal nodes
	warpNodes2D(&r, &s, N)

	return r, s
}

// warpNodes2D applies warping to 2D triangle nodes
func warpNodes2D(r, s *utils.Vector, N int) {
	// This is a simplified version - implement full warping if needed
	// For now, just return the equidistant nodes
}

// Vandermonde2D computes 2D Vandermonde matrix on triangle
func Vandermonde2D(N int, r, s utils.Vector) utils.Matrix {
	Np := (N + 1) * (N + 2) / 2
	V := utils.NewMatrix(r.Len(), Np)

	// Convert to collapsed coordinates
	a, b := RStoAB(r, s)

	sk := 0
	for i := 0; i <= N; i++ {
		for j := 0; j <= N-i; j++ {
			V.SetCol(sk, Simplex2DP(a, b, i, j))
			sk++
		}
	}

	return V
}

// Simplex2DP evaluates 2D polynomial on simplex
func Simplex2DP(a, b utils.Vector, i, j int) []float64 {
	n := a.Len()
	P := make([]float64, n)

	h1 := DG1D.JacobiP(a, 0, 0, i)
	h2 := DG1D.JacobiP(b, float64(2*i+1), 0, j)

	for idx := 0; idx < n; idx++ {
		P[idx] = math.Sqrt(2.0) * h1[idx] * h2[idx] * math.Pow((1-b.At(idx))/2, float64(i))
	}

	return P
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
// Note: The PKD basis is orthogonal but NOT orthonormal, so this will not be identity
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
