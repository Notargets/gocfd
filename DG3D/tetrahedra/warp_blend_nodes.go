package tetrahedra

import (
	"github.com/notargets/gocfd/utils"
	"math"
)

// Nodes3D computes Warp & Blend nodes in equilateral tetrahedron coordinates
// This matches the C++ implementation exactly
func Nodes3D(N int) (X, Y, Z utils.Vector) {
	// Choose optimized blending parameter
	alphastore := []float64{
		0.0000, 0.0000, 0.0000, 0.1002, 1.1332,
		1.5608, 1.3413, 1.2577, 1.1603, 1.10153,
		0.6080, 0.4523, 0.8856, 0.8717, 0.9655,
	}

	alpha := 1.0
	if N <= 14 {
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
	Xdata := make([]float64, Np)
	Ydata := make([]float64, Np)
	Zdata := make([]float64, Np)

	for i := 0; i < Np; i++ {
		Xdata[i] = L3[i]*v1[0] + L4[i]*v2[0] + L2[i]*v3[0] + L1[i]*v4[0]
		Ydata[i] = L3[i]*v1[1] + L4[i]*v2[1] + L2[i]*v3[1] + L1[i]*v4[1]
		Zdata[i] = L3[i]*v1[2] + L4[i]*v2[2] + L2[i]*v3[2] + L1[i]*v4[2]
	}

	// Shift coordinates for warping
	shiftX := make([]float64, Np)
	shiftY := make([]float64, Np)
	shiftZ := make([]float64, Np)

	// Orthogonal axis tangents on faces 1-4
	t1 := make([][]float64, 4)
	t2 := make([][]float64, 4)

	// C++ uses 1-based face indexing, we use 0-based
	// Face 0 (C++ face 1)
	t1[0] = []float64{v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]}
	t2[0] = []float64{
		v3[0] - 0.5*(v1[0]+v2[0]),
		v3[1] - 0.5*(v1[1]+v2[1]),
		v3[2] - 0.5*(v1[2]+v2[2]),
	}

	// Face 1 (C++ face 2)
	t1[1] = []float64{v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]}
	t2[1] = []float64{
		v4[0] - 0.5*(v1[0]+v2[0]),
		v4[1] - 0.5*(v1[1]+v2[1]),
		v4[2] - 0.5*(v1[2]+v2[2]),
	}

	// Face 2 (C++ face 3)
	t1[2] = []float64{v3[0] - v2[0], v3[1] - v2[1], v3[2] - v2[2]}
	t2[2] = []float64{
		v4[0] - 0.5*(v2[0]+v3[0]),
		v4[1] - 0.5*(v2[1]+v3[1]),
		v4[2] - 0.5*(v2[2]+v3[2]),
	}

	// Face 3 (C++ face 4)
	t1[3] = []float64{v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]}
	t2[3] = []float64{
		v4[0] - 0.5*(v1[0]+v3[0]),
		v4[1] - 0.5*(v1[1]+v3[1]),
		v4[2] - 0.5*(v1[2]+v3[2]),
	}

	// Normalize tangents
	for n := 0; n < 4; n++ {
		norm1 := math.Sqrt(t1[n][0]*t1[n][0] + t1[n][1]*t1[n][1] + t1[n][2]*t1[n][2])
		norm2 := math.Sqrt(t2[n][0]*t2[n][0] + t2[n][1]*t2[n][1] + t2[n][2]*t2[n][2])

		for i := 0; i < 3; i++ {
			t1[n][i] /= norm1
			t2[n][i] /= norm2
		}
	}

	// Warp and blend for each face
	for face := 0; face < 4; face++ {
		var La, Lb, Lc, Ld []float64

		switch face {
		case 0: // C++ face 1
			La, Lb, Lc, Ld = L1, L2, L3, L4
		case 1: // C++ face 2
			La, Lb, Lc, Ld = L2, L1, L3, L4
		case 2: // C++ face 3
			La, Lb, Lc, Ld = L3, L1, L4, L2
		case 3: // C++ face 4
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
		}

		// Modify linear blend
		for i := 0; i < Np; i++ {
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
		Xdata[i] += shiftX[i]
		Ydata[i] += shiftY[i]
		Zdata[i] += shiftZ[i]
	}

	X = utils.NewVector(Np, Xdata)
	Y = utils.NewVector(Np, Ydata)
	Z = utils.NewVector(Np, Zdata)

	return X, Y, Z
}

// xyztorst transfers from (x,y,z) in equilateral tetrahedron
// to (r,s,t) coordinates in standard tetrahedron
// This is a direct translation of the C++ xyztorst function
func xyztorst(X, Y, Z utils.Vector) (r, s, t utils.Vector) {
	Nc := X.Len()
	sqrt3 := math.Sqrt(3.0)
	sqrt6 := math.Sqrt(6.0)

	// Vertices of equilateral tetrahedron
	v1 := []float64{-1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v2 := []float64{1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v3 := []float64{0.0, 2.0 / sqrt3, -1.0 / sqrt6}
	v4 := []float64{0.0, 0.0, 3.0 / sqrt6}

	// Build the transformation matrix A
	A := utils.NewMatrix(3, 3)
	A.Set(0, 0, 0.5*(v2[0]-v1[0]))
	A.Set(0, 1, 0.5*(v3[0]-v1[0]))
	A.Set(0, 2, 0.5*(v4[0]-v1[0]))
	A.Set(1, 0, 0.5*(v2[1]-v1[1]))
	A.Set(1, 1, 0.5*(v3[1]-v1[1]))
	A.Set(1, 2, 0.5*(v4[1]-v1[1]))
	A.Set(2, 0, 0.5*(v2[2]-v1[2]))
	A.Set(2, 1, 0.5*(v3[2]-v1[2]))
	A.Set(2, 2, 0.5*(v4[2]-v1[2]))

	// Prepare right-hand side
	// rhs = tmat1 - 0.5*outer(v2+v3+v4-v1, ones(Nc))
	offset := []float64{
		0.5 * (v2[0] + v3[0] + v4[0] - v1[0]),
		0.5 * (v2[1] + v3[1] + v4[1] - v1[1]),
		0.5 * (v2[2] + v3[2] + v4[2] - v1[2]),
	}

	rData := make([]float64, Nc)
	sData := make([]float64, Nc)
	tData := make([]float64, Nc)

	// Solve for each point
	for i := 0; i < Nc; i++ {
		// Build RHS for this point
		rhs := utils.NewMatrix(3, 1)
		rhs.Set(0, 0, X.At(i)-offset[0])
		rhs.Set(1, 0, Y.At(i)-offset[1])
		rhs.Set(2, 0, Z.At(i)-offset[2])

		// Solve A * RST = rhs
		RST := A.LUSolve(rhs)

		rData[i] = RST.At(0, 0)
		sData[i] = RST.At(1, 0)
		tData[i] = RST.At(2, 0)
	}

	r = utils.NewVector(Nc, rData)
	s = utils.NewVector(Nc, sData)
	t = utils.NewVector(Nc, tData)

	return r, s, t
}

// Update NewTetBasis to match C++ flow
func NewTetBasis(N int) (tb *TetBasis) {
	tb = &TetBasis{
		N:  N,
		Np: (N + 1) * (N + 2) * (N + 3) / 6,
	}

	// Get nodes in equilateral tetrahedron coordinates
	X, Y, Z := Nodes3D(N)

	// Convert to reference tetrahedron coordinates
	tb.R, tb.S, tb.T = xyztorst(X, Y, Z)

	// Build basis matrices
	tb.V = tb.ComputeVandermonde(tb.R, tb.S, tb.T)
	tb.VInv = tb.V.InverseWithCheck()
	tb.LIFT = tb.Lift3D()

	// Build differentiation matrices
	Vr, Vs, Vt := tb.ComputeGradVandermonde(tb.R, tb.S, tb.T)
	tb.Dr = Vr.Mul(tb.VInv)
	tb.Ds = Vs.Mul(tb.VInv)
	tb.Dt = Vt.Mul(tb.VInv)

	// Build mass matrix
	tb.M = tb.VInv.Transpose().Mul(tb.VInv)
	inverse, err := tb.M.Inverse()
	if err != nil {
		return nil
	}
	tb.MInv = inverse

	return
}
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

func WarpShiftFace3D(p int, pval, pval2 float64, L1, L2, L3, L4 []float64) (warpx, warpy []float64) {
	// Compute warp factors using evalshift
	warpx, warpy = evalshift(p, pval, L2, L3, L4)
	return warpx, warpy
}

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
