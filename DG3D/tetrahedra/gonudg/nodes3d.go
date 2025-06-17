package gonudg

// INDEXING NOTE: Original C++ code uses 1-based indexing to emulate Matlab behavior.
// This Go port uses standard 0-based indexing. Example conversions:
//   C++: sk = 1; V3D(All,sk) = ...    ->    Go: sk = 0; V3D.SetCol(sk, ...)
//   C++: Fmask[1] (first face)        ->    Go: Fmask[0] (first face)
// The indexing has been correctly translated throughout this port.

import (
	"math"
)

// Nodes3D computes Warp & Blend nodes
// Input: p = polynomial order of interpolant
// Output: X,Y,Z vectors of node coordinates in equilateral tetrahedron
// This is the 0-based index version of the C++ Nodes3D function
func Nodes3D(p int) (X, Y, Z []float64) {
	// Choose optimized blending parameter
	alphastore := []float64{
		0.0000, 0.0000, 0.0000, 0.1002, 1.1332,
		1.5608, 1.3413, 1.2577, 1.1603, 1.10153,
		0.6080, 0.4523, 0.8856, 0.8717, 0.9655,
	}

	alpha := 1.0
	if p < 15 && p >= 0 {
		alpha = alphastore[p]
	}

	// Total number of nodes and tolerance
	Np := (p + 1) * (p + 2) * (p + 3) / 6
	tol := 1e-10
	sqrt3 := math.Sqrt(3.0)
	sqrt6 := math.Sqrt(6.0)

	// Create equidistributed nodes
	r, s, t := EquiNodes3D(p)

	// Compute barycentric coordinates
	L1 := make([]float64, Np)
	L2 := make([]float64, Np)
	L3 := make([]float64, Np)
	L4 := make([]float64, Np)

	for i := 0; i < Np; i++ {
		L1[i] = (1.0 + t[i]) / 2.0
		L2[i] = (1.0 + s[i]) / 2.0
		L3[i] = -(1.0 + r[i] + s[i] + t[i]) / 2.0
		L4[i] = (1.0 + r[i]) / 2.0
	}

	// Set vertices of tetrahedron
	v1 := []float64{-1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v2 := []float64{1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v3 := []float64{0.0, 2.0 / sqrt3, -1.0 / sqrt6}
	v4 := []float64{0.0, 0.0, 3.0 / sqrt6}

	// Orthogonal axis tangents on faces 1-4
	t1 := make([][]float64, 4)
	t2 := make([][]float64, 4)

	// Face 1
	t1[0] = vecSub(v2, v1)
	t2[0] = vecSub(v3, vecScale(0.5, vecAdd(v1, v2)))

	// Face 2
	t1[1] = vecSub(v2, v1)
	t2[1] = vecSub(v4, vecScale(0.5, vecAdd(v1, v2)))

	// Face 3
	t1[2] = vecSub(v3, v2)
	t2[2] = vecSub(v4, vecScale(0.5, vecAdd(v2, v3)))

	// Face 4
	t1[3] = vecSub(v3, v1)
	t2[3] = vecSub(v4, vecScale(0.5, vecAdd(v1, v3)))

	// Normalize tangents
	for n := 0; n < 4; n++ {
		t1[n] = vecNormalize(t1[n])
		t2[n] = vecNormalize(t2[n])
	}

	// Form undeformed coordinates
	X = make([]float64, Np)
	Y = make([]float64, Np)
	Z = make([]float64, Np)

	for i := 0; i < Np; i++ {
		X[i] = L3[i]*v1[0] + L4[i]*v2[0] + L2[i]*v3[0] + L1[i]*v4[0]
		Y[i] = L3[i]*v1[1] + L4[i]*v2[1] + L2[i]*v3[1] + L1[i]*v4[1]
		Z[i] = L3[i]*v1[2] + L4[i]*v2[2] + L2[i]*v3[2] + L1[i]*v4[2]
	}

	// Initialize shift
	shiftX := make([]float64, Np)
	shiftY := make([]float64, Np)
	shiftZ := make([]float64, Np)

	// Warp and blend for each face
	for face := 0; face < 4; face++ {
		var La, Lb, Lc, Ld []float64

		// C++ uses 1-based indexing, adjust to 0-based
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
		warp1, warp2 := WarpShiftFace3D(p, alpha, alpha, La, Lb, Lc, Ld)

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

				// If node is on face but not on 3 edges
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

	return X, Y, Z
}

// EquiNodes3D creates equidistributed nodes on the reference tetrahedron
func EquiNodes3D(p int) (r, s, t []float64) {
	Np := (p + 1) * (p + 2) * (p + 3) / 6
	r = make([]float64, Np)
	s = make([]float64, Np)
	t = make([]float64, Np)

	// Special case for p=0: single node at centroid
	if p == 0 {
		r[0] = -0.5
		s[0] = -0.5
		t[0] = -0.5
		return r, s, t
	}

	sk := 0
	for n := 0; n <= p; n++ {
		for m := 0; m <= p-n; m++ {
			for l := 0; l <= p-n-m; l++ {
				r[sk] = -1.0 + 2.0*float64(l)/float64(p)
				s[sk] = -1.0 + 2.0*float64(m)/float64(p)
				t[sk] = -1.0 + 2.0*float64(n)/float64(p)
				sk++
			}
		}
	}

	return r, s, t
}

// EquidistributedNodes3D creates evenly spaced nodes in the reference tetrahedron
// Alias for EquiNodes3D for compatibility
func EquidistributedNodes3D(N int) (r, s, t []float64) {
	return EquiNodes3D(N)
}

// WarpShiftFace3D computes warp shift for a face
func WarpShiftFace3D(p int, pval, pval2 float64, La, Lb, Lc, Ld []float64) (warpx, warpy []float64) {
	// Use evalshift with the appropriate parameters
	// The function uses Lb, Lc, Ld for the warping
	warpx, warpy = evalshift(p, pval, Lb, Lc, Ld)
	return warpx, warpy
}

// evalshift evaluates the warp shift for a face
func evalshift(p int, pval float64, L1, L2, L3 []float64) (dx, dy []float64) {
	n := len(L1)
	dx = make([]float64, n)
	dy = make([]float64, n)

	// Compute Gauss-Lobatto-Legendre node distribution
	gaussX := JacobiGL(0, 0, p)
	// Negate the values (matching C++: gaussX = -JacobiGL(0,0,p))
	for i := range gaussX {
		gaussX[i] = -gaussX[i]
	}

	// Compute blending function at each node for each edge
	blend1 := make([]float64, n)
	blend2 := make([]float64, n)
	blend3 := make([]float64, n)

	for i := 0; i < n; i++ {
		blend1[i] = L2[i] * L3[i]
		blend2[i] = L1[i] * L3[i]
		blend3[i] = L1[i] * L2[i]
	}

	// Amount of warp for each node, for each edge
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

	// Combine blend & warp
	warp1 := make([]float64, n)
	warp2 := make([]float64, n)
	warp3 := make([]float64, n)

	for i := 0; i < n; i++ {
		warp1[i] = blend1[i] * warpfactor1[i] * (1.0 + pval*pval*L1[i]*L1[i])
		warp2[i] = blend2[i] * warpfactor2[i] * (1.0 + pval*pval*L2[i]*L2[i])
		warp3[i] = blend3[i] * warpfactor3[i] * (1.0 + pval*pval*L3[i]*L3[i])
	}

	// Evaluate shift in equilateral triangle
	const TWOPI = 2.0 * math.Pi
	const FOURPI = 4.0 * math.Pi

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

// Vector helper functions

func vecSub(a, b []float64) []float64 {
	result := make([]float64, len(a))
	for i := range a {
		result[i] = a[i] - b[i]
	}
	return result
}

func vecAdd(a, b []float64) []float64 {
	result := make([]float64, len(a))
	for i := range a {
		result[i] = a[i] + b[i]
	}
	return result
}

func vecScale(s float64, v []float64) []float64 {
	result := make([]float64, len(v))
	for i := range v {
		result[i] = s * v[i]
	}
	return result
}

func vecNorm(v []float64) float64 {
	sum := 0.0
	for _, val := range v {
		sum += val * val
	}
	return math.Sqrt(sum)
}

func vecNormalize(v []float64) []float64 {
	norm := vecNorm(v)
	result := make([]float64, len(v))
	for i := range v {
		result[i] = v[i] / norm
	}
	return result
}
