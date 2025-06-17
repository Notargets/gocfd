package gonudg

import "math"

// XYZtoRST transfers from (x,y,z) in equilateral tetrahedron
// to (r,s,t) coordinates in standard tetrahedron
// This is the 0-based index version of the C++ xyztorst function
func XYZtoRST(X, Y, Z []float64) (r, s, t []float64) {
	sqrt3 := math.Sqrt(3.0)
	sqrt6 := math.Sqrt(6.0)
	Nc := len(X)
	
	// Vertices of the equilateral tetrahedron
	v1 := []float64{-1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v2 := []float64{1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v3 := []float64{0.0, 2.0 / sqrt3, -1.0 / sqrt6}
	v4 := []float64{0.0, 0.0, 3.0 / sqrt6}
	
	// Build transformation matrix A
	// A = 0.5 * [v2-v1, v3-v1, v4-v1]
	A := make([][]float64, 3)
	for i := 0; i < 3; i++ {
		A[i] = make([]float64, 3)
	}
	
	A[0][0] = 0.5 * (v2[0] - v1[0])
	A[0][1] = 0.5 * (v3[0] - v1[0])
	A[0][2] = 0.5 * (v4[0] - v1[0])
	A[1][0] = 0.5 * (v2[1] - v1[1])
	A[1][1] = 0.5 * (v3[1] - v1[1])
	A[1][2] = 0.5 * (v4[1] - v1[1])
	A[2][0] = 0.5 * (v2[2] - v1[2])
	A[2][1] = 0.5 * (v3[2] - v1[2])
	A[2][2] = 0.5 * (v4[2] - v1[2])
	
	// Compute offset vector
	// offset = 0.5*(v2+v3+v4-v1)
	offset := []float64{
		0.5 * (v2[0] + v3[0] + v4[0] - v1[0]),
		0.5 * (v2[1] + v3[1] + v4[1] - v1[1]),
		0.5 * (v2[2] + v3[2] + v4[2] - v1[2]),
	}
	
	// Invert matrix A
	Ainv := inverse3x3(A)
	
	// Initialize output
	r = make([]float64, Nc)
	s = make([]float64, Nc)
	t = make([]float64, Nc)
	
	// Solve for each point
	for n := 0; n < Nc; n++ {
		// rhs = [X,Y,Z] - offset
		rhs := []float64{
			X[n] - offset[0],
			Y[n] - offset[1],
			Z[n] - offset[2],
		}
		
		// [r,s,t] = A^{-1} * rhs
		r[n] = Ainv[0][0]*rhs[0] + Ainv[0][1]*rhs[1] + Ainv[0][2]*rhs[2]
		s[n] = Ainv[1][0]*rhs[0] + Ainv[1][1]*rhs[1] + Ainv[1][2]*rhs[2]
		t[n] = Ainv[2][0]*rhs[0] + Ainv[2][1]*rhs[1] + Ainv[2][2]*rhs[2]
	}
	
	return r, s, t
}

// RSTtoXYZ is the inverse transformation from (r,s,t) to (x,y,z)
// in the equilateral tetrahedron
func RSTtoXYZ(r, s, t []float64) (X, Y, Z []float64) {
	sqrt3 := math.Sqrt(3.0)
	sqrt6 := math.Sqrt(6.0)
	Nc := len(r)
	
	// Vertices of the equilateral tetrahedron
	v1 := []float64{-1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v2 := []float64{1.0, -1.0 / sqrt3, -1.0 / sqrt6}
	v3 := []float64{0.0, 2.0 / sqrt3, -1.0 / sqrt6}
	v4 := []float64{0.0, 0.0, 3.0 / sqrt6}
	
	// Initialize output
	X = make([]float64, Nc)
	Y = make([]float64, Nc)
	Z = make([]float64, Nc)
	
	// Transform each point using barycentric coordinates
	// The transformation is:
	// [x,y,z] = 0.5*(-(1+r+s+t)*v1 + (1+r)*v2 + (1+s)*v3 + (1+t)*v4)
	for n := 0; n < Nc; n++ {
		coeff1 := 0.5 * (-(1.0 + r[n] + s[n] + t[n]))
		coeff2 := 0.5 * (1.0 + r[n])
		coeff3 := 0.5 * (1.0 + s[n])
		coeff4 := 0.5 * (1.0 + t[n])
		
		X[n] = coeff1*v1[0] + coeff2*v2[0] + coeff3*v3[0] + coeff4*v4[0]
		Y[n] = coeff1*v1[1] + coeff2*v2[1] + coeff3*v3[1] + coeff4*v4[1]
		Z[n] = coeff1*v1[2] + coeff2*v2[2] + coeff3*v3[2] + coeff4*v4[2]
	}
	
	return X, Y, Z
}

// Helper function to invert a 3x3 matrix
func inverse3x3(A [][]float64) [][]float64 {
	// Compute determinant
	det := A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) -
		A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]) +
		A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
	
	// Compute inverse
	inv := make([][]float64, 3)
	for i := 0; i < 3; i++ {
		inv[i] = make([]float64, 3)
	}
	
	inv[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1]) / det
	inv[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2]) / det
	inv[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1]) / det
	inv[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2]) / det
	inv[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0]) / det
	inv[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2]) / det
	inv[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0]) / det
	inv[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1]) / det
	inv[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0]) / det
	
	return inv
}

