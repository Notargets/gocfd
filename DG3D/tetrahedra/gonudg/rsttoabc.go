package gonudg

import "math"

// RSTtoABC transfers from (r,s,t) to (a,b,c) coordinates in tetrahedron
// This matches the C++ implementation exactly - no special cases
func RSTtoABC(r, s, t []float64) (a, b, c []float64) {
	Np := len(r)
	a = make([]float64, Np)
	b = make([]float64, Np)
	c = make([]float64, Np)

	const tol = 1e-12

	for n := 0; n < Np; n++ {
		// Handle a coordinate
		if math.Abs(s[n]+t[n]) < tol {
			a[n] = -1.0
		} else {
			a[n] = 2.0*(1.0+r[n])/(-s[n]-t[n]) - 1.0
		}

		// Handle b coordinate
		if math.Abs(t[n]-1.0) < tol {
			b[n] = -1.0
		} else {
			b[n] = 2.0*(1.0+s[n])/(1.0-t[n]) - 1.0
		}

		c[n] = t[n]
	}

	return
}

// RSTtoABCSingle transfers a single point from (r,s,t) to (a,b,c) coordinates
func RSTtoABCSingle(r, s, t float64) (a, b, c float64) {
	const tol = 1e-12

	// Handle a coordinate
	if math.Abs(s+t) < tol {
		a = -1.0
	} else {
		a = 2.0*(1.0+r)/(-s-t) - 1.0
	}

	// Handle b coordinate
	if math.Abs(t-1.0) < tol {
		b = -1.0
	} else {
		b = 2.0*(1.0+s)/(1.0-t) - 1.0
	}

	c = t
	return
}

// isTetrahedronVertex checks if a point is one of the 4 vertices of the reference tetrahedron
func isTetrahedronVertex(r, s, t float64) bool {
	const tol = 1e-12
	// Vertex 1: (-1, -1, -1)
	if math.Abs(r+1) < tol && math.Abs(s+1) < tol && math.Abs(t+1) < tol {
		return true
	}
	// Vertex 2: (1, -1, -1)
	if math.Abs(r-1) < tol && math.Abs(s+1) < tol && math.Abs(t+1) < tol {
		return true
	}
	// Vertex 3: (-1, 1, -1)
	if math.Abs(r+1) < tol && math.Abs(s-1) < tol && math.Abs(t+1) < tol {
		return true
	}
	// Vertex 4: (-1, -1, 1)
	if math.Abs(r+1) < tol && math.Abs(s+1) < tol && math.Abs(t-1) < tol {
		return true
	}
	return false
}

// ABCtoRST is the inverse transformation from (a,b,c) to (r,s,t) coordinates
func ABCtoRST(a, b, c []float64) (r, s, t []float64) {
	Np := len(a)
	r = make([]float64, Np)
	s = make([]float64, Np)
	t = make([]float64, Np)

	for n := 0; n < Np; n++ {
		t[n] = c[n]
		s[n] = 0.5*(1+b[n])*(1-c[n]) - 1.0
		r[n] = 0.5*(1+a[n])*(1-b[n])*(1-c[n]) - 1.0
	}
	return
}

// ABCtoRSTSingle is the inverse transformation for a single point
func ABCtoRSTSingle(a, b, c float64) (r, s, t float64) {
	t = c
	s = 0.5*(1+b)*(1-c) - 1.0
	r = 0.5*(1+a)*(1-b)*(1-c) - 1.0
	return
}

// IsValidRST checks if (r,s,t) coordinates are within the reference tetrahedron
func IsValidRST(r, s, t float64) bool {
	const tol = 1e-10
	return r >= -1-tol && s >= -1-tol && t >= -1-tol && r+s+t <= -1+tol
}

// IsValidABC checks if (a,b,c) coordinates are within the valid range
func IsValidABC(a, b, c float64) bool {
	const tol = 1e-10
	return math.Abs(a) <= 1+tol && math.Abs(b) <= 1+tol && math.Abs(c) <= 1+tol
}
