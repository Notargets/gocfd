package gonudg

// INDEXING NOTE: Original C++ code uses 1-based indexing to emulate Matlab behavior.
// This Go port uses standard 0-based indexing. Example conversions:
//   C++: sk = 1; V3D(All,sk) = ...    ->    Go: sk = 0; V3D.SetCol(sk, ...)
//   C++: Fmask[1] (first face)        ->    Go: Fmask[0] (first face)
// The indexing has been correctly translated throughout this port.

import (
	"github.com/notargets/gocfd/DG1D"
)

// JacobiGL computes the Gauss-Lobatto quadrature points for Jacobi polynomials
// These are the zeros of (1-X^2)*P'_N^{alpha,beta}(X)
func JacobiGL(alpha, beta float64, N int) []float64 {
	if N == 0 {
		return []float64{0.0}
	}

	if N == 1 {
		return []float64{-1.0, 1.0}
	}

	// Get the interior Gauss-Jacobi points
	// FIXED: Use N-2 to match C++ implementation
	// This gives us N-1 interior points, plus 2 endpoints = N+1 total points
	xint := JacobiGQ(alpha+1, beta+1, N-2)

	// Combine with endpoints
	x := make([]float64, N+1)
	x[0] = -1.0
	copy(x[1:N], xint) // Copy N-1 interior points
	x[N] = 1.0

	return x
}

// JacobiGQ computes the Gauss quadrature points for Jacobi polynomials
// These are the zeros of P_N^{alpha,beta}(X)
func JacobiGQ(alpha, beta float64, N int) []float64 {
	if N == 0 {
		return []float64{-(alpha - beta) / (alpha + beta + 2)}
	}

	// Use DG1D.JacobiGQ if available
	x, _ := DG1D.JacobiGQ(alpha, beta, N)
	return x.DataP
}
