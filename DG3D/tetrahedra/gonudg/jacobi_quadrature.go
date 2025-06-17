package gonudg

import (
	"github.com/notargets/gocfd/DG1D"
)

// JacobiGL computes the Gauss-Lobatto quadrature points for Jacobi polynomials
// These are the zeros of (1-x^2)*P'_N^{alpha,beta}(x)
func JacobiGL(alpha, beta float64, N int) []float64 {
	if N == 0 {
		return []float64{0.0}
	}

	if N == 1 {
		return []float64{-1.0, 1.0}
	}

	// Get the interior Gauss-Jacobi points
	xint := JacobiGQ(alpha+1, beta+1, N-1)

	// Combine with endpoints
	x := make([]float64, N+1)
	x[0] = -1.0
	copy(x[1:N], xint)
	x[N] = 1.0

	return x
}

// JacobiGQ computes the Gauss quadrature points for Jacobi polynomials
// These are the zeros of P_N^{alpha,beta}(x)
func JacobiGQ(alpha, beta float64, N int) []float64 {
	if N == 0 {
		return []float64{-(alpha - beta) / (alpha + beta + 2)}
	}

	// Use DG1D.JacobiGQ if available
	x, _ := DG1D.JacobiGQ(alpha, beta, N)
	return x.DataP
}
