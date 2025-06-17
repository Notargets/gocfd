package gonudg

import (
	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

// JacobiP evaluates Jacobi polynomial of type (alpha, beta) at points x for order N
// This wraps the DG1D.JacobiP function to work with plain float64 slices
func JacobiP(x []float64, alpha, beta float64, N int) []float64 {
	// Convert to utils.Vector
	xVec := utils.NewVector(len(x), x)

	// Call DG1D function
	result := DG1D.JacobiP(xVec, alpha, beta, N)

	// Convert back to slice
	return result
}

// JacobiPSingle evaluates Jacobi polynomial at a single point
func JacobiPSingle(x, alpha, beta float64, N int) float64 {
	// Create single element vectors
	xVec := utils.NewVector(1, []float64{x})
	result := DG1D.JacobiP(xVec, alpha, beta, N)
	return result[0]
}

// GradJacobiP evaluates the derivative of Jacobi polynomial
// This wraps the DG1D.GradJacobiP function to work with plain float64 slices
func GradJacobiP(x []float64, alpha, beta float64, N int) []float64 {
	// Convert to utils.Vector
	xVec := utils.NewVector(len(x), x)

	// Call DG1D function
	result := DG1D.GradJacobiP(xVec, alpha, beta, N)

	// Convert back to slice
	return result
}

// GradJacobiPSingle evaluates the derivative of Jacobi polynomial at a single point
func GradJacobiPSingle(x, alpha, beta float64, N int) float64 {
	// Create single element vectors
	xVec := utils.NewVector(1, []float64{x})
	result := DG1D.GradJacobiP(xVec, alpha, beta, N)
	return result[0]
}

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

	// Use DG2D.JacobiGQ if available
	x, _ := DG1D.JacobiGQ(alpha, beta, N)
	return x.DataP
}
