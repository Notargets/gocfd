package utils

import (
	"math"
	"sort"

	"gonum.org/v1/gonum/mat"
)

// Add this method to your utils.Matrix type

func (m Matrix) ConditionNumber() float64 {
	// Get the underlying mat.Dense
	dense := m.M // Assuming M is the mat.Dense field

	// Compute SVD to get singular values
	var svd mat.SVD
	if !svd.Factorize(dense, mat.SVDThin) {
		// If SVD fails, return a large number indicating poor conditioning
		return 1e16
	}

	// Get singular values
	values := svd.Values(nil)

	if len(values) == 0 {
		return 1e16
	}

	// Condition number is ratio of largest to smallest singular value
	minVal := values[len(values)-1] // Singular values are in descending order
	maxVal := values[0]

	// Handle near-zero singular values
	if minVal < 1e-16 {
		return 1e16
	}

	return maxVal / minVal
}

// Alternative implementation using QR decomposition (faster for square matrices)
func (m Matrix) ConditionNumberQR() float64 {
	dense := m.M
	rows, cols := dense.Dims()

	// Only works for square matrices
	if rows != cols {
		return m.ConditionNumber() // Fall back to SVD
	}

	// Compute QR decomposition
	var qr mat.QR
	qr.Factorize(dense)

	// Extract R matrix
	var r mat.Dense
	qr.RTo(&r)

	// Condition number approximation from R diagonal elements
	var minDiag, maxDiag float64 = math.Inf(1), 0

	for i := 0; i < rows; i++ {
		diagVal := math.Abs(r.At(i, i))
		if diagVal > maxDiag {
			maxDiag = diagVal
		}
		if diagVal < minDiag {
			minDiag = diagVal
		}
	}

	if minDiag < 1e-16 {
		return 1e16
	}

	return maxDiag / minDiag
}

// For getting eigenvalues (needed for mass matrix analysis)
func (m Matrix) Eigenvalues() []float64 {
	dense := m.M
	rows, cols := dense.Dims()

	if rows != cols {
		panic("Eigenvalues only defined for square matrices")
	}

	// Use eigen decomposition
	var eigen mat.Eigen
	if !eigen.Factorize(dense, mat.EigenRight) {
		return nil
	}

	// Get eigenvalues (real parts only for now)
	values := eigen.Values(nil)
	realValues := make([]float64, len(values))

	for i, val := range values {
		realValues[i] = real(val)
	}

	// Sort eigenvalues
	sort.Float64s(realValues)

	return realValues
}

// For getting singular values (useful for debugging)
func (m Matrix) SingularValues() (min, max float64) {
	dense := m.M

	var svd mat.SVD
	if !svd.Factorize(dense, mat.SVDThin) {
		return 0, 1e16
	}

	values := svd.Values(nil)
	if len(values) == 0 {
		return 0, 1e16
	}

	return values[len(values)-1], values[0]
}
