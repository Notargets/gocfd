package utils

import (
	"math"
	"sort"

	"gonum.org/v1/gonum/mat"
)

// Add this method to your utils.Matrix type

func (m Matrix) Cond() float64 {
	return m.ConditionNumber()
}

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

// FrobeniusNorm computes the Frobenius norm of the matrix
// The Frobenius norm is defined as: ||A||_F = sqrt(sum(a_ij^2))
// It's equivalent to the Euclidean norm of the matrix when viewed as a vector
func (m Matrix) FrobeniusNorm() float64 {
	sumSquares := 0.0

	// If your Matrix type has direct access to the underlying data array
	for _, val := range m.DataP {
		sumSquares += val * val
	}

	return math.Sqrt(sumSquares)
}

// FrobeniusNormSquared returns the square of the Frobenius norm
// This is more efficient when you only need to compare norms
func (m Matrix) FrobeniusNormSquared() float64 {
	sumSquares := 0.0

	for _, val := range m.DataP {
		sumSquares += val * val
	}

	return sumSquares
}

// For sparse matrices or when you want to compute weighted Frobenius norm
func (m Matrix) WeightedFrobeniusNorm(weights Matrix) float64 {
	rows, cols := m.Dims()
	wRows, wCols := weights.Dims()

	if rows != wRows || cols != wCols {
		panic("Matrix dimensions must match for weighted Frobenius norm")
	}

	sumSquares := 0.0

	for i, val := range m.DataP {
		weight := weights.DataP[i]
		sumSquares += weight * val * val
	}

	return math.Sqrt(sumSquares)
}

// Utility function to compute the Frobenius norm of the difference between two matrices
func (m Matrix) FrobeniusDistance(B Matrix) float64 {
	rowsA, colsA := m.Dims()
	rowsB, colsB := B.Dims()

	if rowsA != rowsB || colsA != colsB {
		panic("Matrix dimensions must match for Frobenius distance")
	}

	sumSquares := 0.0

	for i, val := range m.DataP {
		diff := val - B.DataP[i]
		sumSquares += diff * diff
	}

	return math.Sqrt(sumSquares)
}

// NormalizedFrobeniusNorm returns the Frobenius norm divided by sqrt(rows*cols)
// This gives a scale-invariant measure of the "average" matrix entry magnitude
func (m Matrix) NormalizedFrobeniusNorm() float64 {
	rows, cols := m.Dims()
	if rows == 0 || cols == 0 {
		return 0.0
	}

	return m.FrobeniusNorm() / math.Sqrt(float64(rows*cols))
}

// SVDResult holds the results of singular value decomposition
// A = U * S * V^T where U and V are orthogonal and S is diagonal
type SVDResult struct {
	U Matrix // Left singular vectors (m x m)
	S Vector // Singular values (min(m,n) x 1)
	V Matrix // Right singular vectors (n x n)
}

// SVD computes the singular value decomposition of the matrix
func (m Matrix) SVD() SVDResult {
	// Get dimensions
	rows, cols := m.Dims()

	// Create SVD struct from gonum
	var svd mat.SVD
	ok := svd.Factorize(m.M, mat.SVDFull)
	if !ok {
		panic("SVD factorization failed")
	}

	// Extract U matrix
	var u mat.Dense
	svd.UTo(&u)
	U := WrapDense(&u)

	// Extract V matrix (note: gonum returns V, not V^T)
	var v mat.Dense
	svd.VTo(&v)
	V := WrapDense(&v)

	// Extract singular values
	sVals := make([]float64, minInt(rows, cols))
	svd.Values(sVals)
	S := NewVector(len(sVals), sVals)

	return SVDResult{
		U: U,
		S: S,
		V: V,
	}
}

// SVDThin computes the thin/economy SVD
// For m x n matrix with m > n, returns U as m x n instead of m x m
func (m Matrix) SVDThin() SVDResult {
	// Get dimensions
	rows, cols := m.Dims()

	// Create SVD struct from gonum
	var svd mat.SVD
	ok := svd.Factorize(m.M, mat.SVDThin)
	if !ok {
		panic("SVD factorization failed")
	}

	// Extract U matrix (thin version)
	var u mat.Dense
	svd.UTo(&u)
	U := WrapDense(&u)

	// Extract V matrix
	var v mat.Dense
	svd.VTo(&v)
	V := WrapDense(&v)

	// Extract singular values
	sVals := make([]float64, minInt(rows, cols))
	svd.Values(sVals)
	S := NewVector(len(sVals), sVals)

	return SVDResult{
		U: U,
		S: S,
		V: V,
	}
}

// SingularValues returns just the singular values (more efficient if you don't need U and V)
func (m Matrix) SingularValues() []float64 {
	var svd mat.SVD
	ok := svd.Factorize(m.M, mat.SVDNone)
	if !ok {
		panic("SVD factorization failed")
	}

	rows, cols := m.Dims()
	sVals := make([]float64, minInt(rows, cols))
	svd.Values(sVals)

	return sVals
}

// Rank computes the numerical rank of the matrix using SVD
func (m Matrix) Rank(tolerance float64) int {
	sVals := m.SingularValues()

	// If tolerance is negative, use machine epsilon * max dimension * largest singular value
	if tolerance < 0 {
		rows, cols := m.Dims()
		maxDim := float64(maxInt(rows, cols))
		tolerance = maxDim * machineEpsilon * sVals[0]
	}

	rank := 0
	for _, s := range sVals {
		if s > tolerance {
			rank++
		} else {
			break
		}
	}

	return rank
}

// PseudoInverse computes the Moore-Penrose pseudoinverse using SVD
func (m Matrix) PseudoInverse(tolerance float64) Matrix {
	svd := m.SVD()
	rows, cols := m.Dims()

	// Compute reciprocal of singular values
	sInv := make([]float64, len(svd.S.DataP))

	// Use default tolerance if not specified
	if tolerance < 0 {
		maxDim := float64(maxInt(rows, cols))
		tolerance = maxDim * machineEpsilon * svd.S.AtVec(0)
	}

	for i, s := range svd.S.DataP {
		if s > tolerance {
			sInv[i] = 1.0 / s
		} else {
			sInv[i] = 0.0
		}
	}

	// Create diagonal matrix with reciprocal singular values
	SInv := NewMatrix(cols, rows)
	minDim := minInt(rows, cols)
	for i := 0; i < minDim; i++ {
		SInv.Set(i, i, sInv[i])
	}

	// Compute pseudoinverse: A+ = V * S+ * U^T
	temp := svd.V.Mul(SInv)
	return temp.Mul(svd.U.Transpose())
}

// NuclearNorm computes the nuclear norm (sum of singular values)
func (m Matrix) NuclearNorm() float64 {
	sVals := m.SingularValues()
	sum := 0.0
	for _, s := range sVals {
		sum += s
	}
	return sum
}

// SpectralNorm computes the spectral norm (largest singular value)
func (m Matrix) SpectralNorm() float64 {
	sVals := m.SingularValues()
	if len(sVals) == 0 {
		return 0.0
	}
	return sVals[0]
}

// Methods for SVDResult

// Reconstruct builds the original matrix from the SVD components
func (svd SVDResult) Reconstruct() Matrix {
	// Create diagonal matrix from singular values
	m := len(svd.U.DataP) / svd.U.Cols()
	n := len(svd.V.DataP) / svd.V.Cols()
	S := NewMatrix(m, n)

	for i := 0; i < len(svd.S.DataP); i++ {
		S.Set(i, i, svd.S.AtVec(i))
	}

	// Compute U * S * V^T
	temp := svd.U.Mul(S)
	return temp.Mul(svd.V.Transpose())
}

// LowRankApproximation returns the best rank-k approximation of the original matrix
func (svd SVDResult) LowRankApproximation(k int) Matrix {
	if k > len(svd.S.DataP) {
		k = len(svd.S.DataP)
	}

	// Get dimensions
	uRows, _ := svd.U.Dims()
	_, vCols := svd.V.Dims()

	// Extract first k columns of U
	Uk := NewMatrix(uRows, k)
	for i := 0; i < uRows; i++ {
		for j := 0; j < k; j++ {
			Uk.Set(i, j, svd.U.At(i, j))
		}
	}

	// Create diagonal matrix with first k singular values
	Sk := NewMatrix(k, k)
	for i := 0; i < k; i++ {
		Sk.Set(i, i, svd.S.AtVec(i))
	}

	// Extract first k columns of V
	Vk := NewMatrix(vCols, k)
	for i := 0; i < vCols; i++ {
		for j := 0; j < k; j++ {
			Vk.Set(i, j, svd.V.At(i, j))
		}
	}

	// Compute Uk * Sk * Vk^T
	temp := Uk.Mul(Sk)
	return temp.Mul(Vk.Transpose())
}

// Helper functions

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

const machineEpsilon = 2.220446049250313e-16

// WrapDense creates a utils.Matrix from a gonum mat.Dense
// You'll need to implement this based on your Matrix type definition
func WrapDense(d *mat.Dense) Matrix {
	rows, cols := d.Dims()
	data := make([]float64, rows*cols)

	// Copy data from dense matrix
	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			data[i*cols+j] = d.At(i, j)
		}
	}

	return NewMatrix(rows, cols, data)
}
