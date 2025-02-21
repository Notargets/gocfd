package utils

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

// almostEqual returns true if a and b differ by less than tol.
func almostEqual(a, b, tol float64) bool {
	return math.Abs(a-b) < tol
}

func TestNewMatrix(t *testing.T) {
	// Create a 2x3 matrix using NewMatrix.
	m := NewMatrix(2, 3)
	r, c := m.Dims()
	if r != 2 || c != 3 {
		t.Errorf("NewMatrix dims expected (2,3), got (%d,%d)", r, c)
	}
	if len(m.DataP) < 6 {
		t.Errorf("NewMatrix underlying data length expected at least 6, got %d", len(m.DataP))
	}
}

func TestTranspose(t *testing.T) {
	// Create a 2x3 matrix and fill it with known values.
	m := NewMatrix(2, 3)
	m.M.Set(0, 0, 1)
	m.M.Set(0, 1, 2)
	m.M.Set(0, 2, 3)
	m.M.Set(1, 0, 4)
	m.M.Set(1, 1, 5)
	m.M.Set(1, 2, 6)

	// Transpose m.
	tm := m.Transpose()
	r, c := tm.Dims()
	if r != 3 || c != 2 {
		t.Errorf("Transpose dims expected (3,2), got (%d,%d)", r, c)
	}
	// Expected transpose is:
	// [1,4]
	// [2,5]
	// [3,6]
	expected := [][]float64{
		{1, 4},
		{2, 5},
		{3, 6},
	}
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			got := tm.M.At(i, j)
			if !almostEqual(got, expected[i][j], 1e-6) {
				t.Errorf("Transpose: at (%d,%d), got %v, want %v", i, j, got, expected[i][j])
			}
		}
	}
}

func TestMul(t *testing.T) {
	// Test multiplication: m (2x3) * A (3x2) should yield a 2x2 matrix.
	m := NewMatrix(2, 3)
	A := NewMatrix(3, 2)
	// Fill m with [1,2,3; 4,5,6]
	m.M.Set(0, 0, 1)
	m.M.Set(0, 1, 2)
	m.M.Set(0, 2, 3)
	m.M.Set(1, 0, 4)
	m.M.Set(1, 1, 5)
	m.M.Set(1, 2, 6)
	// Fill A with [7,8; 9,10; 11,12]
	A.M.Set(0, 0, 7)
	A.M.Set(0, 1, 8)
	A.M.Set(1, 0, 9)
	A.M.Set(1, 1, 10)
	A.M.Set(2, 0, 11)
	A.M.Set(2, 1, 12)

	prod := m.Mul(A)
	// Expected product is:
	// [1*7+2*9+3*11, 1*8+2*10+3*12] = [58, 64]
	// [4*7+5*9+6*11, 4*8+5*10+6*12] = [139, 154]
	expected := [][]float64{
		{58, 64},
		{139, 154},
	}
	r, c := prod.Dims()
	if r != 2 || c != 2 {
		t.Errorf("Mul dims expected (2,2), got (%d,%d)", r, c)
	}
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			got := prod.M.At(i, j)
			if !almostEqual(got, expected[i][j], 1e-6) {
				t.Errorf("Mul: at (%d,%d), got %v, want %v", i, j, got, expected[i][j])
			}
		}
	}
}

func TestAdd(t *testing.T) {
	// Test addition: add two 2x2 matrices.
	m := NewMatrix(2, 2)
	A := NewMatrix(2, 2)
	// m = [1,2; 3,4]
	m.M.Set(0, 0, 1)
	m.M.Set(0, 1, 2)
	m.M.Set(1, 0, 3)
	m.M.Set(1, 1, 4)
	// A = [5,6; 7,8]
	A.M.Set(0, 0, 5)
	A.M.Set(0, 1, 6)
	A.M.Set(1, 0, 7)
	A.M.Set(1, 1, 8)

	sum := m.Add(A)
	// Expected sum = [6,8;10,12]
	expected := [][]float64{
		{6, 8},
		{10, 12},
	}
	r, c := sum.Dims()
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			got := sum.M.At(i, j)
			if !almostEqual(got, expected[i][j], 1e-6) {
				t.Errorf("Add: at (%d,%d), got %v, want %v", i, j, got, expected[i][j])
			}
		}
	}
}

func TestSubMatrix(t *testing.T) {
	// Create a 4x4 matrix with values 1..16 (row-major).
	m := NewMatrix(4, 4)
	for i := 0; i < 16; i++ {
		m.DataP[i] = float64(i + 1)
	}
	// Extract a 2x2 submatrix starting at row 1, col 1.
	// For the matrix:
	//  1  2  3  4
	//  5  6  7  8
	//  9 10 11 12
	// 13 14 15 16
	// The submatrix should be:
	// [6,7;10,11]
	sub := m.SubMatrix(1, 1, 2, 2)
	expected := [][]float64{
		{6, 7},
		{10, 11},
	}
	r, c := sub.Dims()
	if r != 2 || c != 2 {
		t.Errorf("SubMatrix dims expected (2,2), got (%d,%d)", r, c)
	}
	for i := 0; i < r; i++ {
		for j := 0; j < c; j++ {
			got := sub.M.At(i, j)
			if !almostEqual(got, expected[i][j], 1e-6) {
				t.Errorf("SubMatrix: at (%d,%d), got %v, want %v", i, j, got, expected[i][j])
			}
		}
	}
}

func TestNewMatrixFromData(t *testing.T) {
	tol := 0.0000001
	// Preallocate a large buffer for a 2x2 block.
	data := make([]float64, 4)
	// Fill the buffer with known values.
	data[0] = 1.1
	data[1] = 1.2
	data[2] = 2.1
	data[3] = 2.2

	// Create a Matrix that uses this buffer.
	m := NewMatrix(2, 2)
	_ = m.ResetView(data)
	if testing.Verbose() {
		m.Print("Matrix m")
	}
	assert.InDeltaSlicef(t, data, m.DataP, tol, "")

	// Now modify the underlying slice.
	data[0] = 9.9
	data[3] = 8.8
	// The change is visible in m.
	if testing.Verbose() {
		m.Print("Matrix m after modifying underlying data")
	}
	assert.InDeltaSlicef(t, data, m.DataP, tol, "")
	newMat := NewMatrix(2, 2, data)
	data[0] = 12
	data[1] = 13
	if testing.Verbose() {
		newMat.Print("New with old data")
	}
	assert.InDeltaSlicef(t, data, m.DataP, tol, "")
}

// TestGMRESIdentity constructs a simple block‐system with A = I and b given, and verifies
// that GMRES returns x such that ||b − A·x|| is below a tolerance.
func TestGMRESIdentity(t *testing.T) {
	// --- Create system matrix A as a block identity.
	// A is a 3×3 block matrix (each block is 2×2), with nonzero blocks only on the diagonal.
	addressesA := [][2]int{
		{0, 0}, {1, 1}, {2, 2},
	}
	A := NewBlockSparse(3, 3, 2, 2, addressesA)
	// Fill the allocated blocks in A with 2×2 identity.
	for _, addr := range addressesA {
		block := A.GetBlockView(addr[0], addr[1])
		block.M.Set(0, 0, 1.0)
		block.M.Set(0, 1, 0.0)
		block.M.Set(1, 0, 0.0)
		block.M.Set(1, 1, 1.0)
	}

	// --- Create right-hand side vector b as a dense block vector.
	// b is a 3×1 block matrix (each block is 2×2). For a dense block vector we allocate every row.
	addressesB := make([][2]int, 3)
	for i := 0; i < 3; i++ {
		addressesB[i] = [2]int{i, 0}
	}
	b := NewBlockSparse(3, 1, 2, 2, addressesB)
	// Fill b with known values.
	// Let block 0 = [1,2;3,4], block 1 = [5,6;7,8], block 2 = [9,10;11,12]
	valuesB := [][]float64{
		{1, 2, 3, 4},
		{5, 6, 7, 8},
		{9, 10, 11, 12},
	}
	for i := 0; i < 3; i++ {
		block := b.GetBlockView(i, 0)
		// Fill in row-major order.
		for j, v := range valuesB[i] {
			row := j / 2
			col := j % 2
			block.M.Set(row, col, v)
		}
	}

	// --- Run GMRES to solve A x = b.
	// Since A is the identity, the solution should equal b.
	tol := 1e-6
	maxIter := 5
	x := A.GMRES(b, tol, maxIter)

	// --- Compute the residual: r = b - A*x.
	Ax := A.Mul(x)
	rPool := b.Subtract(Ax)
	resNorm := rPool.FrobNorm()

	// For the identity, the residual should be (nearly) zero.
	if resNorm > tol {
		t.Errorf("GMRES residual norm too high: got %v, expected below %v", resNorm, tol)
	}
}

// TestSparseBlockMatrixSparseMul tests the multiplication of an "actually sparse" block matrix A
// by a dense block vector b. A is a 4x4 block matrix (blocks are 2x2) with nonzero blocks only at:
// (0,0), (0,2), (1,1), (2,0), (3,3). b is a dense block vector (4x1).
// The expected result is computed manually.
func TestSparseBlockMatrixSparseMul(t *testing.T) {
	// --- Create A: 4x4 block matrix with 2x2 blocks.
	// Nonzero blocks allocated at: (0,0), (0,2), (1,1), (2,0), (3,3).
	addressesA := [][2]int{
		{0, 0},
		{0, 2},
		{1, 1},
		{2, 0},
		{3, 3},
	}
	A := NewBlockSparse(4, 4, 2, 2, addressesA)

	// Fill A's allocated blocks with known values.
	// GetBlockView (0,0): Identity.
	block := A.GetBlockView(0, 0)
	block.M.Set(0, 0, 1.0)
	block.M.Set(0, 1, 0.0)
	block.M.Set(1, 0, 0.0)
	block.M.Set(1, 1, 1.0)

	// GetBlockView (0,2): Constant matrix with all entries = 2.
	block = A.GetBlockView(0, 2)
	block.M.Set(0, 0, 2.0)
	block.M.Set(0, 1, 2.0)
	block.M.Set(1, 0, 2.0)
	block.M.Set(1, 1, 2.0)

	// GetBlockView (1,1): Identity.
	block = A.GetBlockView(1, 1)
	block.M.Set(0, 0, 1.0)
	block.M.Set(0, 1, 0.0)
	block.M.Set(1, 0, 0.0)
	block.M.Set(1, 1, 1.0)

	// GetBlockView (2,0): Constant matrix with all entries = 3.
	block = A.GetBlockView(2, 0)
	block.M.Set(0, 0, 3.0)
	block.M.Set(0, 1, 3.0)
	block.M.Set(1, 0, 3.0)
	block.M.Set(1, 1, 3.0)

	// GetBlockView (3,3): Identity.
	block = A.GetBlockView(3, 3)
	block.M.Set(0, 0, 1.0)
	block.M.Set(0, 1, 0.0)
	block.M.Set(1, 0, 0.0)
	block.M.Set(1, 1, 1.0)

	// --- Create b: Dense block vector of dimension 4x1 (each block 2x2).
	addressesB := make([][2]int, 4)
	for i := 0; i < 4; i++ {
		addressesB[i] = [2]int{i, 0}
	}
	b := NewBlockSparse(4, 1, 2, 2, addressesB)
	// Fill b with known values:
	// GetBlockView (0,0): [1,1;1,1]
	block = b.GetBlockView(0, 0)
	block.M.Set(0, 0, 1.0)
	block.M.Set(0, 1, 1.0)
	block.M.Set(1, 0, 1.0)
	block.M.Set(1, 1, 1.0)
	// GetBlockView (1,0): [2,2;2,2]
	block = b.GetBlockView(1, 0)
	block.M.Set(0, 0, 2.0)
	block.M.Set(0, 1, 2.0)
	block.M.Set(1, 0, 2.0)
	block.M.Set(1, 1, 2.0)
	// GetBlockView (2,0): [3,3;3,3]
	block = b.GetBlockView(2, 0)
	block.M.Set(0, 0, 3.0)
	block.M.Set(0, 1, 3.0)
	block.M.Set(1, 0, 3.0)
	block.M.Set(1, 1, 3.0)
	// GetBlockView (3,0): [4,4;4,4]
	block = b.GetBlockView(3, 0)
	block.M.Set(0, 0, 4.0)
	block.M.Set(0, 1, 4.0)
	block.M.Set(1, 0, 4.0)
	block.M.Set(1, 1, 4.0)

	// --- Multiply: result = A * b.
	result := A.Mul(b)

	// --- Expected result:
	// r[0] = (0,0)*b[0] + (0,2)*b[2]
	//      = I * [1,1;1,1] + [2,2;2,2] * [3,3;3,3]
	//      = [1,1;1,1] + [12,12;12,12] = [13,13;13,13]
	//
	// r[1] = (1,1)*b[1] = I * [2,2;2,2] = [2,2;2,2]
	//
	// r[2] = (2,0)*b[0] = [3,3;3,3] * [1,1;1,1] = [6,6;6,6]
	//
	// r[3] = (3,3)*b[3] = I * [4,4;4,4] = [4,4;4,4]
	expected := map[[2]int][]float64{
		[2]int{0, 0}: {13, 13, 13, 13},
		[2]int{1, 0}: {2, 2, 2, 2},
		[2]int{2, 0}: {6, 6, 6, 6},
		[2]int{3, 0}: {4, 4, 4, 4},
	}
	tol := 1e-6
	// Verify each expected block in the result.
	for key, expVals := range expected {
		rblock := result.GetBlockView(key[0], key[1])
		for i := 0; i < 2; i++ {
			for j := 0; j < 2; j++ {
				got := rblock.M.At(i, j)
				exp := expVals[i*2+j]
				if math.Abs(got-exp) > tol {
					t.Errorf("Result block %v at (%d,%d): got %v, expected %v", key, i, j, got, exp)
				}
			}
		}
	}
}

// TestGMRESBlockTridiagonalKnownSolution creates a 5x5 block tridiagonal system where:
//   - Each block of A is 2x2.
//   - Diagonal blocks are 4·I and off-diagonals (immediately above and below) are –I.
//
// The true solution x_true is chosen as a block vector with every block equal to a 2x1 vector of ones.
// Then b = A*x_true is computed, GMRES is used to solve A*x = b, and the computed x is compared against x_true.
func TestGMRESBlockTridiagonalKnownSolution(t *testing.T) {
	// --- Dimensions ---
	numBlockRows := 5 // A is 5x5 in block dimensions.
	numBlockCols := 5
	// A's blocks are 2x2.
	A_blockRows, A_blockCols := 2, 2
	// We want the dense vector x_true (and b) to be block vectors with blocks of size 2x1.
	x_blockRows, x_blockCols := 2, 1

	// --- Build addresses for A's nonzero blocks in a block-tridiagonal pattern.
	// For block row i, allocate:
	//    if i == 0: (0,0) and (0,1)
	//    if 0 < i < numBlockRows-1: (i,i-1), (i,i), (i,i+1)
	//    if i == numBlockRows-1: (i,i-1) and (i,i)
	var addressesA [][2]int
	for i := 0; i < numBlockRows; i++ {
		if i == 0 {
			addressesA = append(addressesA, [2]int{i, i})     // diagonal
			addressesA = append(addressesA, [2]int{i, i + 1}) // upper
		} else if i == numBlockRows-1 {
			addressesA = append(addressesA, [2]int{i, i - 1}) // lower
			addressesA = append(addressesA, [2]int{i, i})     // diagonal
		} else {
			addressesA = append(addressesA, [2]int{i, i - 1}) // lower
			addressesA = append(addressesA, [2]int{i, i})     // diagonal
			addressesA = append(addressesA, [2]int{i, i + 1}) // upper
		}
	}

	// Create A as a BlockPool.
	A := NewBlockSparse(numBlockRows, numBlockCols, A_blockRows, A_blockCols,
		addressesA)

	// Fill in A's allocated blocks.
	// For each allocated block at (i,j):
	//   if i == j (diagonal): 4*I.
	//   if |i - j| == 1 (off-diagonal): -I.
	for _, addr := range addressesA {
		i, j := addr[0], addr[1]
		block := A.GetBlockView(i, j)
		if i == j {
			// Diagonal: 4*I.
			block.M.Set(0, 0, 4)
			block.M.Set(0, 1, 0)
			block.M.Set(1, 0, 0)
			block.M.Set(1, 1, 4)
		} else {
			// Off-diagonal: -I.
			block.M.Set(0, 0, -1)
			block.M.Set(0, 1, 0)
			block.M.Set(1, 0, 0)
			block.M.Set(1, 1, -1)
		}
	}

	// --- Create x_true as a dense block vector.
	// x_true is represented as a BlockPool with numBlockRows block rows and 1 block column.
	// Each block has dimensions x_blockRows x x_blockCols (2x1).
	var addressesX [][2]int
	for i := 0; i < numBlockRows; i++ {
		addressesX = append(addressesX, [2]int{i, 0})
	}
	x_true := NewBlockSparse(numBlockRows, 1, x_blockRows, x_blockCols, addressesX)
	// Fill each block of x_true with ones.
	for i := 0; i < numBlockRows; i++ {
		block := x_true.GetBlockView(i, 0)
		// For a 2x1 block.
		block.M.Set(0, 0, 1)
		block.M.Set(1, 0, 1)
	}

	// --- Compute b = A * x_true.
	// Note: Multiplication here involves A (5x5, blocks 2x2) times x_true (5x1, blocks 2x1).
	// The result b will be a BlockPool representing a dense block vector with block size 2x1.
	b := A.Mul(x_true)

	// --- Run GMRES to solve A * x = b.
	tol := 1e-6
	maxIter := 50
	x_approx := A.GMRES(b, tol, maxIter)

	// --- Compare x_approx to x_true.
	// We can compute the block-wise error.
	errSum := 0.0
	totalElements := 0
	for i := 0; i < numBlockRows; i++ {
		x_trueBlock := x_true.GetBlockView(i, 0)
		x_approxBlock := x_approx.GetBlockView(i, 0)
		r, c := x_trueBlock.Dims()
		for p := 0; p < r; p++ {
			for q := 0; q < c; q++ {
				diff := x_trueBlock.M.At(p, q) - x_approxBlock.M.At(p, q)
				errSum += diff * diff
				totalElements++
			}
		}
	}
	avgErr := math.Sqrt(errSum / float64(totalElements))
	if avgErr > tol {
		t.Errorf("GMRES solution error too high: avg error %v, tol %v", avgErr, tol)
	}
}
