package utils

import (
	"fmt"
	"math"
	"testing"
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

func TestSparseBlockMatrixMul(t *testing.T) {
	// Build a small sparse block matrix with 2 block rows and 2 block columns.
	// Each block is 2x2.
	nBlocks := 2
	blockRows, blockCols := 2, 2
	// We'll construct a block diagonal matrix for simplicity.
	// CSR: rowPtr = [0,1,2] means block row 0 has 1 block, block row 1 has 1 block.
	rowPtr := []int{0, 1, 2}
	colIdx := []int{0, 1} // block row 0: block in column 0; block row 1: block in column 1.

	// Create two blocks:
	// Block0: 2x2 identity.
	block0 := NewMatrix(blockRows, blockCols, []float64{1, 0, 0, 1})
	// Block1: 2x2 matrix [2,3;4,5].
	block1 := NewMatrix(blockRows, blockCols, []float64{2, 3, 4, 5})
	blocks := []Matrix{block0, block1}

	// Build the SparseBlockMatrix.
	sbm := SparseBlockMatrix{
		NrBlocks:  nBlocks,
		NcBlocks:  nBlocks,
		BlockRows: blockRows,
		BlockCols: blockCols,
		RowPtr:    rowPtr,
		ColIdx:    colIdx,
		Blocks:    blocks,
	}

	// Create a block vector x with dimension (NcBlocks*BlockCols) x 1.
	// Let x be [1,2,3,4]^T, interpreted as two blocks: [1,2] and [3,4].
	xData := []float64{1, 2, 3, 4}
	x := NewMatrix(nBlocks*blockCols, 1, xData)

	// Multiply: y = sbm * x.
	// Expected:
	// For block row 0: identity * [1,2] = [1,2].
	// For block row 1: [2,3;4,5] * [3,4] = [2*3+3*4, 4*3+5*4] = [18,32].
	y := sbm.Mul(x)

	expected := []float64{1, 2, 18, 32}
	if y.M.RawMatrix().Rows != nBlocks*blockRows {
		t.Errorf("Result vector has wrong number of rows")
	}
	for i := 0; i < len(expected); i++ {
		got := y.M.At(i, 0)
		if !almostEqual(got, expected[i], 1e-6) {
			t.Errorf("SparseBlockMatrix Mul: at index %d, got %v, expected %v", i, got, expected[i])
		}
	}
}

// TestSparseBlockMatrixPattern creates a 5x5 sparse block matrix (with each block 2x2),
// verifies that each block is correctly stored, and multiplies by a block vector.
func TestSparseBlockMatrixPattern(t *testing.T) {
	// --- Setup parameters ---
	nBlocks := 5                 // 5 block rows and 5 block columns.
	blockRows, blockCols := 2, 2 // Each block is 2x2.
	// Each block row has 1 diagonal + 3 off-diagonals = 4 nonzero blocks.
	totalNonzeroBlocks := nBlocks * 4

	// Build the CSR row pointer.
	rowPtr := make([]int, nBlocks+1)
	for i := 0; i <= nBlocks; i++ {
		rowPtr[i] = i * 4
	}

	// Build the CSR column indices.
	// For each block row i: first entry is i (diagonal),
	// then the first three indices from {0,1,2,3,4} excluding i.
	colIdx := make([]int, 0, totalNonzeroBlocks)
	for i := 0; i < nBlocks; i++ {
		colIdx = append(colIdx, i) // diagonal
		count := 0
		for j := 0; j < nBlocks && count < 3; j++ {
			if j == i {
				continue
			}
			colIdx = append(colIdx, j)
			count++
		}
	}
	if len(colIdx) != totalNonzeroBlocks {
		t.Fatalf("expected %d column indices, got %d", totalNonzeroBlocks, len(colIdx))
	}

	// Create the Blocks slice.
	// For each block row:
	//   Diagonal block: defined as [d,0;0,d] with d = 10*(i+1).
	//   Off-diagonals: for each off-diagonal at column j, candidate block:
	//         [10*j, 10*j+1; 10*j+2, 10*j+3].
	blocks := make([]Matrix, 0, totalNonzeroBlocks)
	for i := 0; i < nBlocks; i++ {
		// Diagonal block.
		d := 10.0 * float64(i+1)
		diagBlock := NewMatrix(blockRows, blockCols, []float64{
			d, 0,
			0, d,
		})
		blocks = append(blocks, diagBlock)
		// Off-diagonals.
		count := 0
		for j := 0; j < nBlocks && count < 3; j++ {
			if j == i {
				continue
			}
			candidate := NewMatrix(blockRows, blockCols, []float64{
				10 * float64(j), 10*float64(j) + 1,
				10*float64(j) + 2, 10*float64(j) + 3,
			})
			blocks = append(blocks, candidate)
			count++
		}
	}
	if len(blocks) != totalNonzeroBlocks {
		t.Fatalf("expected %d blocks, got %d", totalNonzeroBlocks, len(blocks))
	}

	// Build the SparseBlockMatrix.
	sbm := SparseBlockMatrix{
		NrBlocks:  nBlocks,
		NcBlocks:  nBlocks,
		BlockRows: blockRows,
		BlockCols: blockCols,
		RowPtr:    rowPtr,
		ColIdx:    colIdx,
		Blocks:    blocks,
	}

	// --- Validate addressing ---
	// Check that each block's (0,0) entry is as expected.
	for i := 0; i < nBlocks; i++ {
		start := sbm.RowPtr[i]
		end := sbm.RowPtr[i+1]
		for k := start; k < end; k++ {
			col := sbm.ColIdx[k]
			block := sbm.Blocks[k]
			val := block.M.At(0, 0)
			if k == start { // diagonal block.
				expected := 10.0 * float64(i+1)
				if !almostEqual(val, expected, 1e-6) {
					t.Errorf("Row %d, diagonal block: expected (0,0) = %v, got %v", i, expected, val)
				}
			} else {
				expected := 10.0 * float64(col)
				if !almostEqual(val, expected, 1e-6) {
					t.Errorf("Row %d, off-diagonal block (col %d): expected (0,0) = %v, got %v", i, col, expected, val)
				}
			}
		}
	}

	// --- Build the block vector x ---
	// x consists of 5 blocks; each block is a 2x2 matrix.
	// Overall, x has dimensions (nBlocks*blockCols) x p.
	p := 2
	totalRows := nBlocks * blockCols // 5*2 = 10 rows.
	xData := make([]float64, totalRows*p)
	// For each block j, define:
	// x_j = [ j, j+1; j+2, j+3 ]
	for j := 0; j < nBlocks; j++ {
		r := j * blockCols
		xData[r*p+0] = float64(j)
		xData[r*p+1] = float64(j + 1)
		xData[(r+1)*p+0] = float64(j + 2)
		xData[(r+1)*p+1] = float64(j + 3)
	}
	x := NewMatrix(totalRows, p, xData)

	// --- Multiply: y = A * x ---
	// y will have dimensions (NrBlocks*BlockRows) x p, i.e. 10x2.
	y := sbm.Mul(x)

	// --- Verify result for block row 0 ---
	// Block row 0 has 4 blocks:
	// Diagonal (col 0): D0 = [10,0;0,10]
	// Off-diagonals for row 0 are for columns: 1, 2, and 3.
	// Their candidate blocks are:
	// For col 1: C1 = [10,11;12,13]
	// For col 2: C2 = [20,21;22,23]
	// For col 3: C3 = [30,31;32,33]
	// Also, the corresponding x blocks are:
	// x0 = [0,1;2,3], x1 = [1,2;3,4], x2 = [2,3;4,5], x3 = [3,4;5,6].
	// Compute:
	// D0*x0 = [10,0;0,10]*[0,1;2,3] = [0,10;20,30].
	// C1*x1 = [10,11;12,13]*[1,2;3,4] = [10+33,20+44;12+39,24+52] = [43,64;51,76].
	// C2*x2 = [20,21;22,23]*[2,3;4,5] = [40+84,60+105;44+92,66+115] = [124,165;136,181].
	// C3*x3 = [30,31;32,33]*[3,4;5,6] = [90+155,120+186;96+165,128+198] = [245,306;261,326].
	// Sum these results:
	// Top-left: 0 + 43 + 124 + 245 = 412.
	// Top-right: 10 + 64 + 165 + 306 = 545.
	// Bottom-left: 20 + 51 + 136 + 261 = 468.
	// Bottom-right: 30 + 76 + 181 + 326 = 613.
	expectedBlock0 := [][]float64{
		{412, 545},
		{468, 613},
	}
	// Extract block row 0 of y: rows 0..1.
	yBlock0 := y.SubMatrix(0, 0, blockRows, p)
	for i := 0; i < blockRows; i++ {
		for j := 0; j < p; j++ {
			got := yBlock0.M.At(i, j)
			if !almostEqual(got, expectedBlock0[i][j], 1e-6) {
				t.Errorf("y block0 element (%d,%d): got %v, expected %v", i, j, got, expectedBlock0[i][j])
			}
		}
	}
}

// TestAccessHelpers creates a 5x5 sparse block matrix (with each block 2x2)
// using a pattern: each block row has a diagonal and 3 off-diagonals chosen from
// the set of 5 candidate blocks (excluding the diagonal).
// It then uses GetBlock/SetBlock to write and read back values and verifies
// that a subsequent matrix/vector multiply produces the expected result.
func TestAccessHelpers(t *testing.T) {
	// --- Parameters: 5 block rows/columns, each block 2x2.
	nBlocks := 5
	blockRows, blockCols := 2, 2
	// Each block row has 1 diagonal + 3 off-diagonals = 4 blocks per row.
	totalNonzeroBlocks := nBlocks * 4

	// Build the CSR row pointer.
	rowPtr := make([]int, nBlocks+1)
	for i := 0; i <= nBlocks; i++ {
		rowPtr[i] = i * 4
	}

	// Build the CSR column indices.
	// For each block row i: first entry is the diagonal (i),
	// then choose the first three indices from {0,1,2,3,4} excluding i.
	colIdx := make([]int, 0, totalNonzeroBlocks)
	for i := 0; i < nBlocks; i++ {
		colIdx = append(colIdx, i) // diagonal
		count := 0
		for j := 0; j < nBlocks && count < 3; j++ {
			if j == i {
				continue
			}
			colIdx = append(colIdx, j)
			count++
		}
	}

	// Allocate blocks.
	// For a diagonal block in block row i, define it as:
	//   D_i = [10*(i+1), 0; 0, 10*(i+1)]
	// For an off-diagonal block corresponding to column j:
	//   C_j = [10*j, 10*j+1; 10*j+2, 10*j+3]
	blocks := make([]Matrix, 0, totalNonzeroBlocks)
	for i := 0; i < nBlocks; i++ {
		// Diagonal block.
		d := 10.0 * float64(i+1)
		diagBlock := NewMatrix(blockRows, blockCols, []float64{
			d, 0,
			0, d,
		})
		blocks = append(blocks, diagBlock)
		// Off-diagonals.
		count := 0
		for j := 0; j < nBlocks && count < 3; j++ {
			if j == i {
				continue
			}
			offBlock := NewMatrix(blockRows, blockCols, []float64{
				10 * float64(j), 10*float64(j) + 1,
				10*float64(j) + 2, 10*float64(j) + 3,
			})
			blocks = append(blocks, offBlock)
			count++
		}
	}

	// Build the SparseBlockMatrix.
	sbm := SparseBlockMatrix{
		NrBlocks:  nBlocks,
		NcBlocks:  nBlocks,
		BlockRows: blockRows,
		BlockCols: blockCols,
		RowPtr:    rowPtr,
		ColIdx:    colIdx,
		Blocks:    blocks,
	}

	// --- Test GetBlock and SetBlock ---
	// For block row 2, block column 3 (an off-diagonal location), retrieve the block.
	block, found := sbm.GetBlock(2, 3)
	if !found {
		t.Fatalf("GetBlock: block (2,3) not found")
	}
	// Check that the (0,0) element is as expected:
	// For off-diagonals, we defined the candidate as [10*j, ...] so expected is 10*3 = 30.
	if !almostEqual(block.M.At(0, 0), 30, 1e-6) {
		t.Errorf("GetBlock: expected (0,0) element 30, got %v", block.M.At(0, 0))
	}

	// Now update that block using SetBlock.
	// Write a new 2x2 block: [100,101;102,103]
	newBlock := NewMatrix(blockRows, blockCols, []float64{
		100, 101,
		102, 103,
	})
	if err := sbm.SetBlock(2, 3, newBlock); err != nil {
		t.Fatalf("SetBlock failed: %v", err)
	}
	// Retrieve again and verify the change.
	updatedBlock, found := sbm.GetBlock(2, 3)
	if !found {
		t.Fatalf("GetBlock: block (2,3) not found after update")
	}
	if !almostEqual(updatedBlock.M.At(0, 0), 100, 1e-6) ||
		!almostEqual(updatedBlock.M.At(1, 1), 103, 1e-6) {
		t.Errorf("SetBlock: update failed, got block = %v", updatedBlock.Print(""))
	}

	// --- Test a Matrix/Vector Multiply ---
	// Build a block vector x with 5 blocks, each 2x2.
	// Let x_j = [ j, j+1; j+2, j+3 ].
	p := 2 // number of columns in x.
	totalRows := nBlocks * blockCols
	xData := make([]float64, totalRows*p)
	for j := 0; j < nBlocks; j++ {
		r := j * blockCols
		xData[r*p+0] = float64(j)
		xData[r*p+1] = float64(j + 1)
		xData[(r+1)*p+0] = float64(j + 2)
		xData[(r+1)*p+1] = float64(j + 3)
	}
	x := NewMatrix(totalRows, p, xData)

	// Multiply: y = A * x.
	y := sbm.Mul(x)

	// For this test, we simply print out the first block of y.
	// (In your real tests you would compare against an expected result.)
	yBlock0 := y.SubMatrix(0, 0, blockRows, p)
	t.Logf("Block row 0 of y:\n%v", yBlock0.Print(""))
	// For example, if you had manually computed the expected result for block row 0,
	// you could compare each element here.
}

func TestNewMatrixFromData(t *testing.T) {
	// Preallocate a large buffer for a 2x2 block.
	data := make([]float64, 4)
	// Fill the buffer with known values.
	data[0] = 1.1
	data[1] = 1.2
	data[2] = 2.1
	data[3] = 2.2

	// Create a Matrix that uses this buffer.
	m := NewMatrixFromData(2, 2, data)
	fmt.Println("Matrix m:")
	m.M.Apply(func(i, j int, v float64) float64 {
		return v
	}, m.M)
	// You can use your existing Print method instead.
	m.Print("Matrix m")

	// Now modify the underlying slice.
	data[0] = 9.9
	data[3] = 8.8
	// The change is visible in m.
	m.Print("Matrix m after modifying underlying data")
}
