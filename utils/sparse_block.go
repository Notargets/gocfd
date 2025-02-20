package utils

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

// SubMatrix returns a view corresponding to the submatrix of m starting at
// (row, col) of size (nr x nc). This uses mat.Dense.Slice.
func (m Matrix) SubMatrix(row, col, nr, nc int) Matrix {
	sub := m.M.Slice(row, row+nr, col, col+nc).(*mat.Dense)
	return Matrix{M: sub, DataP: sub.RawMatrix().Data}
}

// --------------------------------------------------------------------------
// SparseBlockMatrix definition and updated Mul method.
// --------------------------------------------------------------------------

// SparseBlockMatrix represents a sparse matrix whose nonzero entries are full blocks.
type SparseBlockMatrix struct {
	// Number of block rows and block columns.
	NrBlocks, NcBlocks int
	// Dimensions of each block. For simplicity we assume all blocks are of the same size.
	BlockRows, BlockCols int
	// CSR storage for the block pattern.
	RowPtr []int    // length = NrBlocks+1; for each block row, indices into ColIdx and Blocks
	ColIdx []int    // for each nonzero block, the block column index
	Blocks []Matrix // nonzero blocks stored in row-major order.
}

// Mul multiplies the sparse block matrix by a dense matrix x.
// x can have an arbitrary number of columns. In our setting, x is a "block vector"
// with dimension (NcBlocks*BlockCols) x p. The result will have dimensions (NrBlocks*BlockRows) x p.
func (sbm *SparseBlockMatrix) Mul(x Matrix) Matrix {
	_, xCols := x.Dims()
	totalRows := sbm.NrBlocks * sbm.BlockRows
	result := NewMatrix(totalRows, xCols)
	// For each block row...
	for i := 0; i < sbm.NrBlocks; i++ {
		// Get view for the i-th output block.
		yBlock := result.SubMatrix(i*sbm.BlockRows, 0, sbm.BlockRows, xCols)
		// Initialize yBlock to zero (NewMatrix zeros the data).
		// For each nonzero block in block row i:
		for k := sbm.RowPtr[i]; k < sbm.RowPtr[i+1]; k++ {
			j := sbm.ColIdx[k]
			// Get the corresponding block from x.
			xBlock := x.SubMatrix(j*sbm.BlockCols, 0, sbm.BlockCols, xCols)
			// Multiply the block: temp = A_ij * xBlock.
			temp := sbm.Blocks[k].Mul(xBlock)
			// Add the result into yBlock.
			yBlock = yBlock.Add(temp)
		}
		// (Since SubMatrix returns a view, changes are reflected in result.)
	}
	return result
}

// ------------------------------------------------------------
// GMRES (a very simplified sketch)
// ------------------------------------------------------------

// GMRES solves A x = b using the GMRES algorithm.
// x and b are assumed to be block vectors (i.e. column matrices whose height is a multiple of the block size).
// This is only a skeleton implementation.
func GMRES(A *SparseBlockMatrix, b Matrix, tol float64, maxIter int) Matrix {
	n, _ := b.Dims() // total number of rows in the full vector
	// Start with an initial guess x0 = 0.
	x := NewMatrix(n, 1)
	// Compute initial residual r0 = b - A*x. (x is zero so r0 = b.)
	r := b // in a full implementation, you’d copy b
	beta := mat.Norm(r.M, 2)
	if beta < tol {
		return x
	}

	// Allocate space for Arnoldi vectors. V[0] is r normalized.
	V := make([]Matrix, maxIter+1)
	V[0] = NewMatrix(n, 1)
	// Normalize r into V[0]:
	// (Here we create a temporary copy; a production version should avoid extra allocations.)
	for i := 0; i < n; i++ {
		V[0].M.Set(i, 0, r.M.At(i, 0)/beta)
	}

	// Hessenberg matrix H with dimensions (maxIter+1) x maxIter.
	H := mat.NewDense(maxIter+1, maxIter, nil)
	// g vector for the least squares problem.
	g := mat.NewVecDense(maxIter+1, nil)
	g.SetVec(0, beta)

	// GMRES iterations (very skeletal)
	var j int
	for j = 0; j < maxIter; j++ {
		// w = A * V[j]
		w := A.Mul(V[j])
		// Gram-Schmidt orthogonalization against V[0:j]
		for i := 0; i <= j; i++ {
			hij := mat.Dot(V[i].M.ColView(0), w.M.ColView(0))
			H.Set(i, j, hij)
			// Create a temporary matrix to hold hij * V[i]
			var scaled mat.Dense
			scaled.Scale(hij, V[i].M) // scaled = hij * V[i]
			// Subtract the scaled matrix from w: w = w - scaled
			w.M.Sub(w.M, &scaled)
		}
		// Compute h_{j+1,j} and so on...
	}

	// Solve the least squares problem: minimize || g - H y ||.
	// (Here you’d use a QR factorization of the (j+1)xj submatrix of H.)
	// For simplicity, we pretend y is computed and update x accordingly:
	// x = x + sum_{i=0}^{j-1} y_i * V[i]
	// (This step is left as an exercise to fill in the details.)

	// Return the computed x.
	return x
}

// GetBlock returns a Matrix handle for the block at the given blockRow and blockCol.
// If the block is not stored (i.e. is zero), it returns (Matrix{}, false).
func (sbm *SparseBlockMatrix) GetBlock(blockRow, blockCol int) (Matrix, bool) {
	if blockRow < 0 || blockRow >= sbm.NrBlocks {
		return Matrix{}, false
	}
	// For block row i, the nonzeros are stored from RowPtr[i] to RowPtr[i+1]-1.
	start := sbm.RowPtr[blockRow]
	end := sbm.RowPtr[blockRow+1]
	for i := start; i < end; i++ {
		if sbm.ColIdx[i] == blockCol {
			// Return the Matrix that wraps the block's allocated memory.
			return sbm.Blocks[i], true
		}
	}
	return Matrix{}, false
}

// SetBlock writes the provided matrix m into the block at (blockRow, blockCol).
// m must have the same dimensions as the blocks (BlockRows x BlockCols).
// It returns an error if the block location is not present or dimensions do not match.
func (sbm *SparseBlockMatrix) SetBlock(blockRow, blockCol int, m Matrix) error {
	if blockRow < 0 || blockRow >= sbm.NrBlocks {
		return fmt.Errorf("SetBlock: blockRow %d out of range", blockRow)
	}
	start := sbm.RowPtr[blockRow]
	end := sbm.RowPtr[blockRow+1]
	for i := start; i < end; i++ {
		if sbm.ColIdx[i] == blockCol {
			// Check dimensions.
			r, c := m.Dims()
			if r != sbm.BlockRows || c != sbm.BlockCols {
				return fmt.Errorf("SetBlock: dimensions mismatch; expected %dx%d, got %dx%d",
					sbm.BlockRows, sbm.BlockCols, r, c)
			}
			// Overwrite the existing block with m's data.
			// We assume that sbm.Blocks[i].M is already allocated to the right size.
			sbm.Blocks[i].M.Copy(m.M)
			return nil
		}
	}
	return fmt.Errorf("SetBlock: block at row %d, col %d not found", blockRow, blockCol)
}

// BlockPool allocates a contiguous pool of float64 values that will back a set of block matrices.
// Each block is of size (blockRows x blockCols) and the pool holds 'numBlocks' such blocks.
// The addresses slice indicates which [i,j] locations (block addresses) are allocated,
// and its length must equal numBlocks.
type BlockPool struct {
	data      []float64
	blockRows int
	blockCols int
	// addresses maps a block address (i,j) to the offset (in floats) of that block within data.
	addresses map[[2]int]int
}

// NewBlockPool allocates a new BlockPool. The addresses slice contains each allocated block’s address.
func NewBlockPool(numBlocks, blockRows, blockCols int, addresses [][2]int) *BlockPool {
	if len(addresses) != numBlocks {
		panic("number of addresses must equal numBlocks")
	}
	totalFloats := numBlocks * blockRows * blockCols
	data := make([]float64, totalFloats)
	addrMap := make(map[[2]int]int, numBlocks)
	// Each block gets a contiguous slice of length blockRows*blockCols.
	for i, addr := range addresses {
		offset := i * blockRows * blockCols
		addrMap[addr] = offset
	}
	return &BlockPool{
		data:      data,
		blockRows: blockRows,
		blockCols: blockCols,
		addresses: addrMap,
	}
}

// Block returns a Matrix view for the block at the given address (i, j).
func (bp *BlockPool) Block(i, j int) Matrix {
	offset, ok := bp.addresses[[2]int{i, j}]
	if !ok {
		panic(fmt.Sprintf("Block (%d,%d) not allocated", i, j))
	}
	// Create a subslice of the pool's data corresponding to this block.
	subData := bp.data[offset : offset+bp.blockRows*bp.blockCols]
	// Create a new Matrix that wraps subData.
	// This does not copy subData.
	m := NewMatrix(2, 2) // dummy allocation
	_ = m.ResetView(subData)
	return m
}
