package utils

import (
	"fmt"
	"math"
)

// BlockSparse represents a sparse block matrix. Only blocks provided via addresses are allocated;
// all other blocks are implicitly zero.
type BlockSparse struct {
	// Global block-matrix dimensions (in block counts).
	NrBlocks, NcBlocks int

	// Each block has dimensions blockRows x blockCols.
	blockRows, blockCols int

	// Contiguous storage for all allocated (nonzero) blocks.
	data []float64

	// addresses maps a block coordinate [i,j] to the offset (in floats) within data.
	addresses map[[2]int]int
}

// NewBlockSparse creates a new BlockSparse for a sparse block matrix.
// The input parameter addresses is a slice of [2]int specifying the coordinates
// of each nonzero block. The total number of allocated blocks is len(addresses).
// All other blocks (not in addresses) are implicitly zero.
func NewBlockSparse(nrBlocks, ncBlocks, blockRows, blockCols int, addresses [][2]int) *BlockSparse {
	totalBlocks := len(addresses)
	totalFloats := totalBlocks * blockRows * blockCols
	data := make([]float64, totalFloats)
	addrMap := make(map[[2]int]int, totalBlocks)
	// Each nonzero block gets a contiguous slice of length blockRows*blockCols.
	for i, addr := range addresses {
		offset := i * blockRows * blockCols
		addrMap[addr] = offset
	}
	return &BlockSparse{
		NrBlocks:  nrBlocks,
		NcBlocks:  ncBlocks,
		blockRows: blockRows,
		blockCols: blockCols,
		data:      data,
		addresses: addrMap,
	}
}

// GetBlockView returns a Matrix view for the block at coordinate (i, j).
// It uses your existing Matrix API (ResetView) so that the returned Matrix
// wraps the appropriate region of the BlockSparse’s contiguous data.
// If (i,j) is not allocated in this sparse matrix, the function panics.
func (bs *BlockSparse) GetBlockView(i, j int) Matrix {
	key := [2]int{i, j}
	offset, ok := bs.addresses[key]
	if !ok {
		panic(fmt.Sprintf("GetBlockView (%d,%d) not allocated", i, j))
	}
	// Create a view over bs.data for this block.
	subData := bs.data[offset : offset+bs.blockRows*bs.blockCols]
	// Create a dummy Matrix (which will be rebound via ResetView).
	m := NewMatrix(bs.blockRows, bs.blockCols)
	if err := m.ResetView(subData); err != nil {
		panic(err)
	}
	return m
}

// Mul multiplies bp (of dimensions m x n blocks) by other (of dimensions n x p blocks),
// using only the allocated blocks in each. That is, for each allocated block in bp
// at (i,k) and each allocated block in other at (k,j), the product is accumulated
// into the output block at (i,j).
//
// Since the sparsity pattern of the result is not known in advance, we accumulate results
// in a map and then build a new BlockSparse for the result. (For performance, one might eventually
// preallocate a workspace and reuse objects.)
func (bs *BlockSparse) Mul(other *BlockSparse) *BlockSparse {
	// Check block-matrix compatibility.
	if bs.NcBlocks != other.NrBlocks {
		panic("block matrix dimensions mismatch: bs.NcBlocks must equal other.NrBlocks")
	}
	// Let m = bs.NrBlocks, n = bs.NcBlocks (and other.NrBlocks), p = other.NcBlocks.
	// m, n, p := bs.NrBlocks, bs.NcBlocks, other.NcBlocks
	// Also check that the inner block dimensions are compatible.
	if bs.blockCols != other.blockRows {
		panic("block size mismatch: bs.blockCols must equal other.blockRows")
	}
	// The result's blocks have dimensions: blockRows (from bs) x other.blockCols.
	resBlockRows := bs.blockRows
	resBlockCols := other.blockCols

	// resultMap will accumulate computed blocks: key [i,j] -> yourmatrix.Matrix.
	resultMap := make(map[[2]int]Matrix)

	// Iterate over allocated blocks in bs.
	for key1 := range bs.addresses {
		// key1 = (i,k)
		i, k := key1[0], key1[1]
		// For each allocated block in other with row == k.
		for key2 := range other.addresses {
			if key2[0] != k {
				continue
			}
			j := key2[1]
			// Multiply bs.GetBlockView(i,k) by other.GetBlockView(k,j)
			A := bs.GetBlockView(i, k)
			B := other.GetBlockView(k, j)
			prod := A.Mul(B) // using your Matrix.Mul API that returns a new Matrix
			resKey := [2]int{i, j}
			if existing, ok := resultMap[resKey]; ok {
				// Accumulate: existing = existing + prod
				resultMap[resKey] = existing.Add(prod)
			} else {
				resultMap[resKey] = prod
			}
		}
	}

	// Now resultMap contains the nonzero result blocks.
	// Build the list of addresses for the result.
	var resAddresses [][2]int
	for key := range resultMap {
		resAddresses = append(resAddresses, key)
	}
	// Allocate a contiguous slice for the result blocks.
	numResBlocks := len(resAddresses)
	totalResFloats := numResBlocks * resBlockRows * resBlockCols
	resData := make([]float64, totalResFloats)
	// Build a new addresses map for the result.
	resAddrMap := make(map[[2]int]int, numResBlocks)
	// Now copy each result block's data into resData.
	// The order will be the order of resAddresses.
	for i, key := range resAddresses {
		offset := i * resBlockRows * resBlockCols
		resAddrMap[key] = offset
		// Get the block from resultMap.
		block := resultMap[key]
		// Copy its underlying data into the correct region of resData.
		copy(resData[offset:offset+resBlockRows*resBlockCols], block.DataP)
	}

	// Create and return the result BlockSparse.
	return &BlockSparse{
		NrBlocks:  bs.NrBlocks,    // result has the same block row dimension as bs
		NcBlocks:  other.NcBlocks, // result has the same block col dimension as other
		blockRows: resBlockRows,
		blockCols: resBlockCols,
		data:      resData,
		addresses: resAddrMap,
	}
}

// -----------------------------------------------------------------------------
// Helper routines for block vectors
// -----------------------------------------------------------------------------

// BlockFrobNorm computes the Frobenius norm of a block vector bp.
// It assumes bp is a BlockSparse with NcBlocks == 1.
func (bp *BlockSparse) FrobNorm() (norm float64) {
	sum := 0.0
	// Iterate over all allocated blocks in the dense vector.
	// For a dense block vector, every block row (i from 0 to NrBlocks-1) is allocated.
	for i := 0; i < bp.NrBlocks; i++ {
		// In a dense vector, the block address is (i,0)
		block := bp.GetBlockView(i, 0)
		// Compute the Frobenius norm of the 2x2 block.
		r, c := block.Dims()
		for i := 0; i < r; i++ {
			for j := 0; j < c; j++ {
				v := block.M.At(i, j)
				sum += v * v
			}
		}
	}
	norm = math.Sqrt(sum)
	return
}

// BlockInnerProduct computes the inner product between two block vectors x and y.
// The inner product is defined as the sum over blocks of the Frobenius inner product.
func (x *BlockSparse) InnerProduct(y *BlockSparse) float64 {
	if x.NrBlocks != y.NrBlocks {
		panic("InnerProduct: dimension mismatch")
	}
	sum := 0.0
	for i := 0; i < x.NrBlocks; i++ {
		xb := x.GetBlockView(i, 0)
		yb := y.GetBlockView(i, 0)
		r, c := xb.Dims()
		for i := 0; i < r; i++ {
			for j := 0; j < c; j++ {
				sum += xb.M.At(i, j) * yb.M.At(i, j)
			}
		}
	}
	return sum
}

// ScaleBlockPool scales every element of the block vector bp by alpha, in place.
func (bp *BlockSparse) Scale(alpha float64) {
	for i := 0; i < bp.NrBlocks; i++ {
		// For dense block vector, column is always 0.
		block := bp.GetBlockView(i, 0)
		r, c := block.Dims()
		for i := 0; i < r; i++ {
			for j := 0; j < c; j++ {
				v := block.M.At(i, j)
				block.M.Set(i, j, alpha*v)
			}
		}
	}
}

// Add returns a new block vector equal to x + y.
// It assumes x and y are both dense block vectors (NcBlocks==1) with identical dimensions.
func (x *BlockSparse) Add(y *BlockSparse) *BlockSparse {
	if x.NrBlocks != y.NrBlocks || x.blockRows != y.blockRows || x.blockCols != y.blockCols {
		panic("Add: dimension mismatch")
	}
	// Allocate a new dense block vector with the same dimensions.
	newAddrs := make([][2]int, x.NrBlocks)
	for i := 0; i < x.NrBlocks; i++ {
		newAddrs[i] = [2]int{i, 0}
	}
	res := NewBlockSparse(x.NrBlocks, 1, x.blockRows, x.blockCols, newAddrs)
	// For each block, add corresponding blocks.
	for i := 0; i < x.NrBlocks; i++ {
		xb := x.GetBlockView(i, 0)
		yb := y.GetBlockView(i, 0)
		r, c := xb.Dims()
		resBlock := res.GetBlockView(i, 0)
		for i := 0; i < r; i++ {
			for j := 0; j < c; j++ {
				sum := xb.M.At(i, j) + yb.M.At(i, j)
				resBlock.M.Set(i, j, sum)
			}
		}
	}
	return res
}

// SubtractBlockPool returns a new block vector equal to x - y.
func (x *BlockSparse) Subtract(y *BlockSparse) *BlockSparse {
	if x.NrBlocks != y.NrBlocks || x.blockRows != y.blockRows || x.blockCols != y.blockCols {
		panic("Subtract: dimension mismatch")
	}
	newAddrs := make([][2]int, x.NrBlocks)
	for i := 0; i < x.NrBlocks; i++ {
		newAddrs[i] = [2]int{i, 0}
	}
	res := NewBlockSparse(x.NrBlocks, 1, x.blockRows, x.blockCols, newAddrs)
	for i := 0; i < x.NrBlocks; i++ {
		xb := x.GetBlockView(i, 0)
		yb := y.GetBlockView(i, 0)
		r, c := xb.Dims()
		resBlock := res.GetBlockView(i, 0)
		for i := 0; i < r; i++ {
			for j := 0; j < c; j++ {
				diff := xb.M.At(i, j) - yb.M.At(i, j)
				resBlock.M.Set(i, j, diff)
			}
		}
	}
	return res
}

// CopyBlockPool makes a deep copy of a dense block vector.
func (bp *BlockSparse) Copy() *BlockSparse {
	newAddrs := make([][2]int, bp.NrBlocks)
	for i := 0; i < bp.NrBlocks; i++ {
		newAddrs[i] = [2]int{i, 0}
	}
	res := NewBlockSparse(bp.NrBlocks, 1, bp.blockRows, bp.blockCols, newAddrs)
	// Copy each block.
	for i := 0; i < bp.NrBlocks; i++ {
		src := bp.GetBlockView(i, 0)
		dst := res.GetBlockView(i, 0)
		r, c := src.Dims()
		for i := 0; i < r; i++ {
			for j := 0; j < c; j++ {
				dst.M.Set(i, j, src.M.At(i, j))
			}
		}
	}
	return res
}

// SolveLeastSquares is a stub for solving the small dense least-squares problem
// arising in GMRES (i.e. min||g - H y||). For our purposes, we assume H is a small
// (j+1)xj matrix stored as a slice of slices, and we return a vector y of length j.
// In a real implementation, you would call a LAPACK routine.
func SolveLeastSquares(H [][]float64, rows, cols int, g []float64) []float64 {
	// For simplicity, we assume j = cols and perform a dummy solve.
	// This is a placeholder.
	y := make([]float64, cols)
	// For testing, set y[i] = g[0] / float64(cols)
	for i := 0; i < cols; i++ {
		y[i] = g[0] / float64(cols)
	}
	return y
}

// -----------------------------------------------------------------------------
// GMRES Implementation using the BlockSparse API
// -----------------------------------------------------------------------------

// GMRES solves the linear system A x = b using the GMRES algorithm,
// where A is a block matrix (BlockSparse) and b is a dense block vector
// (a BlockSparse with NcBlocks == 1). The output x is a dense block vector.
func (A *BlockSparse) GMRES(b *BlockSparse, tol float64,
	maxIter int) *BlockSparse {
	// Assume x0 = 0.
	// For a dense block vector, we require that b.NcBlocks == 1.
	if b.NcBlocks != 1 {
		panic("GMRES requires b to be a dense block vector (NcBlocks == 1)")
	}
	m := A.NrBlocks // number of block rows in A
	// x0 = zero vector.
	x0 := NewBlockSparse(m, 1, b.blockRows, b.blockCols, generateDenseAddresses(m))
	// r0 = b - A * x0. Since x0 is zero, r0 = b.
	r0 := b
	beta := r0.FrobNorm()
	if beta < tol {
		return x0
	}
	// v0 = r0 / beta.
	v0 := r0.Copy()
	v0.Scale(1.0 / beta)
	V := make([]*BlockSparse, maxIter+1)
	V[0] = v0
	// H will be stored as a (maxIter+1) x maxIter dense matrix.
	H := make([][]float64, maxIter+1)
	for i := 0; i < maxIter+1; i++ {
		H[i] = make([]float64, maxIter)
	}
	g := make([]float64, maxIter+1)
	g[0] = beta

	var j int
	for j = 0; j < maxIter; j++ {
		// w = A * v_j. Here, v_j is a dense block vector.
		w := A.Mul(V[j])
		// Orthogonalize w against V[0] ... V[j].
		for i := 0; i <= j; i++ {
			hij := V[i].InnerProduct(w)
			H[i][j] = hij
			// w = w - hij * V[i].
			temp := V[i].Copy()
			temp.Scale(hij)
			w = w.Subtract(temp)
		}
		hj1j := w.FrobNorm()
		H[j+1][j] = hj1j
		if hj1j < tol {
			j++
			break
		}
		vNext := w.Copy()
		vNext.Scale(1.0 / hj1j)
		V[j+1] = vNext
	}
	// Solve the least-squares problem for y.
	y := SolveLeastSquares(H, j+1, j, g)
	// x = x0 + sum_{i=0}^{j-1} y[i] * V[i].
	x := x0.Copy()
	for i := 0; i < j; i++ {
		temp := V[i].Copy()
		temp.Scale(y[i])
		x = x.Add(temp)
	}
	return x
}

// generateDenseAddresses returns a slice of addresses for a dense block vector with numBlocks blocks.
// Each block is assumed to be at (i, 0) for i = 0..numBlocks-1.
func generateDenseAddresses(numBlocks int) [][2]int {
	addrs := make([][2]int, numBlocks)
	for i := 0; i < numBlocks; i++ {
		addrs[i] = [2]int{i, 0}
	}
	return addrs
}
