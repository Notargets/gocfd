package utils

import (
	"bytes"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

type BlockMatrix struct {
	A      [][]Matrix // First slice points to rows of matrices - Note, the Matrix type allows for scalar matrices
	N      int        // number of rows, columns in the square block matrix consisting of a sub-matrix in each cell
	NB     int        // number of rows, columns in each square sub-matrix
	P      []int      // Permutation "matrix", created during an LUP decomposition, otherwise nil
	Pcount int        // count of number of pivots, used in determining sign of determinant
	tol    float64    // tolerance for reduction operations, like LUP decomposition
}

func NewBlockMatrix(N, NB int) (R BlockMatrix) {
	R = BlockMatrix{
		N:   N,
		NB:  NB,
		tol: 0.00000001, // Default value
	}
	R.A = make([][]Matrix, N)
	for n := range R.A {
		R.A[n] = make([]Matrix, N)
	}
	return R
}

func (bm BlockMatrix) Print() (out string) {
	var (
		output string
	)
	buf := bytes.Buffer{}
	for n, row := range bm.A {
		for m, Mat := range row {
			label := fmt.Sprintf("[%d:%d]", n, m)
			if Mat.IsEmpty() {
				output = label + " nil "
			} else {
				output = Mat.Print(label)
			}
			buf.WriteString(fmt.Sprintf("%s", output))
		}
		buf.WriteString("\n")
	}
	return buf.String()
}

func (bm *BlockMatrix) LUPDecompose() (err error) {
	/*
	   Factors the current matrix into a lower [L] and upper [U] pair of diagonal matrices such that [A] = [L]x[U]

	   Algorithm from: https://en.wikipedia.org/wiki/LU_decomposition#C_code_example

	   The matrix is factored in place, replacing the current matrices within by a new matrix composed of the
	   [L-E] and [U] matrices, stored in the same original matrix locations. The companion method to LUPD decompose is
	   LUPSolve(), which can be called repeatedly to efficiently produce solutions to the problem:
	                                       [A] * X = B
	   where [A] is this matrix, and B is the known RHS vector and X is the target.

	   Matrix A is changed, it contains a copy of both matrices L-I and U as (L-I)+U such that:
	                                       P * [A] = L * U
	*/

	var (
		imax       int
		absA, maxA float64
		Scratch    Matrix
	)
	if len(bm.P) != 0 {
		err = fmt.Errorf("LUPDecompose already called on this matrix, which has overwritten it")
		return
	}
	bm.P = make([]int, bm.N)
	for i := range bm.P {
		bm.P[i] = i
	}
	// counting pivots starting from N
	bm.Pcount = bm.N // initialize Pcount with N
	for i := 0; i < bm.N; i++ {
		maxA = 0.
		imax = i
		for k := 0; k < bm.N; k++ {
			absA = math.Abs(mat.Det(bm.A[k][i]))
			if absA > maxA {
				maxA = absA
				imax = k
			}
		}
		if maxA < bm.tol {
			err = fmt.Errorf("matrix is degenerate with tolerance %8.5f", bm.tol)
			return
		}
		if imax != i {
			// pivot P
			bm.P[i], bm.P[imax] = bm.P[imax], bm.P[i] // swap
			// pivot rows of A
			bm.A[i], bm.A[imax] = bm.A[imax], bm.A[i]
			// counting pivots starting from N
			bm.Pcount++
		}
		for j := i + 1; j < bm.N; j++ {
			if Scratch, err = bm.A[i][i].Inverse(); err != nil {
				return
			}
			bm.A[j][i] = bm.A[j][i].Mul(Scratch)
			for k := i + 1; k < bm.N; k++ {
				Scratch = bm.A[j][i].Mul(bm.A[i][k])
				bm.A[j][k] = bm.A[j][k].Subtract(Scratch)
			}
		}
	}
	return
}

func (bm BlockMatrix) LUPSolve(b []Matrix) (x []Matrix, err error) {
	/*
	   Provided a solution vector B of size N x NB, calculate X for equation:
	       [A] * X = B
	   where [A] is the block matrix

	   Each sub-matrix within [A] is of size NBxNB
	   Each of the X and B vectors are of size NxNB
	*/

	var (
		Scratch Matrix
	)
	if len(bm.P) == 0 {
		err = fmt.Errorf("uninitialized - call LUPDecompose first")
		return
	}
	/*
		Provided a solution vector B of size N x NB, calculate X for equation:
		[A] * X = B
		where [A] is the block matrix

		Each sub-matrix within [A] is of size NBxNB
		Each of the X and B vectors are of size NxNB
	*/
	{ // validate the size of the input vector B
		var fail bool
		if len(b) != bm.N {
			fail = true
		}
		var nr, nc int
		for i := 0; i < bm.N; i++ {
			nr, nc = b[i].Dims()
			if nr != bm.NB || nc != bm.NB {
				fail = true
			}
		}
		if fail {
			err = fmt.Errorf("solution vector must have size %dx%d, provided size is %dx%d",
				bm.NB, bm.NB, nr, nc)
			return
		}
	}
	// Allocate solution X
	x = make([]Matrix, bm.N)
	for i := 0; i < bm.N; i++ {
		//		fmt.Printf("i = %d, P[i] = %d\n", i, bm.P[i])
		//		fmt.Printf("b[P[i]] = %s\n", b[bm.P[i]].Print())
		x[i] = b[bm.P[i]].Copy()
		for k := 0; k < i; k++ {
			Scratch = bm.A[i][k].Mul(x[k])
			x[i] = x[i].Subtract(Scratch)
		}
	}
	for i := bm.N - 1; i >= 0; i-- {
		for k := i + 1; k < bm.N; k++ {
			Scratch = bm.A[i][k].Mul(x[k])
			x[i] = x[i].Subtract(Scratch)
		}
		if Scratch, err = bm.A[i][i].Inverse(); err != nil {
			panic(err)
		}
		x[i] = x[i].Transpose().Mul(Scratch)
	}
	return
}

func (bm BlockMatrix) Mul(ba BlockMatrix) (R BlockMatrix) {
	var (
		err   error
		N, NB = bm.N, bm.NB
	)
	// Validate input
	if N != ba.N || NB != ba.NB {
		err = fmt.Errorf("dimensions of input matrix (%d,%d) do not match target (%d,%d)",
			N, NB, ba.N, ba.NB)
		panic(err)
	}
	R = NewBlockMatrix(N, NB)
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			R.A[i][j] = NewMatrix(NB, NB)
			for ii := 0; ii < N; ii++ {
				R.A[i][j] = R.A[i][j].Add(bm.A[ii][j].Mul(ba.A[j][ii]))
			}
		}
	}
	return
}
func (bm BlockMatrix) Add(ba BlockMatrix) {
	var (
		err   error
		N, NB = bm.N, bm.NB
	)
	// Validate input
	if N != ba.N || NB != ba.NB {
		err = fmt.Errorf("dimensions of input matrix (%d,%d) do not match target (%d,%d)",
			N, NB, ba.N, ba.NB)
		panic(err)
	}
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			bm.A[i][j].Add(ba.A[i][j])
		}
	}
	return
}
