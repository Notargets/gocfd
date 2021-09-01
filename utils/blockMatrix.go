package utils

import (
	"bytes"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

type BlockMatrix struct {
	A      [][]Matrix // First slice points to rows of matrices - Note, the Matrix type allows for scalar matrices
	Nr, Nc int        // number of rows, columns in the square block matrix consisting of a sub-matrix in each cell
	P      []int      // Permutation "matrix", created during an LUP decomposition, otherwise nil
	Pcount int        // count of number of pivots, used in determining sign of determinant
	tol    float64    // tolerance for reduction operations, like LUP decomposition
}

func NewBlockMatrix(Nr, Nc int) (R BlockMatrix) {
	R = BlockMatrix{
		Nr:  Nr,
		Nc:  Nc,
		tol: 0.00000001, // Default value
	}
	R.A = make([][]Matrix, Nr)
	for n := range R.A {
		R.A[n] = make([]Matrix, Nc)
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
	bm.P = make([]int, bm.Nr)
	for i := range bm.P {
		bm.P[i] = i
	}
	// counting pivots starting from N
	bm.Pcount = bm.Nr // initialize Pcount with N
	for i := 0; i < bm.Nr; i++ {
		maxA = 0.
		imax = i
		for k := 0; k < bm.Nr; k++ {
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
		for j := i + 1; j < bm.Nr; j++ {
			if Scratch, err = bm.A[i][i].Inverse(); err != nil {
				return
			}
			bm.A[j][i] = bm.A[j][i].Mul(Scratch)
			for k := i + 1; k < bm.Nr; k++ {
				Scratch = bm.A[j][i].Mul(bm.A[i][k])
				bm.A[j][k] = bm.A[j][k].Subtract(Scratch)
			}
		}
	}
	return
}

func (bm BlockMatrix) LUPSolve(b []Matrix) (Bx BlockMatrix, err error) {
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
	// Allocate solution X
	Bx = NewBlockMatrix(bm.Nr, 1)
	for i := 0; i < bm.Nr; i++ {
		//		fmt.Printf("i = %d, P[i] = %d\n", i, bm.P[i])
		//		fmt.Printf("b[P[i]] = %s\n", b[bm.P[i]].Print())
		Bx.A[i][0] = b[bm.P[i]].Copy()
		for k := 0; k < i; k++ {
			Scratch = bm.A[i][k].Mul(Bx.A[k][0])
			Bx.A[i][0] = Bx.A[i][0].Subtract(Scratch)
		}
	}
	for i := bm.Nr - 1; i >= 0; i-- {
		for k := i + 1; k < bm.Nr; k++ {
			Scratch = bm.A[i][k].Mul(Bx.A[k][0])
			Bx.A[i][0] = Bx.A[i][0].Subtract(Scratch)
		}
		if Scratch, err = bm.A[i][i].Inverse(); err != nil {
			panic(err)
		}
		Bx.A[i][0] = Bx.A[i][0].Transpose().Mul(Scratch)
	}
	return
}

func (bm BlockMatrix) Mul(ba BlockMatrix) (R BlockMatrix) {
	var (
		Nr, Nc   = bm.Nr, bm.Nc
		Nra, Nca = ba.Nr, ba.Nc
	)
	R = NewBlockMatrix(Nca, Nr)
	for i := 0; i < Nca; i++ {
		for j := 0; j < Nc; j++ {
			for ii := 0; ii < Nra; ii++ {
				if ii == 0 {
					R.A[i][j] = bm.A[ii][j].Mul(ba.A[j][ii])
				} else {
					R.A[i][j] = R.A[i][j].Add(bm.A[ii][j].Mul(ba.A[j][ii]))
				}
			}
		}
	}
	return
}

func (bm BlockMatrix) Add(ba BlockMatrix) {
	var (
		Nr, Nc = bm.Nr, bm.Nc
	)
	for i := 0; i < Nr; i++ {
		for j := 0; j < Nc; j++ {
			bm.A[i][j].Add(ba.A[i][j])
		}
	}
	return
}
