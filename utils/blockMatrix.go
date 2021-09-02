package utils

import (
	"bytes"
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

type BlockMatrix struct {
	M      [][]Matrix // First slice points to rows of matrices - Note, the Matrix type allows for scalar matrices
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
	R.M = make([][]Matrix, Nr)
	for n := range R.M {
		R.M[n] = make([]Matrix, Nc)
	}
	return R
}

func (bm BlockMatrix) Print() (out string) {
	var (
		output string
	)
	buf := bytes.Buffer{}
	for n, row := range bm.M {
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

func (bm BlockMatrix) IsSquare() bool {
	return bm.Nr == bm.Nc
}

func (bm *BlockMatrix) LUPDecompose() (err error) {
	/*
	   Factors the current matrix into a lower [L] and upper [U] pair of diagonal matrices such that [M] = [L]x[U]

	   Algorithm from: https://en.wikipedia.org/wiki/LU_decomposition#C_code_example

	   The matrix is factored in place, replacing the current matrices within by a new matrix composed of the
	   [L-E] and [U] matrices, stored in the same original matrix locations. The companion method to LUPD decompose is
	   LUPSolve(), which can be called repeatedly to efficiently produce solutions to the problem:
	                                       [M] * X = B
	   where [M] is this matrix, and B is the known RHS vector and X is the target.

	   Matrix M is changed, it contains a copy of both matrices L-I and U as (L-I)+U such that:
	                                       P * [M] = L * U
	*/

	var (
		imax       int
		absA, maxA float64
		Scratch    Matrix
	)
	if !bm.IsSquare() {
		err = fmt.Errorf("Matrix must be square")
		return
	}
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
			absA = math.Abs(mat.Det(bm.M[k][i]))
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
			// pivot rows of M
			bm.M[i], bm.M[imax] = bm.M[imax], bm.M[i]
			// counting pivots starting from N
			bm.Pcount++
		}
		for j := i + 1; j < bm.Nr; j++ {
			if Scratch, err = bm.M[i][i].Inverse(); err != nil {
				return
			}
			bm.M[j][i] = bm.M[j][i].Mul(Scratch)
			for k := i + 1; k < bm.Nr; k++ {
				Scratch = bm.M[j][i].Mul(bm.M[i][k])
				bm.M[j][k] = bm.M[j][k].Subtract(Scratch)
			}
		}
	}
	return
}

func (bm BlockMatrix) LUPSolve(b []Matrix) (Bx BlockMatrix, err error) {
	/*
	   Provided a solution vector B of size N x NB, calculate X for equation:
	       [M] * X = B
	   where [M] is the block matrix

	   Each sub-matrix within [M] is of size NBxNB
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
		[M] * X = B
		where [M] is the block matrix

		Each sub-matrix within [M] is of size NBxNB
		Each of the X and B vectors are of size NxNB
	*/
	// Allocate solution X
	Bx = NewBlockMatrix(bm.Nr, 1)
	for i := 0; i < bm.Nr; i++ {
		//		fmt.Printf("i = %d, P[i] = %d\n", i, bm.P[i])
		//		fmt.Printf("b[P[i]] = %s\n", b[bm.P[i]].Print())
		Bx.M[i][0] = b[bm.P[i]].Copy()
		for k := 0; k < i; k++ {
			Scratch = bm.M[i][k].Mul(Bx.M[k][0])
			Bx.M[i][0] = Bx.M[i][0].Subtract(Scratch)
		}
	}
	for i := bm.Nr - 1; i >= 0; i-- {
		for k := i + 1; k < bm.Nr; k++ {
			Scratch = bm.M[i][k].Mul(Bx.M[k][0])
			Bx.M[i][0] = Bx.M[i][0].Subtract(Scratch)
		}
		if Scratch, err = bm.M[i][i].Inverse(); err != nil {
			panic(err)
		}
		Bx.M[i][0] = Bx.M[i][0].Transpose().Mul(Scratch)
	}
	return
}

func (bm BlockMatrix) Mul(ba BlockMatrix) (R BlockMatrix) {
	var (
		Left, Right        = bm.M, ba.M
		NrLeft, NcLeft     = bm.Nr, bm.Nc
		NrRight, NcRight   = ba.Nr, ba.Nc
		NrTarget, NcTarget = NcRight, NrLeft
		Scratch            Matrix
	)
	/*
		fmt.Printf("Left dimensions: [%d,%d]\n", NrLeft, NcLeft)
		fmt.Printf("Right dimensions: [%d,%d]\n", NrRight, NcRight)
		fmt.Printf("Result dimensions: [%d,%d]\n", NcRight, NrLeft)
	*/
	if NrRight != NcLeft {
		panic(fmt.Errorf("number of rows in right Matrix should be %d, is %d", NcLeft, NrRight))
	}
	R = NewBlockMatrix(NrTarget, NcTarget)
	for j := 0; j < NcRight; j++ {
		for i := 0; i < NrLeft; i++ {
			// Iterate across columns of left and rows of right (NcLeft == NrRight) for sum at column j:0-NrLeft
			for ii := 0; ii < NcLeft; ii++ { // For each column in left, or row in right
				/*
					fmt.Printf("NrLeft,NcLeft,NrRight,NcRight - i,j,ii = [%d,%d],[%d,%d] %d,%d,%d\n",
						NrLeft, NcLeft, NrRight, NcRight, i, j, ii)
					fmt.Printf(Left[i][ii].Print("Left[i][ii]"))
					fmt.Printf(Right[ii][j].Print("Right[ii][j]"))
				*/
				Scratch = Left[i][ii].Mul(Right[ii][j])
				if ii == 0 {
					R.M[j][i] = Scratch
				} else {
					R.M[j][i] = R.M[j][i].Add(Scratch)
				}
				//fmt.Printf(R.M[i][j].Print("R[i][j]"))
				//os.Exit(1)
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
			bm.M[i][j].Add(ba.M[i][j])
		}
	}
	return
}

func (bm BlockMatrix) Copy() (R BlockMatrix) {
	var (
		Nr, Nc = bm.Nr, bm.Nc
	)
	R = NewBlockMatrix(Nr, Nc)
	for j := 0; j < Nc; j++ {
		for i := 0; i < Nr; i++ {
			if !bm.M[i][j].IsEmpty() {
				R.M[i][j] = bm.M[i][j].Copy()
			}
		}
	}
	return
}
