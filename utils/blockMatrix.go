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

func NewBlockMatrixFromScalar(A Matrix) (R BlockMatrix) {
	var (
		Nr, Nc = A.Dims()
	)
	R = BlockMatrix{
		Nr:  Nr,
		Nc:  Nc,
		tol: 0.00000001, // Default value
	}
	R.M = make([][]Matrix, Nr)
	for n := range R.M {
		R.M[n] = make([]Matrix, Nc)
	}
	for j := 0; j < Nr; j++ {
		for i := 0; i < Nc; i++ {
			R.M[i][j] = NewMatrix(1, 1, []float64{A.At(i, j)})
		}
	}
	return R
}

func (bm BlockMatrix) Print() (out string) {
	var (
		output string
		A      = bm.M
	)
	buf := bytes.Buffer{}
	for n, row := range A {
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
	                                       [M] * R = B
	   where [M] is this matrix, and B is the known RHS vector and R is the target.

	   Matrix M is changed, it contains a copy of both matrices L-I and U as (L-I)+U such that:
	                                       P * [M] = L * U
	*/

	var (
		imax       int
		absA, maxA float64
		Scratch    Matrix
		N          = bm.Nr
		A          = bm.M
	)
	if !bm.IsSquare() {
		err = fmt.Errorf("Matrix must be square")
		return
	}
	if len(bm.P) != 0 {
		err = fmt.Errorf("LUPDecompose already called on this matrix, which has overwritten it")
		return
	}
	bm.P = make([]int, N)
	for i := range bm.P {
		bm.P[i] = i
	}
	// counting pivots starting from N
	bm.Pcount = N // initialize Pcount with N
	for i := 0; i < N; i++ {
		maxA = 0.
		imax = i
		for k := 0; k < N; k++ {
			absA = math.Abs(mat.Det(A[k][i]))
			if absA > maxA {
				maxA = absA
				imax = k
			}
		}
		if maxA < bm.tol {
			err = fmt.Errorf("matrix is degenerate with tolerance %8.5e", bm.tol)
			return
		}
		if imax != i {
			// pivot P
			bm.P[i], bm.P[imax] = bm.P[imax], bm.P[i] // swap
			// pivot rows of M
			A[i], A[imax] = A[imax], A[i]
			// counting pivots starting from N
			bm.Pcount++
		}
		for j := i + 1; j < N; j++ {
			if Scratch, err = A[i][i].Inverse(); err != nil {
				return
			}
			A[j][i] = A[j][i].Mul(Scratch)
			for k := i + 1; k < N; k++ {
				A[j][k] = A[j][k].Subtract(A[j][i].Mul(A[i][k]))
			}
		}
	}
	return
}

func (bm BlockMatrix) LUPSolve(b []Matrix) (Bx BlockMatrix, err error) {
	/*
	   Provided a solution vector B of size N x NB, calculate R for equation:
	       [M] * R = B
	   where [M] is the block matrix

	   Each sub-matrix within [M] is of size NBxNB
	   Each of the R and B vectors are of size NxNB
	*/

	var (
		Scratch Matrix
		P       = bm.P
		N       = bm.Nr
		A       = bm.M
	)
	if len(P) == 0 {
		err = fmt.Errorf("uninitialized - call LUPDecompose first")
		return
	}
	/*
		Provided a solution vector B of size N x NB, calculate R for equation:
		[M] * R = B
		where [M] is the block matrix

		Each sub-matrix within [M] is of size NBxNB
		Each of the R and B vectors are of size NxNB
	*/
	// Allocate solution R
	Bx = NewBlockMatrix(N, 1)
	X := Bx.M
	for i := 0; i < N; i++ {
		X[i][0] = b[P[i]].Copy()
		for k := 0; k < i; k++ {
			X[i][0] = X[i][0].Subtract(A[i][k].Mul(X[k][0]))
		}
	}
	cDims := func(i, k int, a, x Matrix) {
		var (
			NrA, NcA = a.Dims()
			NrX, NcX = x.Dims()
		)
		fmt.Printf("[i,k] = [%d,%d], [NrA,NcA] = [%d,%d], [NrX,NcX] = [%d,%d]\n",
			i, k, NrA, NcA, NrX, NcX)
	}
	_ = cDims
	for i := N - 1; i >= 0; i-- {
		for k := i + 1; k < N; k++ {
			// cDims(i, k, A[i][k], R[k][0])
			// R[i][0] = R[i][0].Subtract(A[i][k].Mul(R[k][0]))
			X[i][0] = X[i][0].Subtract(A[i][k].Mul(X[k][0].Transpose()))
		}
		if Scratch, err = A[i][i].Inverse(); err != nil {
			panic(err)
		}
		X[i][0] = X[i][0].Transpose().Mul(Scratch)
	}
	for i := 0; i < N; i++ {
		X[i][0] = X[i][0].Transpose()
	}
	return
}

func (bm BlockMatrix) LUPInvert() (R BlockMatrix, err error) {
	var (
		N       = bm.Nr
		P       = bm.P
		A       = bm.M
		Scratch Matrix
	)
	if len(bm.P) == 0 {
		err = fmt.Errorf("uninitialized - call LUPDecompose first")
		return
	}
	zero := NewMatrix(1, 1, []float64{0.})
	one := NewMatrix(1, 1, []float64{1.})
	R = NewBlockMatrix(N, N)
	IA := R.M
	for j := 0; j < N; j++ {
		for i := 0; i < N; i++ {
			if P[i] == j {
				IA[i][j] = one.Copy()
			} else {
				IA[i][j] = zero.Copy()
			}
			for k := 0; k < i; k++ {
				IA[i][j] = IA[i][j].Subtract(A[i][k].Mul(IA[k][j]))
			}
		}
		for i := N - 1; i >= 0; i-- {
			for k := i + 1; k < N; k++ {
				IA[i][j] = IA[i][j].Subtract(A[i][k].Mul(IA[k][j]))
			}
			if Scratch, err = A[i][i].Inverse(); err != nil {
				panic(err)
			}
			IA[i][j] = IA[i][j].Mul(Scratch)
		}
	}
	return
}

func (bm BlockMatrix) LUPDeterminant() (det float64, err error) {
	var (
		N      = bm.Nr
		Pcount = bm.Pcount
		A      = bm.M
		P      = bm.P
	)
	if len(P) == 0 {
		err = fmt.Errorf("uninitialized - call LUPDecompose first")
		return
	}
	det = mat.Det(A[0][0])
	for i := 1; i < N; i++ {
		det *= mat.Det(A[i][i])
	}
	if (Pcount-N)%2 != 0 {
		det = -det
	}
	return
}

func (bm BlockMatrix) GetTol() (tol float64) {
	return bm.tol
}

func (bm BlockMatrix) Mul(ba BlockMatrix) (R BlockMatrix) {
	var (
		Left, Right        = bm.M, ba.M
		NrLeft, NcLeft     = bm.Nr, bm.Nc
		NrRight, NcRight   = ba.Nr, ba.Nc
		NrTarget, NcTarget = NcRight, NrLeft
		Scratch            Matrix
	)
	if NrRight != NcLeft {
		panic(fmt.Errorf("number of rows in right Matrix should be %d, is %d", NcLeft, NrRight))
	}
	R = NewBlockMatrix(NrTarget, NcTarget)
	R.tol = bm.tol
	for j := 0; j < NcRight; j++ {
		for i := 0; i < NrLeft; i++ {
			// Iterate across columns of left and rows of right (NcLeft == NrRight) for sum at column j:0-NrLeft
			for ii := 0; ii < NcLeft; ii++ { // For each column in left, or row in right
				if (Left[i][ii].IsEmpty() || Right[ii][j].IsEmpty()) ||
					(Left[i][ii].IsScalar() && Left[i][ii].DataP[0] == 0.) ||
					(Right[ii][j].IsScalar() && Right[ii][j].DataP[0] == 0.) {
					Scratch = NewMatrix(1, 1, []float64{0.})
				} else {
					Scratch = Left[i][ii].Mul(Right[ii][j])
				}
				if ii == 0 {
					R.M[j][i] = Scratch
				} else {
					R.M[j][i] = R.M[j][i].Add(Scratch)
				}
			}
		}
	}
	return
}

func (bm BlockMatrix) Add(ba BlockMatrix) {
	var (
		Nr, Nc = bm.Nr, bm.Nc
		A      = bm.M
	)
	for i := 0; i < Nr; i++ {
		for j := 0; j < Nc; j++ {
			A[i][j].Add(A[i][j])
		}
	}
	return
}

func (bm BlockMatrix) Copy() (R BlockMatrix) {
	var (
		Nr, Nc = bm.Nr, bm.Nc
		A      = bm.M
	)
	R = NewBlockMatrix(Nr, Nc)
	for j := 0; j < Nc; j++ {
		for i := 0; i < Nr; i++ {
			if !A[i][j].IsEmpty() {
				R.M[i][j] = A[i][j].Copy()
			}
		}
	}
	return
}

func (bm BlockMatrix) Transpose() (R BlockMatrix) {
	var (
		Nr, Nc = bm.Nr, bm.Nc
		A      = bm.M
	)
	R = NewBlockMatrix(Nc, Nr)
	for j := 0; j < Nc; j++ {
		for i := 0; i < Nr; i++ {
			if !A[i][j].IsEmpty() {
				R.M[j][i] = A[i][j].Copy()
			}
		}
	}
	return
}

func (bm BlockMatrix) Scale(val float64) (R BlockMatrix) {
	var (
		Nr, Nc = bm.Nr, bm.Nc
		A      = bm.M
	)
	for j := 0; j < Nc; j++ {
		for i := 0; i < Nr; i++ {
			A[i][j].Scale(val)
		}
	}
	return bm
}
