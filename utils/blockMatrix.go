package utils

import (
	"bytes"
	"fmt"
)

type BlockMatrix struct {
	M [][]Matrix // First slice points to rows of matrices - Note, the Matrix type allows for scalar matrices
	N int        // number of rows, columns, square block matrix only
	P []int      // Permutation "matrix", created during an LUP decomposition, otherwise nil
}

func NewBlockMatrix(N int) (R BlockMatrix) {
	R = BlockMatrix{
		N: N,
	}
	R.M = make([][]Matrix, N)
	for n := range R.M {
		R.M[n] = make([]Matrix, N)
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
