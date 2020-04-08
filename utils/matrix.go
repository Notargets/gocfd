package utils

import "gonum.org/v1/gonum/mat"

func MatFind(MI mat.Matrix, op EvalOp, val float64) (I Index2D) {
	var (
		nr, nc         = MI.Dims()
		rowInd, colInd Index
	)
	switch op {
	case Equal:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				if MI.At(i, j) == val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case Less:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				if MI.At(i, j) < val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case LessOrEqual:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				if MI.At(i, j) <= val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case Greater:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				if MI.At(i, j) > val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case GreaterOrEqual:
		for j := 0; j < nc; j++ {
			for i := 0; i < nr; i++ {
				if MI.At(i, j) >= val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	}
	I, _ = NewIndex2D(nr, nc, rowInd, colInd)
	return
}

func NewSymTriDiagonal(d0, d1 []float64) (Tri *mat.SymDense) {
	dd := make([]float64, len(d0)*len(d0))
	var p1, p2 int
	for j := 0; j < len(d0); j++ {
		for i := 0; i < len(d0); i++ {
			if i == j {
				dd[i+j*len(d0)] = d0[p1]
				p1++
				if i != len(d0)-1 {
					dd[+1+i+j*len(d0)] = d1[p2]
					p2++
				}
			}
		}
	}
	Tri = mat.NewSymDense(len(d0), dd)
	return
}
