package utils

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

func MatElementInvert(M mat.Matrix) (O *mat.Dense) {
	var (
		d, s = M.Dims()
	)
	O = mat.NewDense(d, s, nil)
	O.CloneFrom(M)
	for j := 0; j < s; j++ {
		for i := 0; i < d; i++ {
			O.Set(i, j, 1./O.At(i, j))
		}
	}
	return
}

func MatSubRow(MI mat.Matrix, RowIndices *mat.VecDense) (R *mat.Dense) {
	// RowIndices should contain a list of row indices into M
	var (
		nr, nc = MI.Dims()
		nrI    = RowIndices.Len()
		rI     = RowIndices.RawVector().Data
	)
	R = mat.NewDense(nrI, nc, nil)
	for i, val := range rI {
		valI := int(val)
		if valI > nr-1 || valI < 0 {
			fmt.Printf("index out of bounds: index = %d, max_bounds = %d\n", valI, nr-1)
			panic("unable to subset row from matrix")
		}
		var rowSlice []float64
		if M, ok := MI.(*mat.Dense); ok {
			rowSlice = M.RawRowView(valI)
		} else {
			rowSlice = make([]float64, nc)
			for j := 0; j < nc; j++ {
				rowSlice[j] = MI.At(i, j)
			}
		}
		R.SetRow(i, rowSlice)
	}
	return
}

func MatSubCol(MI mat.Matrix, ColIndices *mat.VecDense) (R *mat.Dense) {
	// ColIndices should contain a list of row indices into M
	var (
		nr, nc = MI.Dims()
		ncI    = ColIndices.Len()
		cI     = ColIndices.RawVector().Data
	)
	R = mat.NewDense(nr, ncI, nil)
	for j, val := range cI {
		valI := int(val)
		if valI > nc-1 || valI < 0 {
			panic("unable to subset row from matrix, index out of bounds")
		}
		var colSlice []float64
		if M, ok := MI.(*mat.Dense); ok {
			colSlice = VecGetF64(M.ColView(valI))
		} else {
			colSlice = make([]float64, nr)
			for i := 0; i < nr; i++ {
				colSlice[i] = MI.At(i, j)
			}
		}
		R.SetCol(j, colSlice)
	}
	return
}

func MatFind(MI mat.Matrix, val float64) (rowInd, colInd Index) {
	var (
		rows, cols = MI.Dims()
	)
	for j := 0; j < cols; j++ {
		for i := 0; i < rows; i++ {
			if MI.At(i, j) == val {
				rowInd = append(rowInd, i)
				colInd = append(colInd, j)
			}
		}
	}
	return
}
