package utils

import (
	"fmt"
	"math"

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

func MatCopyEmpty(M *mat.Dense) (R *mat.Dense) {
	var (
		nr, nc = M.Dims()
	)
	R = mat.NewDense(nr, nc, nil)
	return
}

func MatSubset(M *mat.Dense, I Index) (r *mat.VecDense) {
	/*
		Index should contain a list of indices into MI
		Note: native mat library matrix storage is in column traversal first (row-major) order
	*/
	var (
		Mr     = M.RawMatrix()
		nr, nc = M.Dims()
		data   = make([]float64, len(I))
	)
	for i, ind := range I {
		data[i] = Mr.Data[RowMajorToColMajor(nr, nc, ind)]
	}
	r = mat.NewVecDense(len(I), data)
	return
}

func RowMajorToColMajor(nr, nc, ind int) (cind int) {
	// ind = i + nr * j
	// ind / nr = 0 + j
	j := ind / nr
	i := ind - nr*j
	cind = j + nc*i
	return
}

func MatSubsetRow(MI mat.Matrix, RowIndices *mat.VecDense) (R *mat.Dense) {
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

func MatSubsetCol(MI mat.Matrix, ColIndices *mat.VecDense) (R *mat.Dense) {
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

func MatFind(MI mat.Matrix, op EvalOp, val float64) (I Index2D) {
	var (
		rows, cols     = MI.Dims()
		rowInd, colInd Index
	)
	switch op {
	case Equal:
		for j := 0; j < cols; j++ {
			for i := 0; i < rows; i++ {
				if MI.At(i, j) == val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case Less:
		for j := 0; j < cols; j++ {
			for i := 0; i < rows; i++ {
				if MI.At(i, j) < val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case LessOrEqual:
		for j := 0; j < cols; j++ {
			for i := 0; i < rows; i++ {
				if MI.At(i, j) <= val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case Greater:
		for j := 0; j < cols; j++ {
			for i := 0; i < rows; i++ {
				if MI.At(i, j) > val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	case GreaterOrEqual:
		for j := 0; j < cols; j++ {
			for i := 0; i < rows; i++ {
				if MI.At(i, j) >= val {
					rowInd = append(rowInd, i)
					colInd = append(colInd, j)
				}
			}
		}
	}
	I, _ = NewIndex2D(rowInd, colInd)
	return
}

func MatIndexedAssign(MI *mat.Dense, I2 Index2D, Val Index) (err error) {
	// RI and CI are the row and column indices
	var (
		nr, nc = MI.Dims()
		RI, CI = I2.RI, I2.CI
		N      = len(RI)
	)
	switch {
	case N != len(CI):
		err = fmt.Errorf("dimension mismatch: RI and CI should have the same length")
		return
	case N != len(Val):
		err = fmt.Errorf("dimension mismatch: Val and RI,CI should have the same length")
		return
	}
	for i, val := range Val {
		ri, ci := RI[i], CI[i]
		switch {
		case ri < 0:
			err = fmt.Errorf("dimension bounds error, row index < 0: ri = %v\n", ri)
			return
		case ci < 0:
			err = fmt.Errorf("dimension bounds error, col index < 0: ci = %v\n", ci)
			return
		case ri > nr-1:
			err = fmt.Errorf("dimension bounds error, row index > max: ri = %v\n", ri)
			return
		case ci > nc-1:
			err = fmt.Errorf("dimension bounds error, col index > max: ci = %v\n", ci)
			return
		}
		MI.Set(RI[i], CI[i], float64(val))
	}
	return
}

func MatInPlace(M *mat.Dense, op MathOp, args ...float64) {
	var (
		data = M.RawMatrix().Data
	)
	switch op {
	case Abs:
		for i, val := range data {
			data[i] = math.Abs(val)
		}
	case Sqrt:
		for i, val := range data {
			data[i] = math.Sqrt(val)
		}
	case Pow:
		for i, val := range data {
			p := int(args[0])
			data[i] = POW(val, p)
		}
	}
}
