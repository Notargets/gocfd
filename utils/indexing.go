package utils

import (
	"fmt"
)

type Index2D struct {
	// Ind is the combined index in column-major format
	RI, CI, Ind Index
	Len         int
}

func NewIndex2D(nr, nc int, RI, CI Index, permuteO ...bool) (I2 Index2D, err error) {
	var (
		permute bool
	)
	if len(permuteO) != 0 {
		permute = permuteO[0]
	} else {
		if len(RI) != len(CI) {
			err = fmt.Errorf("lengths of row and column indices must be the same: nr, nc = %v, %v\n", len(RI), len(CI))
			return
		}
	}
	I2 = Index2D{
		RI: RI,
		CI: CI,
	}
	switch permute {
	case false:
		I2.Ind = make(Index, len(RI))
		for i, ci := range CI {
			ri := RI[i]
			if ci > nc-1 || ci < 0 {
				err = fmt.Errorf("index dimension exceeds bounds: ci=%v, ciMax=%v\n", ci, nc-1)
				return
			}
			if ri > nr-1 || ri < 0 {
				err = fmt.Errorf("index dimension exceeds bounds: ri=%v, riMax=%v\n", ri, nr-1)
				return
			}
			I2.Ind[i] = ci + nc*ri
		}
	case true:
		I2.Ind = make(Index, len(RI)*len(CI))
		var ind int
		for _, ri := range RI {
			for _, ci := range CI {
				if ci > nc-1 || ci < 0 {
					err = fmt.Errorf("index dimension exceeds bounds: ci=%v, ciMax=%v\n", ci, nc-1)
					return
				}
				if ri > nr-1 || ri < 0 {
					err = fmt.Errorf("index dimension exceeds bounds: ri=%v, riMax=%v\n", ri, nr-1)
					return
				}
				I2.Ind[ind] = ci + nc*ri
				ind++
			}
		}
	}
	I2.Len = len(I2.Ind)
	return
}

func (i2 *Index2D) ToIndex() (I Index) {
	I = i2.Ind
	return
}

type Index []int

// Chainable methods
func NewIndex(N int, ValIO ...interface{}) (I Index) {
	I = make(Index, N)
	if len(ValIO) != 0 {
		ValI := ValIO[0]
		switch Val := ValI.(type) {
		case []float64:
			for i := range Val {
				I[i] = int(Val[i])
			}
		case []int:
			for i := range Val {
				I[i] = Val[i]
			}
		case Index:
			for i := range Val {
				I[i] = Val[i]
			}
		}
	}
	return
}
func NewRangeOffset(rmin, rmax int) (r Index) {
	// Input range is "1 based" and converted to zero based index
	return NewRangeInclusive(rmin-1, rmax-1)
}
func NewRangeInclusive(rmin, rmax int) (r Index) {
	var (
		size = rmax - rmin + 1 // INCLUSIVE RANGE
	)
	r = make(Index, size)
	for i := range r {
		r[i] = i + rmin
	}
	return
}

func NewRange(rmin, rmax int) (r Index) {
	return NewRangeInclusive(rmin, rmax-1)
}
func NewOnes(N int) (r Index) {
	r = make(Index, N)
	for i := 0; i < N; i++ {
		r[i] = 1
	}
	return
}
func NewFromFloat(IF []float64) (r Index) {
	r = make(Index, len(IF))
	for i, val := range IF {
		r[i] = int(val)
	}
	return
}

func (I Index) Copy() (r Index) {
	r = make(Index, len(I))
	for i, val := range I {
		r[i] = val
	}
	return r
}

func (I Index) Add(val int) Index {
	for i := range I {
		I[i] += val
	}
	return I
}

func (I Index) Subset(J Index) (r Index) {
	r = make(Index, len(J))
	for j, val := range J {
		r[j] = I[val]
	}
	return
}

func (I Index) Apply(f func(val int) int) Index {
	for i, val := range I {
		I[i] = f(val)
	}
	return I
}

func (I Index) Compare(op EvalOp, Values Index) (J Index) {
	/*
		Each element of Values is compared to the corresponding value of I:
		if (Values[i] op I[i]): append i to the output index J
	*/
	switch op {
	case Equal:
		for i, val := range I {
			if val == Values[i] {
				J = append(J, i)
			}
		}
	case Less:
		for i, val := range I {
			if val < Values[i] {
				J = append(J, i)
			}
		}
	case LessOrEqual:
		for i, val := range I {
			if val <= Values[i] {
				J = append(J, i)
			}
		}
	case Greater:
		for i, val := range I {
			if val > Values[i] {
				J = append(J, i)
			}
		}
	case GreaterOrEqual:
		for i, val := range I {
			if val >= Values[i] {
				J = append(J, i)
			}
		}
	}
	return
}

func (I Index) RowMajorToColumnMajor(nr, nc int) Index {
	for i, val := range I {
		I[i] = RowMajorToColMajor(nr, nc, val)
	}
	return I
}

// Non chainable methods
func (I Index) Outer(J Index) (A Matrix) {
	var (
		ni   = len(I)
		nj   = len(J)
		data []float64
		ind  int
	)
	data = make([]float64, ni*nj)
	for i := 0; i < ni; i++ {
		for j := 0; j < nj; j++ {
			data[ind] = float64(J[j] * I[i])
			ind++
		}
	}
	A = NewMatrix(ni, nj, data)
	return
}

func (I *Index) IndexedAssign(J, Val Index) (err error) {
	var (
		nr = len(*I)
		N  = len(J)
	)
	if nr == 0 {
		var jmax int
		for _, val := range J {
			if val > jmax {
				jmax = val
			}
		}
		jmax += 1 // indices are zero based, dimensions not
		*I = make(Index, jmax)
		nr = jmax
		fmt.Println("New jmax = ", jmax)
	}
	switch {
	case N != len(Val):
		err = fmt.Errorf("dimension mismatch: index and values should have the same length")
		return
	}
	for i, val := range Val {
		ji := J[i]
		switch {
		case ji < 0:
			err = fmt.Errorf("dimension bounds error, row index < 0: ji = %v\n", ji)
			return
		case ji > nr-1:
			err = fmt.Errorf("dimension bounds error, row index > max: ji = %v, max = %v\n", ji, nr-1)
			return
		}
		(*I)[ji] = val
	}
	return
}

func (I Index) ToMatrix(nr, nc int) (R Matrix) {
	var (
		dataR = make([]float64, nr*nc)
	)
	if nr*nc != len(I) {
		err := fmt.Errorf("dimensions do not match nr, nc = %v, %v, len(I) = %v\n", nr, nc, len(I))
		panic(err)
	}
	for i, val := range I {
		dataR[i] = float64(val)
	}
	R = NewMatrix(nr, nc, dataR)
	return
}

func (I Index) ToMatrixReversed(nr, nc int) (R Matrix) {
	var (
		dataR = make([]float64, nr*nc)
	)
	if nr*nc != len(I) {
		err := fmt.Errorf("dimensions do not match nr, nc = %v, %v, len(I) = %v\n", nr, nc, len(I))
		panic(err)
	}
	for i, val := range I {
		ind := RowMajorToColMajor(nr, nc, i)
		dataR[ind] = float64(val)
	}
	R = NewMatrix(nr, nc, dataR)
	return
}

func ColMajorIndex(nc, i, j int) (ind int) {
	return j + nc*i
}
