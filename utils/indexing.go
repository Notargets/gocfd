package utils

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type Index2D struct {
	RI, CI Index
	Len    int
}

func NewIndex2D(RI, CI Index) (I2 Index2D, err error) {
	if len(RI) != len(CI) {
		err = fmt.Errorf("lengths of row and column indices must be the same: nr, nc = %v, %v\n", len(RI), len(CI))
		return
	}
	return Index2D{
		RI:  RI,
		CI:  CI,
		Len: len(RI),
	}, nil
}

type Index []int

func NewIndex(N int) (I Index) {
	return make(Index, N)
}
func NewRangeOffset(rmin, rmax int) (r Index) {
	// Input range is "1 based" and converted to zero based index
	return NewRange(rmin-1, rmax-1)
}

func NewRange(rmin, rmax int) (r Index) {
	var (
		size = rmax - rmin + 1 // INCLUSIVE RANGE
	)
	r = make(Index, size)
	for i := range r {
		r[i] = i + rmin
	}
	return
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

func (I Index) Add(val int) (r Index) {
	r = make(Index, len(I))
	for i, ival := range I {
		r[i] = val + ival
	}
	return r
}
func (I Index) AddInPlace(val int) (r Index) {
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

func (I Index) Apply(f func(val int) int) (r Index) {
	r = make(Index, len(I))
	for i, val := range I {
		r[i] = f(val)
	}
	return
}

func (I Index) Outer(J Index) (A *mat.Dense) {
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
	A = mat.NewDense(ni, nj, data)
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

func (I Index) FindVec(op EvalOp, Values Index) (J Index) {
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
