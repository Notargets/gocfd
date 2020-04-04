package utils

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type Index []int

func NewRangeOffset(rmin, rmax int) (r Index) {
	// Input range is "1 based" and converted to zero based index
	return NewRange(rmin-1, rmax-1)
}

func NewRange(rmin, rmax int) (r Index) {
	var (
		size = rmax - rmin + 1
	)
	r = make([]int, size)
	var ind int
	for i := rmin; i <= rmax; i++ {
		r[ind] = i
		ind++
	}
	return
}
func NewOnes(N int) (r Index) {
	r = make([]int, N)
	for i := 0; i < N; i++ {
		r[i] = 1
	}
	return
}

func (I Index) ApplyFunc(f func(val int) int) (r Index) {
	r = make([]int, len(I))
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

func (I Index) IndexedAssign(J, Val Index) (err error) {
	var (
		nr = len(I)
		N  = len(J)
	)
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
			err = fmt.Errorf("dimension bounds error, row index > max: ji = %v\n", ji)
			return
		}
		I[ji] = val
	}
	return
}
