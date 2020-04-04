package utils

import "gonum.org/v1/gonum/mat"

type Index []int

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
