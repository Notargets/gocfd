package utils

import "gonum.org/v1/gonum/mat"

func NewTriDiagonal(d0, d1, dm1 []float64) (Tri *mat.Dense) {
	dd := make([]float64, len(d0)*len(d0))
	var p1, p2, p3 int
	for j := 0; j < len(d0); j++ {
		for i := 0; i < len(d0); i++ {
			if i == j {
				dd[i+j*len(d0)] = d0[p1]
				p1++
				if i != len(d0)-1 {
					dd[+1+i+j*len(d0)] = d1[p2]
					p2++
				}
				if i != 0 {
					dd[-1+i+j*len(d0)] = dm1[p3]
					p3++
				}
			}
		}
	}
	Tri = mat.NewDense(len(d0), len(d0), dd)
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

func ConstArray(val float64, N int) (v []float64) {
	v = make([]float64, N)
	for i := range v {
		v[i] = val
	}
	return
}
