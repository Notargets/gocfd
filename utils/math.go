package utils

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

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

func POW(x float64, pp int) (y float64) {
	var (
		p       = pp
		flipped bool
	)
	if p < 0 {
		p = -pp
		flipped = true
	}
	switch p {
	case 0:
		y = 1
	case 1:
		y = x
	case 2:
		y = x * x
	case 3:
		y = x * x * x
	case 4:
		y = x * x
		y = y * y
	case 5:
		y = x * x
		y = y * y * x
	case 6:
		y = x * x
		y = y * y * y
	case 7:
		y = x * x
		y = y * y * y * x
	case 8:
		y = x * x
		y = y * y * y * y
	default:
		y = math.Pow(x, float64(p))
	}
	if flipped {
		y = 1. / y
	}
	return
}
