package utils

import (
	"math"
)

func ConstArray(N int, val float64) (v []float64) {
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
	if pp > 8 || pp < -8 {
		goto MATHPOW
	}

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
	}
	if flipped {
		y = 1. / y
	}
	return

MATHPOW:
	y = math.Pow(x, float64(p))
	return
}
