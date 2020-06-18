package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/DG1D"

	"github.com/stretchr/testify/assert"
)

func TestElements2D(t *testing.T) {
	{
		N := 1
		x, y := Nodes2D(N)
		assert.True(t, nearVec([]float64{-1, 1, 0}, x.Data(), 0.0001))
		assert.True(t, nearVec([]float64{-0.5774, -0.5774, 1.1547}, y.Data(), 0.0001))
		N = 2
		x, y = Nodes2D(N)
		assert.True(t, nearVec([]float64{-1, 0, 1, -.5, .5, 0}, x.Data(), 0.0001))
		assert.True(t, nearVec([]float64{-0.5774, -0.5774, -0.5774, 0.2887, 0.2887, 1.1547}, y.Data(), 0.0001))
		r, s := XYtoRS(Nodes2D(N))
		assert.True(t, nearVec([]float64{-1, 0, 1, -1, 0, -1}, r.Data(), 0.0001))
		assert.True(t, nearVec([]float64{-1, -1, -1, 0, 0, 1}, s.Data(), 0.0001))
		a, b := RStoAB(r, s)
		assert.True(t, nearVec([]float64{-1, 0, 1, -1, 1, -1}, a.Data(), 0.0001))
		assert.True(t, nearVec([]float64{-1, -1, -1, 0, 0, 1}, b.Data(), 0.0001))

		h1 := DG1D.JacobiP(a, 0, 0, 0)
		assert.True(t, nearVec([]float64{0.7071, 0.7071, 0.7071, 0.7071, 0.7071, 0.7071}, h1, 0.0001))
		h2 := DG1D.JacobiP(b, float64(1), 0, 0)
		assert.True(t, nearVec([]float64{0.7071, 0.7071, 0.7071, 0.7071, 0.7071, 0.7071}, h2, 0.0001))
		P := Simplex2DP(a, b, 0, 0)
		assert.True(t, nearVec([]float64{0.7071, 0.7071, 0.7071, 0.7071, 0.7071, 0.7071}, P, 0.0001))

		V := Vandermonde2D(N, r, s)
		assert.True(t, nearVec([]float64{
			0.7071, -1.0000, 1.2247, -1.7321, 2.1213, 2.7386,
			0.7071, -1.0000, 1.2247, 0, 0, -1.3693,
			0.7071, -1.0000, 1.2247, 1.7321, -2.1213, 2.7386,
			0.7071, 0.5000, -0.6124, -0.8660, -1.5910, 0.6847,
			0.7071, 0.5000, -0.6124, 0.8660, 1.5910, 0.6847,
			0.7071, 2.0000, 3.6742, 0, 0, 0,
		}, V.Data(), 0.0001))
		V2Dr, V2Ds := GradVandermonde2D(N, r, s)
		assert.True(t, nearVec([]float64{
			0, 0, 0, 1.7321, -2.1213, -8.2158,
			0, 0, 0, 1.7321, -2.1213, 0,
			0, 0, 0, 1.7321, -2.1213, 8.2158,
			0, 0, 0, 1.7321, 3.1820, -4.1079,
			0, 0, 0, 1.7321, 3.1820, 4.1079,
			0, 0, 0, 1.7321, 8.4853, 0,
		}, V2Dr.Data(), 0.0001))
		assert.True(t, nearVec([]float64{
			0, 1.5000, -4.8990, 0.8660, -6.3640, -2.7386,
			0, 1.5000, -4.8990, 0.8660, -1.0607, 1.3693,
			0, 1.5000, -4.8990, 0.8660, 4.2426, 5.4772,
			0, 1.5000, 1.2247, 0.8660, -1.0607, -1.3693,
			0, 1.5000, 1.2247, 0.8660, 4.2426, 2.7386,
			0, 1.5000, 7.3485, 0.8660, 4.2426, 0,
		}, V2Ds.Data(), 0.0001))
	}
}

func nearVec(a, b []float64, tol float64) (l bool) {
	for i, val := range a {
		if !near(b[i], val, tol) {
			fmt.Printf("Diff = %v, Left[%d] = %v, Right[%d] = %v\n", math.Abs(val-b[i]), i, val, i, b[i])
			return false
		}
	}
	return true
}

func near(a, b float64, tolI ...float64) (l bool) {
	var (
		tol float64
	)
	if len(tolI) == 0 {
		tol = 1.e-08
	} else {
		tol = tolI[0]
	}
	if math.Abs(a-b) <= tol*math.Abs(a) {
		l = true
	}
	return
}
