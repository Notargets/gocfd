package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

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
	{ // Warpfactor
		N := 3
		Np := (N + 1) * (N + 2) / 2
		warpf := Warpfactor(N, utils.NewVector(Np, []float64{-1.0000, -0.3333, 0.3333, 1.0000, -0.6667, -0.0000, 0.6667, -0.3333, 0.3333, 0}))
		assert.True(t, nearVec([]float64{0, -0.1281, 0.1281, 0, -0.2562, 0.0000, 0.2562, -0.1281, 0.1281, 0.0000}, warpf, 0.0001))
	}
	{ // Nodes2D - distribution
		N := 3
		x, y := Nodes2D(N)
		assert.True(t, nearVec([]float64{-1.0000, -0.4472, 0.4472, 1.0000, -0.7236, -0.0000, 0.7236, -0.2764, 0.2764, 0}, x.Data(), 0.0001))
		assert.True(t, nearVec([]float64{-0.5774, -0.5774, -0.5774, -0.5774, -0.0986, -0.0000, -0.0986, 0.6760, 0.6760, 1.1547}, y.Data(), 0.0001))
		N = 4
		x, y = Nodes2D(N)
		assert.True(t, nearVec([]float64{-1.0000, -0.6547, -0.0000, 0.6547, 1.0000, -0.8273, -0.3274, 0.3274, 0.8273, -0.5000, 0.0000, 0.5000, -0.1727, 0.1727, 0}, x.Data(), 0.0001))
		assert.True(t, nearVec([]float64{-0.5774, -0.5774, -0.5774, -0.5774, -0.5774, -0.2783, -0.1890, -0.1890, -0.2783, 0.2887, 0.3780, 0.2887, 0.8556, 0.8556, 1.1547}, y.Data(), 0.0001))
	}

	{ // Read file to test specific metrics
		// Check N = 1 briefly
		el := NewElements2D(1, "fstepA001.neu", false)
		assert.True(t, nearVec([]float64{
			2.5000, 0.5000, -1.5000, -1.5000, 2.5000, 0.5000,
			0.5000, 2.5000, 2.5000, 0.5000, -1.5000, -1.5000,
			-1.5000, -1.5000, 0.5000, 2.5000, 0.5000, 2.5000,
		}, el.LIFT.Data(), 0.0001))

		// Check N = 2 by comparison with Matlab code
		el = NewElements2D(2, "fstepA001.neu", false)
		assert.True(t, nearVec([]float64{
			4.5000, 2.0000, -0.5000, 1.0000, 4.0000, 1.0000, 4.5000, 2.0000, -0.5000,
			0.5000, 5.0000, 0.5000, -0.6250, -1.5000, 0.6250, -0.6250, -1.5000, 0.6250,
			-0.5000, 2.0000, 4.5000, 4.5000, 2.0000, -0.5000, 1.0000, 4.0000, 1.0000,
			-0.6250, -1.5000, 0.6250, 0.6250, -1.5000, -0.6250, 0.5000, 5.0000, 0.5000,
			0.6250, -1.5000, -0.6250, 0.5000, 5.0000, 0.5000, 0.6250, -1.5000, -0.6250,
			1.0000, 4.0000, 1.0000, -0.5000, 2.0000, 4.5000, -0.5000, 2.0000, 4.5000,
		}, el.LIFT.Data(), 0.0001))
		subset := utils.NewR2(el.Rx.Dims())
		assert.True(t, nearVec([]float64{
			2.7337, -10.4712, 15.5416, 19.2110,
			2.7337, -10.4712, 15.5416, 19.2110,
			2.7337, -10.4712, 15.5416, 19.2110,
			2.7337, -10.4712, 15.5416, 19.2110,
			2.7337, -10.4712, 15.5416, 19.2110,
			2.7337, -10.4712, 15.5416, 19.2110,
		}, el.Rx.Subset(subset.Range(":", "0:4"), 6, 4).Data(), 0.0001))
		assert.True(t, nearVec([]float64{
			-15.9949, 10.7477, 5.2533, -12.0803,
			-15.9949, 10.7477, 5.2533, -12.0803,
			-15.9949, 10.7477, 5.2533, -12.0803,
			-15.9949, 10.7477, 5.2533, -12.0803,
			-15.9949, 10.7477, 5.2533, -12.0803,
			-15.9949, 10.7477, 5.2533, -12.0803,
		}, el.Ry.Subset(subset.Range(":", "0:4"), 6, 4).Data(), 0.0001))
		assert.True(t, nearVec([]float64{
			14.1324, -14.8689, -12.2718, -1.3128,
			14.1324, -14.8689, -12.2718, -1.3128,
			14.1324, -14.8689, -12.2718, -1.3128,
			14.1324, -14.8689, -12.2718, -1.3128,
			14.1324, -14.8689, -12.2718, -1.3128,
			14.1324, -14.8689, -12.2718, -1.3128,
		}, el.Sx.Subset(subset.Range(":", "0:4"), 6, 4).Data(), 0.0001))
		assert.True(t, nearVec([]float64{
			4.4896, -10.7477, 10.6892, 16.5204,
			4.4896, -10.7477, 10.6892, 16.5204,
			4.4896, -10.7477, 10.6892, 16.5204,
			4.4896, -10.7477, 10.6892, 16.5204,
			4.4896, -10.7477, 10.6892, 16.5204,
			4.4896, -10.7477, 10.6892, 16.5204,
		}, el.Sy.Subset(subset.Range(":", "0:4"), 6, 4).Data(), 0.0001))
		subsetFacePts := utils.NewR2(el.Nfp*el.NFaces, el.K)
		assert.True(t, nearVec([]float64{
			-0.9531, 0.8104, 0.7541, 0.0792,
			-0.9531, 0.8104, 0.7541, 0.0792,
			-0.9531, 0.8104, 0.7541, 0.0792,
			0.8261, -1.0000, 0.2009, 0.9706,
			0.8261, -1.0000, 0.2009, 0.9706,
			0.8261, -1.0000, 0.2009, 0.9706,
			-0.1685, 0.6978, -0.9473, -0.8465,
			-0.1685, 0.6978, -0.9473, -0.8465,
			-0.1685, 0.6978, -0.9473, -0.8465,
		}, el.NX.Subset(subsetFacePts.Range(":", "0:4"), 9, 4).Data(), 0.0001))
		assert.True(t, nearVec([]float64{
			-0.3028, 0.5858, -0.6568, -0.9969,
			-0.3028, 0.5858, -0.6568, -0.9969,
			-0.3028, 0.5858, -0.6568, -0.9969,
			-0.5635, 0.0000, 0.9796, 0.2408,
			-0.5635, 0.0000, 0.9796, 0.2408,
			-0.5635, 0.0000, 0.9796, 0.2408,
			0.9857, -0.7163, -0.3202, 0.5323,
			0.9857, -0.7163, -0.3202, 0.5323,
			0.9857, -0.7163, -0.3202, 0.5323,
		}, el.NY.Subset(subsetFacePts.Range(":", "0:4"), 9, 4).Data(), 0.0001))
	}
	{ //Test Interpolation
		el := NewElements2D(1, "test_tris_1.neu", false)
		s := make([]float64, el.Np)
		for i := 0; i < el.Np; i++ {
			s[i] = float64(2 * i)
		}
		for i, rVal := range el.R.Data() {
			//TODO: Implement affine mapping for input coordinates to interpolation
			sVal := el.S.Data()[i]
			sInterp := el.Simplex2DInterpolate(rVal, el.S.Data()[i], s)
			fmt.Printf("fInterp[%8.5f,%8.5f] = %8.5f\n", rVal, sVal, sInterp)
			assert.True(t, near(s[i], sInterp, 0.00001))
		}
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
	bound := math.Max(tol, tol*math.Abs(a))
	if math.Abs(a-b) <= bound {
		l = true
	}
	return
}
