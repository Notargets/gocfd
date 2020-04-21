package DG1D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/stretchr/testify/assert"
)

func TestElements1D(t *testing.T) {
	/*
		type Elements1D struct {
			K, Np, Nfp, NFaces                int
			VX, FMask                         utils.Vector
			EToV, EToE, EToF                  utils.Matrix
			X, Dr, Rx, FScale, NX, LIFT       utils.Matrix
			V, Vinv                           utils.Matrix
			VmapM, VmapP, VmapB, VmapI, VmapO utils.Index
			MapB, MapI, MapO                  utils.Index
	*/
	{
		K := 4
		N := 3
		VX, EToV := SimpleMesh1D(0, 2, K)

		var el *Elements1D
		el = NewElements1D(N, VX, EToV)
		// Verify X
		assert.True(t, near(el.X.At(0, 1), 0.5))
		assert.True(t, near(el.X.At(3, 1), 1.0))
		assert.True(t, near(el.X.At(3, 2), 1.5))
		assert.True(t, near(el.X.At(2, 3), 1.8618033988))
		assert.True(t, near(el.X.At(1, 1), 0.6381966011))
		assert.True(t, near(el.X.SumCols().AtVec(0), 1))
		assert.True(t, near(el.X.SumRows().AtVec(0), 3))
		assert.True(t, near(el.X.SumRows().AtVec(3), 5))

		// Verify LIFT
		assert.True(t, near(el.LIFT.SumRows().AtVec(0), 6))
		assert.True(t, near(el.LIFT.SumRows().AtVec(3), 6))
		assert.True(t, near(el.LIFT.At(2, 0), 0.8944271909))
		assert.True(t, near(el.LIFT.At(2, 1), -0.8944271909))
		assert.True(t, near(el.LIFT.At(1, 0), -0.8944271909))
		assert.True(t, near(el.LIFT.At(1, 1), 0.8944271909))

		// Verify VmapM
		// Row-major
		//assert.Equal(t, utils.Index{0, 3, 4, 7, 8, 11, 12, 15}, el.VmapM)
		//assert.Equal(t, utils.Index{0, 4, 3, 8, 7, 12, 11, 15}, el.VmapP)
		assert.Equal(t, utils.Index{0, 12, 1, 13, 2, 14, 3, 15}, el.VmapM)
		assert.Equal(t, utils.Index{0, 1, 12, 2, 13, 3, 14, 15}, el.VmapP)
	}
}

func near(a, b float64) (l bool) {
	if math.Abs(a-b) < 1.e-08*math.Abs(a) {
		l = true
	}
	return
}
