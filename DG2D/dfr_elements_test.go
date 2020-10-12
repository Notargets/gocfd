package DG2D

import (
	"testing"

	"github.com/notargets/gocfd/utils"
	"github.com/stretchr/testify/assert"
)

func TestDFR2D(t *testing.T) {
	{ //Test Interpolation
		N := 1
		el := NewNDG2D(N, "test_tris_1.neu", false)
		s := make([]float64, el.Np)
		for i := 0; i < el.Np; i++ {
			s[i] = float64(2 * i)
		}
		// For each nodal location, interpolate a value (should equal the nodal function value)
		for i, rVal := range el.R.Data() {
			sVal := el.S.Data()[i]
			sInterp := el.Simplex2DInterpolate(rVal, sVal, s)
			//fmt.Printf("fInterp[%8.5f,%8.5f] = %8.5f\n", rVal, sVal, sInterp)
			assert.True(t, near(s[i], sInterp, 0.00001))
		}
		// Build an interpolating polynomial matrix using the nodal geometry
		interpM := el.Simplex2DInterpolatingPolyMatrix(el.R, el.S)
		values := interpM.Mul(utils.NewMatrix(el.Np, 1, s))
		assert.True(t, nearVec(s, values.Data(), 0.0000001))
	}
}
