package utils

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestMatrix(t *testing.T) {
	// Transpose
	{
		M := NewMatrix(2, 3, []float64{1, 2, 3, 4, 5, 6})
		mNr, mNc := M.Dims()
		A := M.Transpose()
		aNr, aNc := A.Dims()
		assert.Equal(t, aNc, mNr)
		assert.Equal(t, aNr, mNc)
		assert.Equal(t, A.RawMatrix().Data, []float64{1, 4, 2, 5, 3, 6})
	}
}
