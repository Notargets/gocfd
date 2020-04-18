package utils

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestMatrix(t *testing.T) {
	// Transpose
	{
		M := NewMatrix(2, 3, []float64{
			1, 2, 3,
			4, 5, 6,
		})
		mNr, mNc := M.Dims()
		A := M.Transpose()
		aNr, aNc := A.Dims()
		assert.Equal(t, aNc, mNr)
		assert.Equal(t, aNr, mNc)
		assert.Equal(t, A.RawMatrix().Data, []float64{1, 4, 2, 5, 3, 6})
	}
	// SliceRows
	{
		M := NewMatrix(2, 3, []float64{
			1, 2, 3,
			4, 5, 6,
		})
		I := NewIndex(2)
		I[0] = 1
		I[1] = 0
		A := M.SliceRows(I)
		//fmt.Printf("A = \n%v\n", mat.Formatted(A, mat.Squeeze()))
		assert.Equal(t, A, NewMatrix(2, 3, []float64{
			4, 5, 6,
			1, 2, 3,
		}))
	}
	// SliceCols
	{
		M := NewMatrix(2, 3, []float64{
			1, 2, 3,
			4, 5, 6,
		})
		I := NewIndex(2)
		I[0] = 1
		I[1] = 0
		A := M.SliceCols(I)
		//fmt.Printf("A = \n%v\n", mat.Formatted(A, mat.Squeeze()))
		assert.Equal(t, A, NewMatrix(2, 2, []float64{
			2, 1,
			5, 4,
		}))
	}
}
