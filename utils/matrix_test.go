package utils

import (
	"fmt"
	"testing"

	"gonum.org/v1/gonum/mat"

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
	// SetRange
	{
		M := NewMatrix(2, 3, []float64{
			1, 2, 3,
			4, 5, 6,
		})
		A := M.Copy().SetRange(0, -1, -2, -2, 0)
		fmt.Printf("A = \n%v\n", mat.Formatted(A, mat.Squeeze()))
		assert.Equal(t, A, NewMatrix(2, 3, []float64{
			1, 0, 3,
			4, 0, 6,
		}))
		A = M.Copy().SetRange(0, -1, -3, -3, 0)
		fmt.Printf("A = \n%v\n", mat.Formatted(A, mat.Squeeze()))
		assert.Equal(t, A, NewMatrix(2, 3, []float64{
			0, 2, 3,
			0, 5, 6,
		}))
		A = M.Copy().SetRange(-1, -1, -2, -2, 0)
		fmt.Printf("A = \n%v\n", mat.Formatted(A, mat.Squeeze()))
		assert.Equal(t, A, NewMatrix(2, 3, []float64{
			1, 2, 3,
			4, 0, 6,
		}))
	}
	// Sum
	{
		M := NewMatrix(2, 3, []float64{
			1, 2, 3,
			4, 5, 6,
		})
		V := M.SumRows()
		assert.Equal(t, V, NewVector(2, []float64{6, 15}))
		V = M.SumCols()
		assert.Equal(t, V, NewVector(3, []float64{5, 7, 9}))
	}
	// LU Solve
	{
		A := NewMatrix(3, 3, []float64{
			2, 1, 1,
			-1, 1, -1,
			1, 2, 3,
		})
		B := NewMatrix(3, 1, []float64{
			2,
			3,
			-10,
		})
		X := A.LUSolve(B)
		xd := X.Data()
		assert.Equal(t, []float64{3, 1, -5}, xd)
	}
}
