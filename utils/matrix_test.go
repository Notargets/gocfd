package utils

import (
	"fmt"
	"math"
	"testing"

	"gonum.org/v1/gonum/mat"

	"github.com/stretchr/testify/assert"
)

func TestMatrix(t *testing.T) {
	// Basic Index utilities
	{
		nr, nc := 2, 3
		A := NewMatrix(nr, nc, []float64{0, 1, 2, 3, 4, 5})
		index := []int{0, 1, 2, 3, 4, 5}
		for _, ind := range index {
			i, j := indexToIJ(ind, nc)
			assert.Equal(t, A.At(i, j), float64(ind))
		}
		A = NewMatrix(nc, nr, []float64{0, 1, 2, 3, 4, 5})
		for _, ind := range index {
			i, j := indexToIJ(ind, nr)
			assert.Equal(t, A.At(i, j), float64(ind))
		}
	}
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
	// Equate
	/*
		Examples:
			A.Equate(Values, ":", MyIndex) // 2D index, uses all rows permuted with MyIndex columns values
			A.Equate(Values, "0:3", MyIndex) // Same with limited rows
			A.Equate(Values, MyIndex, ":") // Same reversed
			A.Equate(Values, B, ":") // Row index comes from data values in matrix B
			A.Equate(2, B, ":") 	// Equate indexed locations to a constant, example of constant promotion
			A.Equate(2, MyRowColumnIndex) 	// 1D indexed assignment using combined row+column index
			A.Equate(2, ":", ":", "0:3") 	// 3D indexed assignment
			A.Equate(2, ":", ":", ":", "0:3") 	// 4D indexed assignment, etc
	*/
	{
		A := NewMatrix(3, 3, []float64{
			0, 1, 2,
			3, 4, 5,
			6, 7, 8,
		})
		A.Equate(-1, ":", 0)
		assert.True(t, nearVec([]float64{
			-1.0000, 1.0000, 2.0000,
			-1.0000, 4.0000, 5.0000,
			-1.0000, 7.0000, 8.0000,
		}, A.Data(), 0.0001))

		A.Equate(-2, 0, ":")
		assert.True(t, nearVec([]float64{
			-2.0000, -2.0000, -2.0000,
			-1.0000, 4.0000, 5.0000,
			-1.0000, 7.0000, 8.0000,
		}, A.Data(), 0.0001))

		A.Equate(-3, 1, []float64{1})
		assert.True(t, nearVec([]float64{
			-2.0000, -2.0000, -2.0000,
			-1.0000, -3.0000, 5.0000,
			-1.0000, 7.0000, 8.0000,
		}, A.Data(), 0.0001))

		B := NewMatrix(1, 1, []float64{2})
		A.Equate(-4, 2, B)
		assert.True(t, nearVec([]float64{
			-2.0000, -2.0000, -2.0000,
			-1.0000, -3.0000, 5.0000,
			-1.0000, 7.0000, -4.0000,
		}, A.Data(), 0.0001))

		I := Index{1}
		A.Equate(-5, 2, I)
		assert.True(t, nearVec([]float64{
			-2.0000, -2.0000, -2.0000,
			-1.0000, -3.0000, 5.0000,
			-1.0000, -5.0000, -4.0000,
		}, A.Data(), 0.0001))

		A.Equate([]float64{
			0, 1, 2,
			3, 4, 5,
			6, 7, 8,
		}, ":", ":")
		assert.True(t, nearVec([]float64{
			0.0000, 1.0000, 2.0000,
			3.0000, 4.0000, 5.0000,
			6.0000, 7.0000, 8.0000,
		}, A.Data(), 0.0001))
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
