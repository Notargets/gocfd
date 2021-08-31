package utils

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestBlockMatrix(t *testing.T) {
	// [Scalar]: Test LU decomposition and solve
	{
		Bm := NewBlockMatrix(4, 1)
		Bm.A[0][0] = NewMatrix(1, 1, []float64{1.})
		Bm.A[0][1] = NewMatrix(1, 1, []float64{2.})
		Bm.A[0][2] = NewMatrix(1, 1, []float64{3.})
		Bm.A[0][3] = NewMatrix(1, 1, []float64{4.})

		Bm.A[1][0] = NewMatrix(1, 1, []float64{4.})
		Bm.A[1][1] = NewMatrix(1, 1, []float64{1.})
		Bm.A[1][2] = NewMatrix(1, 1, []float64{2.})
		Bm.A[1][3] = NewMatrix(1, 1, []float64{3.})

		Bm.A[2][0] = NewMatrix(1, 1, []float64{3.})
		Bm.A[2][1] = NewMatrix(1, 1, []float64{4.})
		Bm.A[2][2] = NewMatrix(1, 1, []float64{1.})
		Bm.A[2][3] = NewMatrix(1, 1, []float64{2.})

		Bm.A[3][0] = NewMatrix(1, 1, []float64{2.})
		Bm.A[3][1] = NewMatrix(1, 1, []float64{3.})
		Bm.A[3][2] = NewMatrix(1, 1, []float64{4.})
		Bm.A[3][3] = NewMatrix(1, 1, []float64{1.})
		//fmt.Printf(Bm.Print())

		b := make([]Matrix, Bm.N)
		for i := 0; i < Bm.N; i++ {
			b[i] = NewMatrix(1, 1, []float64{float64(i + 1)})
		}
		// Call LUPSolve without first calling LUPDecompose, expect an error
		x, err := Bm.LUPSolve(b)
		assert.NotNil(t, err)

		err = Bm.LUPDecompose()
		assert.Nil(t, err)
		assert.True(t, len(Bm.P) != 0)

		// Call LUPDecompose again, expect an error
		err = Bm.LUPDecompose()
		assert.NotNil(t, err)

		x, err = Bm.LUPSolve(b)
		assert.Nil(t, err)
		// Known answer x = [0.5,0.5,0.5,-0.5]
		assert.InDeltaf(t, 0.5, x[0].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, 0.5, x[1].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, 0.5, x[2].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, -0.5, x[3].DataP[0], 0.000001, "error msg %s")
		assert.Equal(t, []int{1, 2, 3, 0}, Bm.P) // Known permutation matrix, one swap per each row
	}
	// [Matrix]: Test LU decomposition and solve
	{
		one := NewMatrix(4, 4, []float64{
			1., 0., 0., 0.,
			0., 1., 0., 0.,
			0., 0., 1., 0.,
			0., 0., 0., 1.,
		})
		Bm := NewBlockMatrix(4, 4)
		Bm.A[0][0] = one.Copy()
		Bm.A[0][1] = NewMatrix(1, 1, []float64{2.})
		Bm.A[0][2] = NewMatrix(1, 1, []float64{3.})
		Bm.A[0][3] = NewMatrix(1, 1, []float64{4.})

		Bm.A[1][0] = NewMatrix(1, 1, []float64{4.})
		Bm.A[1][1] = one.Copy()
		Bm.A[1][2] = NewMatrix(1, 1, []float64{2.})
		Bm.A[1][3] = NewMatrix(1, 1, []float64{3.})

		Bm.A[2][0] = NewMatrix(1, 1, []float64{3.})
		Bm.A[2][1] = NewMatrix(1, 1, []float64{4.})
		Bm.A[2][2] = one.Copy()
		Bm.A[2][3] = NewMatrix(1, 1, []float64{2.})

		Bm.A[3][0] = NewMatrix(1, 1, []float64{2.})
		Bm.A[3][1] = NewMatrix(1, 1, []float64{3.})
		Bm.A[3][2] = NewMatrix(1, 1, []float64{4.})
		Bm.A[3][3] = one.Copy()
		//fmt.Printf(Bm.Print())

		b := make([]Matrix, Bm.N)
		for i := 0; i < Bm.N; i++ {
			b[i] = one.Copy().Scale(float64(i + 1))
		}
		// Call LUPSolve without first calling LUPDecompose, expect an error
		x, err := Bm.LUPSolve(b)
		assert.NotNil(t, err)

		err = Bm.LUPDecompose()
		assert.Nil(t, err)
		assert.True(t, len(Bm.P) != 0)

		// Call LUPDecompose again, expect an error
		err = Bm.LUPDecompose()
		assert.NotNil(t, err)

		x, err = Bm.LUPSolve(b)
		assert.Nil(t, err)
		// Known answer x = [0.5,0.5,0.5,-0.5]
		assert.InDeltaf(t, 0.5, x[0].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, 0.5, x[1].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, 0.5, x[2].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, -0.5, x[3].DataP[0], 0.000001, "error msg %s")
		mul := 1.
		for i := 0; i < Bm.N; i++ {
			//msg := "x[" + strconv.Itoa(i) + "]"
			//fmt.Printf(x[i].Print(msg))
			if i == 3 {
				mul = -1.
			}
			assert.InDeltaSlicef(t,
				one.Copy().Scale(mul*0.5).DataP,
				x[i].DataP, 0.0000001, "err msg %s")
		}
		assert.Equal(t, []int{1, 2, 3, 0}, Bm.P) // Known permutation matrix, one swap per each row
	}
}
