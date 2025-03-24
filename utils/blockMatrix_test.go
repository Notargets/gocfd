package utils

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestBlockMatrix(t *testing.T) {
	// [Scalar]: Test LU decomposition and solve
	{
		/*
				A = 1 2 3 4
					4 1 2 3
					3 4 1 2
					2 3 4 1
				Known solutions:
				    det(A) = -160
				Ainv =
			   -0.225  0.275  0.025  0.025
			    0.025 -0.225  0.275  0.025
			    0.025  0.025 -0.225  0.275
			    0.275  0.025  0.025 -0.225
		*/
		Bm := NewBlockMatrix(4, 4)
		Bm.M[0][0] = NewMatrix(1, 1, []float64{1.})
		Bm.M[0][1] = NewMatrix(1, 1, []float64{2.})
		Bm.M[0][2] = NewMatrix(1, 1, []float64{3.})
		Bm.M[0][3] = NewMatrix(1, 1, []float64{4.})

		Bm.M[1][0] = NewMatrix(1, 1, []float64{4.})
		Bm.M[1][1] = NewMatrix(1, 1, []float64{1.})
		Bm.M[1][2] = NewMatrix(1, 1, []float64{2.})
		Bm.M[1][3] = NewMatrix(1, 1, []float64{3.})

		Bm.M[2][0] = NewMatrix(1, 1, []float64{3.})
		Bm.M[2][1] = NewMatrix(1, 1, []float64{4.})
		Bm.M[2][2] = NewMatrix(1, 1, []float64{1.})
		Bm.M[2][3] = NewMatrix(1, 1, []float64{2.})

		Bm.M[3][0] = NewMatrix(1, 1, []float64{2.})
		Bm.M[3][1] = NewMatrix(1, 1, []float64{3.})
		Bm.M[3][2] = NewMatrix(1, 1, []float64{4.})
		Bm.M[3][3] = NewMatrix(1, 1, []float64{1.})
		BmOrig := Bm.Copy()
		// fmt.Printf(Bm.String())

		b := make([]Matrix, Bm.Nr)
		for i := 0; i < Bm.Nr; i++ {
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
		assert.InDeltaf(t, 0.5, x.M[0][0].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, 0.5, x.M[1][0].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, 0.5, x.M[2][0].DataP[0], 0.000001, "error msg %s")
		assert.InDeltaf(t, -0.5, x.M[3][0].DataP[0], 0.000001, "error msg %s")
		assert.Equal(t, []int{1, 2, 3, 0}, Bm.P) // Known permutation matrix, one swap per each row

		// Validate solution
		A := BmOrig.Mul(x) // Multiply original block matrix by solution to get b
		At := A.Transpose()
		for i := 0; i < len(b); i++ {
			assert.InDeltaf(t, b[i].DataP[0], A.M[0][i].DataP[0], 0.0000001, "err msg %s")
			assert.InDeltaf(t, b[i].DataP[0], At.M[i][0].DataP[0], 0.0000001, "err msg %s")
		}
		// Determinant
		{
			det, err := Bm.LUPDeterminant()
			assert.Nil(t, err)
			assert.InDeltaf(t, -160, det, 0.0000001, "err msg %s")
		}
		// Inverse
		{
			Ainv, err := Bm.LUPInvert()
			assert.Nil(t, err)
			N := Bm.Nr
			Binv := []float64{
				-0.225, 0.275, 0.025, 0.025,
				0.025, -0.225, 0.275, 0.025,
				0.025, 0.025, -0.225, 0.275,
				0.275, 0.025, 0.025, -0.225,
			}
			var ii int
			for i := 0; i < N; i++ {
				for j := 0; j < N; j++ {
					val := Ainv.M[i][j].DataP[0]
					assert.InDeltaf(t, Binv[ii], val, 0.0000001, "err msg %s")
					ii++
				}
			}
		}
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
		Bm.M[0][0] = one.Copy()
		Bm.M[0][1] = NewMatrix(1, 1, []float64{2.})
		Bm.M[0][2] = NewMatrix(1, 1, []float64{3.})
		Bm.M[0][3] = NewMatrix(1, 1, []float64{4.})

		Bm.M[1][0] = NewMatrix(1, 1, []float64{4.})
		Bm.M[1][1] = one.Copy()
		Bm.M[1][2] = NewMatrix(1, 1, []float64{2.})
		Bm.M[1][3] = NewMatrix(1, 1, []float64{3.})

		Bm.M[2][0] = NewMatrix(1, 1, []float64{3.})
		Bm.M[2][1] = NewMatrix(1, 1, []float64{4.})
		Bm.M[2][2] = one.Copy()
		Bm.M[2][3] = NewMatrix(1, 1, []float64{2.})

		Bm.M[3][0] = NewMatrix(1, 1, []float64{2.})
		Bm.M[3][1] = NewMatrix(1, 1, []float64{3.})
		Bm.M[3][2] = NewMatrix(1, 1, []float64{4.})
		Bm.M[3][3] = one.Copy()
		BmOrig := Bm.Copy()
		// fmt.Printf(Bm.String())

		b := make([]Matrix, Bm.Nr)
		for i := 0; i < Bm.Nr; i++ {
			// b[i] = one.Copy().Scale(float64(i + 1))
			val := float64(i + 1.)
			b[i] = NewMatrix(4, 1, []float64{val, val, val, val})
			// fmt.Printf(b[i].String("b[" + strconv.Itoa(i) + "]"))
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
		mul := 1.
		for i := 0; i < Bm.Nr; i++ {
			// msg := "x[" + strconv.Itoa(i) + "]"
			// fmt.Printf(x.M[i][0].String(msg))
			if i == 3 {
				mul = -1.
			}
			// Known answer x = [0.5[I],0.5[I],0.5[I],-0.5[I]]
			val := mul * 0.5
			assert.InDeltaSlicef(t, []float64{val, val, val, val}, x.M[i][0].DataP, 0.0000001, "err msg %s")
			// assert.InDeltaSlicef(t,
			// one.Copy().Scale(mul*0.5).DataP,
			// x.M[i][0].DataP, 0.0000001, "err msg %s")
		}
		assert.Equal(t, []int{1, 2, 3, 0}, Bm.P) // Known permutation matrix, one swap per each row

		// Bs := Bm.Mul(x)
		// fmt.Printf(Bs.String())

		// Validate solution
		A := BmOrig.Mul(x) // Multiply original block matrix by solution to get b
		At := A.Transpose()
		for i := 0; i < len(b); i++ {
			assert.InDeltaf(t, b[i].DataP[0], A.M[0][i].DataP[0], 0.0000001, "err msg %s")
			assert.InDeltaf(t, b[i].DataP[0], At.M[i][0].DataP[0], 0.0000001, "err msg %s")
		}
		// Inverse
		{
			Ainv, err := Bm.LUPInvert()
			assert.Nil(t, err)
			N := Bm.Nr
			Binv := []float64{
				-0.225, 0.275, 0.025, 0.025,
				0.025, -0.225, 0.275, 0.025,
				0.025, 0.025, -0.225, 0.275,
				0.275, 0.025, 0.025, -0.225,
			}
			var ii int
			for i := 0; i < N; i++ {
				for j := 0; j < N; j++ {
					val := Ainv.M[i][j].DataP[0]
					assert.InDeltaf(t, Binv[ii], val, 0.0000001, "err msg %s")
					ii++
				}
			}
		}
	}
}
