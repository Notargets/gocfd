package utils

import (
	"testing"

	"gonum.org/v1/gonum/mat"

	"github.com/stretchr/testify/assert"
)

func Test_QRFactorization(t *testing.T) {
	tol := 0.000001
	A := NewMatrix(3, 3, []float64{
		1, 2, 3,
		4, 5, 6,
		7, 8, 10,
	})
	if testing.Verbose() {
		A.Print("A")
	}
	AInv, _ := A.Inverse()
	Q, R := A.QRFactorization()
	AInvAlt := R.InverseWithCheck().Mul(Q.Transpose())
	assert.InDeltaSlicef(t, AInv.DataP, AInvAlt.DataP, tol, "")
	if testing.Verbose() {
		Q.Print("Q")
		R.Print("R")
		AInv.Print("A Inverse")
		AInvAlt.Print("Alt A Inverse")
	}
	AltA := Q.Mul(R)
	if testing.Verbose() {
		AltA.Print("A from QR")
	}
	assert.InDeltaSlicef(t, A.DataP, AltA.DataP, tol, "")
}

func TestMatrix(t *testing.T) {
	// Basic Index utilities
	{
		nr, nc := 2, 3
		A := NewMatrix(nr, nc, []float64{0, 1, 2, 3, 4, 5})
		B := NewMatrix(nr, nc, A.Data())
		index := []int{0, 1, 2, 3, 4, 5}
		for _, ind := range index {
			i, j := indexToIJ(ind, nc)
			assert.Equal(t, A.At(i, j), float64(ind))
			i, j = indexToIJColMajor(ind, nc)
			assert.Equal(t, B.Transpose().At(i, j), float64(ind))
		}
		A = NewMatrix(nc, nr, []float64{0, 1, 2, 3, 4, 5})
		B = NewMatrix(nr, nc, A.Data())
		for _, ind := range index {
			i, j := indexToIJ(ind, nr)
			assert.Equal(t, A.At(i, j), float64(ind))
			i, j = indexToIJColMajor(ind, nc)
			assert.Equal(t, B.Transpose().At(i, j), float64(ind))
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
		// fmt.Printf("Ainv = \n%v\n", mat.Formatted(Ainv, mat.Squeeze()))
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
		// fmt.Printf("Ainv = \n%v\n", mat.Formatted(Ainv, mat.Squeeze()))
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
		assert.Equal(t, A, NewMatrix(2, 3, []float64{
			1, 0, 3,
			4, 0, 6,
		}))
		A = M.Copy().SetRange(0, -1, -3, -3, 0)
		assert.Equal(t, A, NewMatrix(2, 3, []float64{
			0, 2, 3,
			0, 5, 6,
		}))
		A = M.Copy().SetRange(-1, -1, -2, -2, 0)
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
			Ainv.Equate(Values, ":", MyIndex) // 2D index, uses all rows permuted with MyIndex columns values
			Ainv.Equate(Values, "0:3", MyIndex) // Same with limited rows
			Ainv.Equate(Values, MyIndex, ":") // Same reversed
			Ainv.Equate(Values, B, ":") // Row index comes from data values in matrix B
			Ainv.Equate(2, B, ":") 	// Equate indexed locations to a constant, example of constant promotion
			Ainv.Equate(2, MyRowColumnIndex) 	// 1D indexed assignment using combined row+column index
			Ainv.Equate(2, ":", ":", "0:3") 	// 3D indexed assignment
			Ainv.Equate(2, ":", ":", ":", "0:3") 	// 4D indexed assignment, etc
	*/
	{
		A := NewMatrix(3, 3, []float64{
			0, 1, 2,
			3, 4, 5,
			6, 7, 8,
		})
		msg := "error msg %s"
		dlt := 0.0001
		R1 := []float64{
			-1.0000, 1.0000, 2.0000,
			-1.0000, 4.0000, 5.0000,
			-1.0000, 7.0000, 8.0000}
		R2 := []float64{
			-2.0000, -2.0000, -2.0000,
			-1.0000, 4.0000, 5.0000,
			-1.0000, 7.0000, 8.0000}
		R3 := []float64{
			-2.0000, -2.0000, -2.0000,
			-1.0000, -3.0000, 5.0000,
			-1.0000, 7.0000, 8.0000}
		R4 := []float64{
			-2.0000, -2.0000, -2.0000,
			-1.0000, -3.0000, 5.0000,
			-1.0000, 7.0000, -4.0000}
		R5 := []float64{
			-2.0000, -2.0000, -2.0000,
			-1.0000, -3.0000, 5.0000,
			-1.0000, -5.0000, -4.0000}
		{
			A.Equate(-1, ":", 0)
			assert.InDeltaSlicef(t, R1, A.DataP, dlt, msg)

			A.Equate(-2, 0, ":")
			assert.InDeltaSlicef(t, R2, A.DataP, dlt, msg)

			A.Equate(-3, 1, []float64{1})
			assert.InDeltaSlicef(t, R3, A.DataP, dlt, msg)

			B := NewMatrix(1, 1, []float64{2})
			A.Equate(-4, 2, B)
			assert.InDeltaSlicef(t, R4, A.DataP, dlt, msg)

			I := Index{1}
			A.Equate(-5, 2, I)
			assert.InDeltaSlicef(t, R5, A.DataP, dlt, msg)

			R := []float64{0, 1, 2, 3, 4, 5, 6, 7, 8}
			A.Equate(R, ":", ":")
			assert.InDeltaSlicef(t, R, A.DataP, dlt, msg)
		}
	}
	// Sparse Equate
	{
		nr, nc := 3, 3
		A := NewDOK(nr, nc)
		A.Equate([]float64{
			0, 1, 2,
			3, 4, 5,
			6, 7, 8,
		}, ":", ":")
		testRowMajor := func(A mat.Matrix) {
			check := []float64{
				0.0000, 1.0000, 2.0000,
				3.0000, 4.0000, 5.0000,
				6.0000, 7.0000, 8.0000,
			}
			// Check to ensure row-major traversal of values
			var ind int
			for j := 0; j < nc; j++ {
				for i := 0; i < nr; i++ {
					assert.Equal(t, check[ind], A.At(i, j))
					ind++
				}
			}
		}
		testRowMajor(A)

		B := A.ToCSR()
		B.Equate([]float64{
			0, 1, 2,
			3, 4, 5,
			6, 7, 8,
		}, ":", ":")
		testRowMajor(B)
	}
	{ // Test parallel multiply
		A := NewMatrix(2, 4, []float64{
			0., 1., 2., 3.,
			4., 5., 6., 7.,
		})
		B := NewMatrix(3, 2, []float64{
			10., 20.,
			30., 40.,
			50., 60.,
		})
		assert.Equal(t, B.Mul(A), B.MulParallel(A, 1))
		assert.Equal(t, B.Mul(A), B.MulParallel(A, 2))
		assert.Equal(t, B.Mul(A), B.MulParallel(A, 4))
		assert.Equal(t, B.Mul(A), B.MulParallel(A, 5))
		// fmt.Println(B.MulParallel(A, 2).Print("BmulA"))
	}
	// Matrix x Matrix
	{
		A := NewMatrix(2, 4, []float64{
			0., 1., 2., 3.,
			4., 5., 6., 7.,
		})
		B := NewMatrix(3, 2, []float64{
			10., 20.,
			30., 40.,
			50., 60.,
		})
		C := NewMatrix(3, 4, []float64{
			80., 110., 140., 170.,
			160., 230., 300., 370.,
			240., 350., 460., 570.,
		})

		assert.Equal(t, C, B.Mul(A))
		R := NewMatrix(3, 4)
		assert.Equal(t, C, B.Mul(A, R))
		assert.Equal(t, C, R)
	}
	// [Matrix,Scalar] <-> [Matrix,Scalar]: Add, Subtract, Mult
	{
		A := NewMatrix(2, 4, []float64{
			0., 1., 2., 3.,
			4., 5., 6., 7.,
		})
		B := NewMatrix(2, 4, []float64{
			10., 11., 12., 13.,
			14., 15., 16., 17.,
		})
		ApB := NewMatrix(2, 4, []float64{
			10., 12., 14., 16.,
			18., 20., 22., 24.,
		})
		Tens := NewMatrix(2, 4, []float64{
			10., 10., 10., 10.,
			10., 10., 10., 10.,
		})
		BmA := Tens.Copy()

		Tensminus100 := NewMatrix(2, 4, []float64{
			-90., 10., 10., 10.,
			10., -90., 10., 10.,
		})
		Tensplus100 := NewMatrix(2, 4, []float64{
			110., 10., 10., 10.,
			10., 110., 10., 10.,
		})
		s100minusTens := NewMatrix(2, 4, []float64{
			90., -10., -10., -10.,
			-10., 90., -10., -10.,
		})

		Bs100 := NewMatrix(1, 1, []float64{100.}) // Scalar matrix
		Cs3 := NewMatrix(1, 1, []float64{3.})     // Scalar matrix
		s100 := NewMatrix(1, 1, []float64{100.})
		s300 := NewMatrix(1, 1, []float64{300.})
		s2000000 := NewMatrix(1, 1, []float64{2000000.})

		// Matrix +/- Matrix
		{
			BB := B.Copy()
			assert.Equal(t, ApB, BB.Add(A))
			R := NewMatrix(2, 4)
			BB = B.Copy()
			assert.Equal(t, ApB, BB.Add(A, R))
			assert.Equal(t, ApB, R)

			BB = B.Copy()
			assert.Equal(t, BmA, BB.Subtract(A))
			R = NewMatrix(2, 4)
			BB = B.Copy()
			assert.Equal(t, BmA, BB.Subtract(A, R))
			assert.Equal(t, BmA, R)
		}
		// Scalar x Matrix, Matrix x Scalar
		{
			Ax100 := A.Copy().Scale(100)

			// [Matrix,Scalar] x [Matrix,Scalar]
			assert.Equal(t, Ax100, Bs100.Mul(A))
			assert.Equal(t, Ax100, A.Mul(Bs100))
			R := NewMatrix(2, 4)
			assert.Equal(t, Ax100, A.Mul(Bs100, R))
			assert.Equal(t, Ax100, R)
		}
		// [Scalar] x [Scalar]
		{
			R := NewMatrix(1, 1)
			assert.Equal(t, s300, Bs100.Mul(Cs3))
			assert.Equal(t, s300, Bs100.Mul(Cs3, R))
			assert.Equal(t, s300, R)
			R = NewMatrix(1, 1)
			assert.Equal(t, s300, Cs3.Mul(Bs100))
			assert.Equal(t, s300, Cs3.Mul(Bs100, R))
			assert.Equal(t, s300, R)

			Cs1000 := Cs3.Copy().Scale(0).AddScalar(1000.)
			Bs2000 := Bs100.Copy().Scale(0).AddScalar(2000.)
			R = NewMatrix(1, 1)
			assert.Equal(t, s2000000, Cs1000.Mul(Bs2000, R))
			assert.Equal(t, s2000000, R)
			R = NewMatrix(1, 1)
			assert.Equal(t, s2000000, Bs2000.Mul(Cs1000, R))
			assert.Equal(t, s2000000, R)
		}
		// [Matrix,Scalar] +/- [Matrix,Scalar], tests upscaling of scalar to diagonal matrix
		{
			// Matrix -/+ Scalar
			RR := Tens.Copy()
			R := NewMatrix(2, 4)
			assert.Equal(t, Tensminus100, RR.Subtract(s100, R)) // Place output in provided storage
			assert.Equal(t, Tensminus100, R)
			assert.Equal(t, Tensminus100, RR.Subtract(s100)) // Overwrite argument

			RR = Tens.Copy()
			R = NewMatrix(2, 4)
			assert.Equal(t, Tensplus100, RR.Add(s100, R)) // Place output in provided storage
			assert.Equal(t, Tensplus100, R)
			assert.Equal(t, Tensplus100, RR.Add(s100)) // Overwrite argument

			// Scalar -/+ Matrix
			RR = Tens.Copy()
			R = NewMatrix(2, 4)
			assert.Equal(t, s100minusTens, s100.Subtract(RR, R)) // Place output in provided storage
			assert.Equal(t, s100minusTens, R)
			R = s100.Subtract(RR)
			assert.Equal(t, s100minusTens, R)
			assert.Equal(t, s100minusTens, s100.Subtract(RR)) // Does not overwrite scalar argument

			RR = Tens.Copy()
			R = NewMatrix(2, 4)
			assert.Equal(t, Tensplus100, s100.Add(RR, R)) // Place output in provided storage
			assert.Equal(t, Tensplus100, R)
			R = s100.Add(RR)
			assert.Equal(t, Tensplus100, R)
			assert.Equal(t, Tensplus100, s100.Add(RR)) // Does not overwrite scalar argument

			// Scalar +/- Scalar
			assert.Equal(t, 103., s100.Add(Cs3).DataP[0])
			R = s100.Add(Cs3)
			assert.True(t, R.IsScalar())
			assert.Equal(t, 103., Cs3.Add(s100).DataP[0])
			R = Cs3.Add(s100)
			assert.True(t, R.IsScalar())

			assert.Equal(t, 97., s100.Subtract(Cs3).DataP[0])
			R = s100.Subtract(Cs3)
			assert.True(t, R.IsScalar())
			assert.Equal(t, -97., Cs3.Subtract(s100).DataP[0])
			R = Cs3.Subtract(s100)
			assert.True(t, R.IsScalar())
		}
	}
	// [Matrix,Scalar]: Inverse
	{
		A := NewMatrix(4, 4, []float64{
			1., 2., 3., 4.,
			4., 1., 2., 3.,
			3., 4., 1., 2.,
			2., 3., 4., 1.,
		})
		Ainv := NewMatrix(4, 4, []float64{
			-0.2250, 0.2750, 0.0250, 0.0250,
			0.0250, -0.2250, 0.2750, 0.0250,
			0.0250, 0.0250, -0.2250, 0.2750,
			0.2750, 0.0250, 0.0250, -0.2250,
		})
		s0 := NewMatrix(1, 1, []float64{0.})
		s100 := NewMatrix(1, 1, []float64{100.})
		// Matrix: Inverse
		{
			R, err := A.Inverse()
			assert.Nil(t, err)
			assert.InDeltaSlicef(t, Ainv.DataP, R.DataP, 0.00000001, "error msg %s")
		}
		// Scalar: Inverse
		{
			R, err := s100.Inverse()
			assert.Nil(t, err)
			assert.InDeltaf(t, 0.01, R.DataP[0], 0.00000001, "error msg %s")

			_, err = s0.Inverse()
			assert.NotNil(t, err)
		}
		// [Matrix,Scalar]: Determinant
		{
			assert.InDeltaf(t, -160., mat.Det(A), 0.00001, "err msg %s")
			assert.InDeltaf(t, -160., mat.Det(A.M), 0.00001, "err msg %s")
			assert.InDeltaf(t, 0., mat.Det(s0), 0.00001, "err msg %s")
			assert.InDeltaf(t, 0., mat.Det(s0.M), 0.00001, "err msg %s")
			assert.InDeltaf(t, 100., mat.Det(s100), 0.00001, "err msg %s")
			assert.InDeltaf(t, 100., mat.Det(s100.M), 0.00001, "err msg %s")
		}
	}
}

func TestMatrix_Inverse(t *testing.T) {
	// Scalar: Inverse
	{
		A := NewMatrix(6, 6, []float64{
			0.70711, -1.00000, 1.22474, -1.73205, 2.12132, 2.73861,
			0.70711, -0.40000, -0.24495, -0.69282, -0.00000, -0.21909,
			0.70711, 0.20000, -0.73485, 0.34641, 0.42426, -0.32863,
			0.70711, 0.80000, -0.24495, 1.38564, 3.39411, 2.40998,
			0.70711, 1.40000, 1.22474, 2.42487, 8.90955, 7.99675,
			0.70711, 2.00000, 3.67423, -0.00000, -0.00000, 0.00000,
		})
		Ainv, err := A.Inverse()
		Ainv.Print("Ainv")
		assert.Nil(t, err)
		I := A.Mul(Ainv)
		assert.InDeltaSlice(t, []float64{
			1, 0, 0, 0, 0, 0,
			0, 1, 0, 0, 0, 0,
			0, 0, 1, 0, 0, 0,
			0, 0, 0, 1, 0, 0,
			0, 0, 0, 0, 1, 0,
			0, 0, 0, 0, 0, 1,
		}, I.DataP, 0.000001, "blah")
	}
}
