package utils

import (
	"fmt"
	"testing"

	"gonum.org/v1/gonum/mat"

	"github.com/stretchr/testify/require"

	"github.com/stretchr/testify/assert"
)

func TestVector(t *testing.T) {
	assert.Equal(t, 123, 123, "should be equal")
	/*
		// x = ones(Np)*VX(va) + 0.5*(r+1.)*sT(vc);
		mm := utils.NewVector(Np).Set(1).ToMatrix().Mul(VX.SubsetIndex(va).Transpose())
		r := utils.Vector{mat.VecDenseCopyOf(R)}
		R = r.AddScalar(1).Scale(0.5).ToMatrix().Mul(sT.Transpose()).AddCopy(mm).M
	*/
	N := 3
	v1 := NewVector(N).Set(1)
	require.Equal(t, 1., v1.V.RawVector().Data[N-1])
	v1.Set(2)
	require.Equal(t, 2., v1.V.RawVector().Data[N-1])

	M := 2
	v2 := NewVector(M).Set(3)
	A := v1.ToMatrix().Mul(v2.Transpose())
	fmt.Printf("Ainv = \n%v\n", mat.Formatted(A, mat.Squeeze()))
	nr, nc := A.Dims()
	require.Equal(t, N, nr)
	require.Equal(t, M, nc)

	v1.V.SetVec(0, 1)
	v1.V.SetVec(1, 2)
	v1.V.SetVec(2, 3)
	v2.V.SetVec(0, 2)
	A = v1.ToMatrix().Mul(v2.Transpose())
	/*
		Ainv =
		⎡2  3⎤
		⎢4  6⎥
		⎣6  9⎦
	*/
	vec := []float64{2, 3, 4, 6, 6, 9} // Column major order
	fmt.Printf("v1, v2 = \n%v\n%v\n", mat.Formatted(v1, mat.Squeeze()), mat.Formatted(v2, mat.Squeeze()))
	fmt.Printf("Ainv = \n%v\n", mat.Formatted(A, mat.Squeeze()))
	require.Equal(t, vec, A.RawMatrix().Data)

	B := v1.Mul(v2)
	require.Equal(t, vec, B.RawMatrix().Data)
	// Linspace
	{
		req := NewVector(2).Linspace(-1, 1)
		assert.Equal(t, -1., req.AtVec(0))
		assert.Equal(t, 1., req.AtVec(1))
		req = NewVector(3).Linspace(-1, 1)
		assert.Equal(t, -1., req.AtVec(0))
		assert.Equal(t, 0., req.AtVec(1))
		assert.Equal(t, 1., req.AtVec(2))
	}
}
