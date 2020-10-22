package DG1D

import (
	"fmt"
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
		assert.Equal(t, utils.Index{0, 1, 2, 3, 12, 13, 14, 15}, el.VmapM)
		assert.Equal(t, utils.Index{0, 12, 13, 14, 1, 2, 3, 15}, el.VmapP)
		/*
			Number of Elements, NFaces =  2 4
			VmapM.Matrix =
			⎡ 0   1   2   3⎤
			⎣12  13  14  15⎦
			VmapM.Matrix.Data = [0 1 2 3 12 13 14 15]
			VmapP.Matrix =
			⎡0  12  13  14⎤
			⎣1   2   3  15⎦
			VmapP.Matrix.Data = [0 12 13 14 1 2 3 15]
		*/
		/*
			fmt.Printf("VmapM.Matrix = \n%v\n", mat.Formatted(el.VmapM.ToMatrixReversed(NF, el.K), mat.Squeeze()))
			fmt.Printf("VmapM.Matrix.Data = %v\n", el.VmapM.ToMatrixReversed(NF, el.K).RawMatrix().Data)
			fmt.Printf("VmapP.Matrix = \n%v\n", mat.Formatted(el.VmapP.ToMatrixReversed(NF, el.K), mat.Squeeze()))
			fmt.Printf("VmapP.Matrix.Data = %v\n", el.VmapP.ToMatrixReversed(NF, el.K).RawMatrix().Data)
		*/
	}
	/*
		Check face mapping
	*/
	{
		K := 4
		for N := 1; N < 10; N++ {
			VX, EToV := SimpleMesh1D(0, 2, K)
			var el *Elements1D
			el = NewElements1D(N, VX, EToV)
			facesMX := el.X.Subset(el.VmapM, 2, K)
			facesPX := el.X.Subset(el.VmapP, 2, K)
			assert.Equal(t, facesMX, facesPX)
		}
	}
	/*
		Check interpolation
	*/
	{
		K := 1
		N := 2
		VX, EToV := SimpleMesh1D(0, 2, K)
		var el, elS *Elements1D
		el = NewElements1D(N+2, VX, EToV)
		elS = NewElements1D(N, VX, EToV, GAUSS)
		tolTight := 1.e-08
		tolLoose := 1.e-03
		li_left_edge := el.LagrangeInterpolant(-1)
		assert.True(t, nearVec(li_left_edge.RawMatrix().Data, []float64{1, 0, 0, 0, 0}, tolTight))
		li_right_edge := el.LagrangeInterpolant(1)
		assert.True(t, nearVec(li_right_edge.RawMatrix().Data, []float64{0, 0, 0, 0, 1}, tolTight))
		assert.True(t, nearVec(elS.LagrangeInterpolant(-1).RawMatrix().Data, []float64{1.9304, -1.3333, 0.4029}, tolLoose))
		assert.True(t, nearVec(elS.LagrangeInterpolant(1).RawMatrix().Data, []float64{0.4029, -1.3333, 1.9304}, tolLoose))

		N = 1
		VX, EToV = SimpleMesh1D(0, 2, K)
		el = NewElements1D(N+2, VX, EToV)
		elS = NewElements1D(N, VX, EToV, GAUSS)
		li_left_edge = elS.LagrangeInterpolant(-1)
		li_right_edge = elS.LagrangeInterpolant(1)
		assert.True(t, nearVec(li_left_edge.RawMatrix().Data, []float64{1.6180, -0.6180}, tolLoose))
		assert.True(t, nearVec(li_right_edge.RawMatrix().Data, []float64{-0.6180, 1.6180}, tolLoose))
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
