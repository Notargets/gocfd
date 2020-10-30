package DG2D

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/utils"
	"github.com/stretchr/testify/assert"
)

func TestDFR2D(t *testing.T) {
	{ // Test Interpolation
		N := 2
		dfr := NewDFR2D(N)
		el := dfr.SolutionElement
		s := make([]float64, el.Np)
		for i := 0; i < el.Np; i++ {
			s[i] = float64(2 * i)
		}
		// For each nodal location, interpolate a value (should equal the nodal function value)
		// Build an interpolating polynomial matrix using the nodal geometry
		interpM := el.Simplex2DInterpolatingPolyMatrix(el.R, el.S)
		values := interpM.Mul(utils.NewMatrix(el.Np, 1, s))
		assert.True(t, nearVec(s, values.Data(), 0.0000001))
	}
	{ // Test point distribution
		N := 1
		dfr := NewDFR2D(N)
		el := dfr.SolutionElement
		rt := dfr.FluxElement
		assert.True(t, nearVec(rt.GetInternalLocations(rt.R), el.R.Data(), 0.000001))
		assert.True(t, nearVec(rt.GetInternalLocations(rt.S), el.S.Data(), 0.000001))
		//fmt.Printf("Edge R = %v\n", rt.GetEdgeLocations(rt.R))
		//fmt.Printf("Edge S = %v\n", rt.GetEdgeLocations(rt.S))
	}
	{ // Test interpolation from solution points to flux points
		N := 1
		dfr := NewDFR2D(N)
		el := dfr.SolutionElement
		rt := dfr.FluxElement
		assert.Equal(t, el.Np, rt.Nint)
		solution := make([]float64, rt.Nint)
		// Load some values into the solution space
		for i := 0; i < rt.Nint; i++ {
			solution[i] = float64(i) / float64(rt.Nint-1)
		}
		// Interpolate from interior to flux points
		sV := utils.NewMatrix(rt.Nint, 1, solution)
		_ = dfr.FluxInterpMatrix.Mul(sV)
		//fmt.Printf("%s\n", fluxInterp.Print("fluxInterp"))
		//fmt.Printf("%s\n", sV.Print("sV"))
	}
	{ // Test packed int for edge labeling
		en := NewEdgeNumber([2]int{1, 0})
		assert.Equal(t, EdgeNumber(1<<32), en)
		assert.Equal(t, [2]int{0, 1}, en.GetVertices())

		en = NewEdgeNumber([2]int{0, 1})
		assert.Equal(t, EdgeNumber(1<<32), en)
		assert.Equal(t, [2]int{0, 1}, en.GetVertices())

		en = NewEdgeNumber([2]int{0, 10})
		assert.Equal(t, EdgeNumber(10*(1<<32)), en)
		assert.Equal(t, [2]int{0, 10}, en.GetVertices())

		en = NewEdgeNumber([2]int{100, 0})
		assert.Equal(t, EdgeNumber(100*(1<<32)), en)
		assert.Equal(t, [2]int{0, 100}, en.GetVertices())

		en = NewEdgeNumber([2]int{100, 1})
		assert.Equal(t, EdgeNumber(100*(1<<32)+1), en)
		assert.Equal(t, [2]int{1, 100}, en.GetVertices())

		en = NewEdgeNumber([2]int{100, 100001})
		assert.Equal(t, EdgeNumber(100001*(1<<32)+100), en)
		assert.Equal(t, [2]int{100, 100001}, en.GetVertices())

		// Test maximum/minimum indices
		en = NewEdgeNumber([2]int{1, 1<<32 - 1})
		assert.Equal(t, EdgeNumber((1<<32-1)<<32+1), en)
		assert.Equal(t, [2]int{1, 1<<32 - 1}, en.GetVertices())

		en = NewEdgeNumber([2]int{1<<32 - 1, 1<<32 - 1})
		assert.Equal(t, EdgeNumber(1<<64-1), en)
		assert.Equal(t, [2]int{1<<32 - 1, 1<<32 - 1}, en.GetVertices())

		en = NewEdgeNumber([2]int{1<<32 - 1, 1})
		assert.Equal(t, EdgeNumber((1<<32-1)<<32+1), en)
		assert.Equal(t, [2]int{1, 1<<32 - 1}, en.GetVertices())
	}
	{ // Test triangulation
		N := 1
		dfr := NewDFR2D(N, "test_tris_5.neu")
		//PlotMesh(dfr.VX, dfr.VY, dfr.Tris.EToV, dfr.BCType, dfr.FluxX, dfr.FluxY, true)
		//utils.SleepFor(50000)
		//dfr := NewDFR2D(N, "fstepA001.neu")
		trn := dfr.Tris
		// Check against known answers for this case
		en := NewEdgeNumber([2]int{0, 2})
		e := trn.Edges[en]
		assert.Equal(t, uint8(2), e.NumConnectedTris)
		assert.Equal(t, uint32(0), e.ConnectedTris[0])
		assert.Equal(t, InternalEdgeNumber(2), e.ConnectedTriEdgeNumber[0])
		assert.Equal(t, uint32(1), e.ConnectedTris[1])
		assert.Equal(t, InternalEdgeNumber(0), e.ConnectedTriEdgeNumber[1])

		en = NewEdgeNumber([2]int{0, 1})
		e = trn.Edges[en]
		assert.Equal(t, BC_In, e.BCType)

		en = NewEdgeNumber([2]int{2, 3})
		e = trn.Edges[en]
		assert.Equal(t, BC_In, e.BCType)

		// Test Piola transform and jacobian
		for en, e := range dfr.Tris.Edges {
			for ii, triNum := range e.ConnectedTris {
				x1, x2 := dfr.Tris.GetEdgeCoordinates(en, e, triNum, dfr.VX, dfr.VY)
				switch e.ConnectedTriEdgeNumber[ii] {
				case Third:
					dx, dy := x2[0]-x1[0], x2[1]-x1[1]
					nx, ny := -dy, dx
					normal := [2]float64{nx, ny}
					J := dfr.J.Row(int(triNum)).Data()[0:4]
					Jdet := dfr.Jdet.Row(int(triNum)).Data()[0]
					nxT := dfr.PiolaTransform(J, Jdet, normal)
					Jinv := dfr.Jinv.Row(int(triNum)).Data()[0:4]
					nxTT := dfr.PiolaTransform(Jinv, 1./Jdet, nxT)
					switch triNum {
					case 0:
						assert.Equal(t, [2]float64{1., -.5}, normal)
						assert.Equal(t, [2]float64{2, 0}, nxT)
					case 1:
						assert.Equal(t, [2]float64{-1., -1.}, normal)
						assert.Equal(t, [2]float64{-2, 0}, nxT)
					}
					assert.Equal(t, normal, nxTT)
				}
			}
		}
	}
	{ // Test divergence
		errorCheck := func(dfr *DFR2D, k int, div, divCheck utils.Matrix) (min, max float64) {
			var (
				error float64
				Np    = dfr.FluxElement.Np
				dcD   = divCheck.Row(k).Data()[0:Np]
			)
			fmt.Println(div.Transpose().Print("div"))
			fmt.Println(divCheck.Row(k).Transpose().Print("divCheck"))
			for i, divVal := range div.Data() {
				error = divVal - dcD[i]
				if i == 0 {
					min, max = error, error
				}
				if min > error {
					error = min
				}
				if max < error {
					error = max
				}
			}
			fmt.Printf("Min, Max Err = %8.5f, %8.5f\n", min, max)
			return
		}
		checkSolution := func(dfr *DFR2D, Order int) (Fx, Fy, divCheck utils.Matrix) {
			var (
				rt  = dfr.FluxElement
				K   = dfr.K
				Np  = rt.Np
				ccf = float64(Order)
			)
			divCheck = utils.NewMatrix(K, Np)
			Fx, Fy = utils.NewMatrix(K, Np), utils.NewMatrix(K, Np)
			for k := 0; k < K; k++ {
				var (
					xD, yD   = dfr.FluxX.Col(k).Data()[0:Np], dfr.FluxY.Col(k).Data()[0:Np]
					fxD, fyD = Fx.Data(), Fy.Data()
					dcD      = divCheck.Data()
				)
				for n := 0; n < Np; n++ {
					x, y := xD[n], yD[n]
					fxD[n+k*Np], fyD[n+k*Np] = utils.POW(x, Order), utils.POW(y, Order)
					dcD[n+k*Np] = ccf * (utils.POW(x, Order-1) + utils.POW(y, Order-1))
				}
			}
			return
		}
		if true { // Check Divergence for polynomial vector fields of order < N against analytical solution
			N := 3
			dfr := NewDFR2D(N, "test_tris_5.neu")
			rt := dfr.FluxElement
			Dr := rt.Dr[0]
			Ds := rt.Ds[1]
			for cOrder := 1; cOrder <= N; cOrder++ {
				fmt.Printf("Check Order = %d, \n", cOrder)
				Fx, Fy, divCheck := checkSolution(dfr, cOrder)
				Fxp, Fyp := dfr.ProjectFluxOntoRTSpace(Fx, Fy)
				var (
					K = dfr.K
				)
				for k := 0; k < K; k++ {
					var (
						Fxpk, Fypk = Fxp.Row(k).ToMatrix(), Fyp.Row(k).ToMatrix()
						Jdet       = dfr.Jdet.Row(k).Data()[0]
					)
					//fmt.Printf("Jdet[%d]=%v\n", k, Jdet)
					//fmt.Printf("%s\n", Fxpk.Transpose().Print("Fxpk"))
					//fmt.Printf("%s\n", Fypk.Transpose().Print("Fypk"))
					div := Dr.Mul(Fxpk).Add(Ds.Mul(Fypk)).Scale(1. / Jdet)
					var divTotal float64
					for _, divVal := range div.Data() {
						divTotal += divVal
					}
					fmt.Printf("divTotal[%d] = %8.5f\n", k, divTotal)
					//div := Dr.Mul(Fxpk).Add(Ds.Mul(Fypk)).Scale(1. / Jdet)
					minErr, maxErr := errorCheck(dfr, k, div, divCheck)
					//assert.True(t, near(minErr, 0.0, 0.00001))
					//assert.True(t, near(maxErr, 0.0, 0.00001))
					_, _, _, _, _ = div, errorCheck, divCheck, minErr, maxErr
				}
			}
		}
	}
	{ // Test face construction
		/*
			N := 1
			dfr := NewDFR2D(N, "test_tris_5.neu")
			//dfr := NewDFR2D(N, "fstepA001.neu")
			fmt.Printf("%s\n", dfr.EToV.Print("EToV"))
		*/
		//fmt.Printf("%s\n", dfr.EToF.Print("EToF"))
		//fmt.Printf("%s\n", dfr.EToE.Print("EToE"))
		//PlotMesh(dfr.VX, dfr.VY, dfr.EToV, dfr.BCType, dfr.FluxX, dfr.FluxY, true)
		//PlotMesh(dfr.VX, dfr.VY, dfr.EToV, dfr.BCType, dfr.SolutionX, dfr.SolutionY, true)
		//utils.SleepFor(50000)
		//el := dfr.SolutionElement
		//rt := dfr.FluxElement
	}
}
