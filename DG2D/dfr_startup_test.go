package DG2D

import (
	"fmt"
	"math"
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
		errorCheck := func(N int, div, divCheck []float64) (minInt, maxInt, minEdge, maxEdge float64) {
			var (
				Npm    = len(div)
				errors = make([]float64, Npm)
			)
			for i := 0; i < Npm; i++ {
				//var ddr, dds float64
				errors[i] = div[i] - divCheck[i]
			}
			minInt, maxInt = errors[0], errors[0]
			Nint := N * (N + 1) / 2
			minEdge, maxEdge = errors[Nint], errors[Nint]
			for i := 0; i < Nint; i++ {
				errAbs := math.Abs(errors[i])
				if minInt > errAbs {
					minInt = errAbs
				}
				if maxInt < errAbs {
					maxInt = errAbs
				}
			}
			for i := Nint; i < Npm; i++ {
				errAbs := math.Abs(errors[i])
				if minEdge > errAbs {
					minEdge = errAbs
				}
				if maxEdge < errAbs {
					maxEdge = errAbs
				}
			}
			fmt.Printf("Order = %d, ", N)
			fmt.Printf("Min, Max Int Err = %8.5f, %8.5f, Min, Max Edge Err = %8.5f, %8.5f\n", minInt, maxInt, minEdge, maxEdge)
			return
		}
		checkSolution := func(rt *RTElement, Order int) (s1, s2, divCheck []float64) {
			var (
				Np = rt.Np
			)
			s1, s2 = make([]float64, Np), make([]float64, Np)
			divCheck = make([]float64, Np)
			var ss1, ss2 float64
			for i := 0; i < Np; i++ {
				r := rt.R.Data()[i]
				s := rt.S.Data()[i]
				ccf := float64(Order)
				s1[i] = utils.POW(r, Order)
				s2[i] = utils.POW(s, Order)
				ss1, ss2 = ccf*utils.POW(r, Order-1), ccf*utils.POW(s, Order-1)
				divCheck[i] = ss1 + ss2
			}
			return
		}
		if false { // Check Divergence for polynomial vector fields of order < N against analytical solution
			Nend := 2
			for N := 1; N < Nend; N++ {
				dfr := NewDFR2D(N, "test_tris_5.neu")
				rt := dfr.FluxElement
				Dr := rt.Dr[0]
				Ds := rt.Ds[1]
				for cOrder := 0; cOrder <= N; cOrder++ {
					fmt.Printf("Check Order = %d, ", cOrder)
					// TODO: Change checkSolution to output divergence and function in the X,Y space of the Flux Element
					s1, s2, divCheck := checkSolution(rt, cOrder)
					s1p, s2p := dfr.ProjectFluxOntoRTSpace(0, s1, s2)
					s1m, s2m := utils.NewMatrix(rt.Np, 1, s1p), utils.NewMatrix(rt.Np, 1, s2p)
					Jdet := dfr.Jdet.Row(0).Data()[0]
					div := Dr.Mul(s1m).Add(Ds.Mul(s2m)).Scale(1. / Jdet).Data()
					//div := rt.Divergence(s1, s2)
					minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, div, divCheck)
					assert.True(t, near(minerrInt, 0.0, 0.00001))
					assert.True(t, near(maxerrInt, 0.0, 0.00001))
					assert.True(t, near(minerrEdge, 0.0, 0.00001))
					assert.True(t, near(maxerrEdge, 0.0, 0.00001))
					// Check the restricted divergence operator for equivalence on the interior points
					div2 := rt.DivergenceInterior(s1, s2)
					Ninterior := N * (N + 1) / 2
					assert.True(t, nearVec(div[0:Ninterior], div2, 0.00001))
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
