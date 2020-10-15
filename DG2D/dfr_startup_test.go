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
		assert.Equal(t, EdgeNumber([2]int{1, 0}), uint64(4294967296))
		assert.Equal(t, EdgeNumber([2]int{0, 1}), uint64(4294967296))
		assert.Equal(t, EdgeNumber([2]int{0, 10}), uint64(42949672960))
		assert.Equal(t, EdgeNumber([2]int{100, 0}), uint64(429496729600))
		assert.Equal(t, EdgeNumber([2]int{100, 1}), uint64(429496729601))
	}
	{ // Test face construction
		N := 1
		dfr := NewDFR2D(N, "test_tris_5.neu")
		//dfr := NewDFR2D(N, "fstepA001.neu")
		fmt.Printf("%s\n", dfr.EToV.Print("EToV"))
		//fmt.Printf("%s\n", dfr.EToF.Print("EToF"))
		//fmt.Printf("%s\n", dfr.EToE.Print("EToE"))
		//PlotMesh(dfr.VX, dfr.VY, dfr.EToV, dfr.BCType, dfr.FluxX, dfr.FluxY, true)
		//PlotMesh(dfr.VX, dfr.VY, dfr.EToV, dfr.BCType, dfr.SolutionX, dfr.SolutionY, true)
		//utils.SleepFor(50000)
		//el := dfr.SolutionElement
		//rt := dfr.FluxElement
	}
}
