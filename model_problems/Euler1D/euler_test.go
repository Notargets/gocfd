package Euler1D

import (
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

func TestFlux(t *testing.T) {
	/*
		Check face mapping
	*/
	// Solution points to flux points mapping
	{
		K := 4
		N := 1
		model := Galerkin_LF
		c := NewEuler(1, 20, 1, N, K, model, SOD_TUBE)
		var (
			el  = c.El
			elS = c.El_S
		)
		assert.Equal(t, utils.Index{0, 1, 2, 3, 4, 5, 6, 7}, el.VmapM)
		assert.Equal(t, utils.Index{0, 4, 5, 6, 1, 2, 3, 7}, el.VmapP)
		assert.Equal(t, el.VmapM, elS.VmapM)
		assert.Equal(t, el.VmapP, elS.VmapP)

		model = DFR_Roe
		c = NewEuler(1, 20, 1, N, K, model, SOD_TUBE)
		el = c.El
		assert.Equal(t, utils.Index{0, 1, 2, 3, 12, 13, 14, 15}, el.VmapM)
		assert.Equal(t, utils.Index{0, 12, 13, 14, 1, 2, 3, 15}, el.VmapP)
		assert.Equal(t, utils.Index{0, 1, 2, 3, 4, 5, 6, 7}, elS.VmapM)
		assert.Equal(t, utils.Index{0, 4, 5, 6, 1, 2, 3, 7}, elS.VmapP)
	}
	// Galerkin Integration
	{
		K := 4
		N := 1
		model := Galerkin_LF
		// model := DFR_Roe
		c := NewEuler(1, 20, 1, N, K, model, SOD_TUBE)
		var (
			el                 = c.El
			s                  = c.State
			fRho, fRhoU, fEner utils.Matrix
			RhoF, RhoUF, EnerF utils.Matrix
		)
		c.MapSolutionSubset()
		_, _, _, RhoF, RhoUF, EnerF = s.Update(c.Rho, c.RhoU, c.Ener, c)
		fRho, fRhoU, fEner = c.RoeFlux(c.Rho, c.RhoU, c.Ener, RhoF, RhoUF, EnerF, el.VmapM, el.VmapP)
		// SetScalar face flux within global flux
		RhoF.AssignVector(el.VmapM, fRho)
		RhoUF.AssignVector(el.VmapM, fRhoU)
		EnerF.AssignVector(el.VmapM, fEner)
		rhofCheck := utils.NewMatrix(el.Np, el.K, []float64{0, 0, 0.5216, 0, 0, 0.5216, 0, 0})
		assert.Less(t, rhofCheck.Subtract(RhoF).Apply(math.Abs).Max(), 0.0001)
		enerfCheck := utils.NewMatrix(el.Np, el.K, []float64{0, 0, 4.0979, 0, 0, 4.0979, 0, 0})
		assert.Less(t, enerfCheck.Subtract(EnerF).Apply(math.Abs).Max(), 0.0001)
	}
	// DFR Integration
	{
		K := 4
		N := 1
		model := DFR_Roe
		// model := DFR_Roe
		c := NewEuler(1, 20, 1, N, K, model, SOD_TUBE)
		var (
			el                 = c.El
			s                  = c.State
			fRho, fRhoU, fEner utils.Matrix
			RhoF, RhoUF, EnerF utils.Matrix
		)
		/*
			K=4, N=1
				[ 0  1  2  3 ]
				[ 4  5  6  7 ]
		*/
		c.MapSolutionSubset()
		var RhoFull, RhoUFull, EnerFull utils.Matrix
		RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF = s.Update(c.Rho, c.RhoU, c.Ener, c)
		c.CopyBoundary(RhoF)
		c.CopyBoundary(RhoUF)
		c.CopyBoundary(EnerF)
		/*
			⎢1.0000  1.0000  0.1000  0.1000⎥
			⎢1.0000  1.0000  0.1000  0.1000⎥
			⎢1.0000  1.0000  0.1000  0.1000⎥
			⎢1.0000  1.0000  0.1000  0.1000⎥
		*/
		rhoufCheck := utils.NewMatrix(el.Np, el.K, []float64{1, 1, 0.1, 0.1, 1, 1, 0.1, 0.1, 1, 1, 0.1, 0.1, 1, 1, 0.1, 0.1})
		assert.Less(t, rhoufCheck.Subtract(RhoUF).Apply(math.Abs).Max(), 0.0001)

		fRho, fRhoU, fEner = c.RoeFlux(RhoFull, RhoUFull, EnerFull, RhoF, RhoUF, EnerF, el.VmapM, el.VmapP)
		fRhoUCheck := utils.NewMatrix(2, el.K, []float64{1, 1, 0.55, 0.1, 1, 0.55, 0.1, 0.1})
		assert.Less(t, fRhoUCheck.Subtract(fRhoU).Apply(math.Abs).Max(), 0.0001)
		// SetScalar face flux within global flux
		RhoF.AssignVector(el.VmapM, fRho)
		RhoUF.AssignVector(el.VmapM, fRhoU)
		EnerF.AssignVector(el.VmapM, fEner)
		fRhoUCheck = utils.NewMatrix(2, el.K, []float64{1, 1, 0.55, 0.1, 1, 0.55, 0.1, 0.1})
		assert.Less(t, fRhoUCheck.Subtract(fRhoU).Apply(math.Abs).Max(), 0.0001)
		rhoufCheck = utils.NewMatrix(el.Np, el.K, []float64{1, 1, 0.55, 0.1, 1, 1, 0.1, 0.1, 1, 1, 0.1, 0.1, 1, 0.55, 0.1, 0.1})
		assert.Less(t, rhoufCheck.Subtract(RhoUF).Apply(math.Abs).Max(), 0.0001)
	}
}
