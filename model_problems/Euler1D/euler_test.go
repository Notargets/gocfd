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
			el = c.El
		)
		assert.Equal(t, utils.Index{0, 1, 2, 3, 4, 5, 6, 7}, el.VmapM)
		assert.Equal(t, utils.Index{0, 4, 5, 6, 1, 2, 3, 7}, el.VmapP)
		assert.Equal(t, el.VmapM, el.VmapMS)
		assert.Equal(t, el.VmapP, el.VmapPS)

		model = Euler_DFR_Roe
		c = NewEuler(1, 20, 1, N, K, model, SOD_TUBE)
		el = c.El
		assert.Equal(t, utils.Index{0, 1, 2, 3, 12, 13, 14, 15}, el.VmapM)
		assert.Equal(t, utils.Index{0, 12, 13, 14, 1, 2, 3, 15}, el.VmapP)
		assert.Equal(t, utils.Index{0, 1, 2, 3, 4, 5, 6, 7}, el.VmapMS)
		assert.Equal(t, utils.Index{0, 4, 5, 6, 1, 2, 3, 7}, el.VmapPS)
	}
	// Galerkin Integration
	{
		K := 4
		N := 1
		model := Galerkin_LF
		//model := Euler_DFR_Roe
		c := NewEuler(1, 20, 1, N, K, model, SOD_TUBE)
		var (
			el                 = c.El
			s                  = c.State
			fRho, fRhoU, fEner utils.Matrix
			RhoF, RhoUF, EnerF utils.Matrix
		)
		c.MapSolutionSubset()
		RhoF, RhoUF, EnerF = s.Update(c.Rho, c.RhoU, c.Ener, c.FluxRanger, c.FluxSubset)
		fRho, fRhoU, fEner = c.RoeFlux(c.Rho, c.RhoU, c.Ener, RhoF, RhoUF, EnerF, el.VmapMS, el.VmapPS, el.VmapM, el.VmapP)
		// Set face flux within global flux
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
		model := Euler_DFR_Roe
		//model := Euler_DFR_Roe
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
		RhoF, RhoUF, EnerF = s.Update(c.Rho, c.RhoU, c.Ener, c.FluxRanger, c.FluxSubset)
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

		fRho, fRhoU, fEner = c.RoeFlux(c.Rho, c.RhoU, c.Ener, RhoF, RhoUF, EnerF, el.VmapMS, el.VmapPS, el.VmapM, el.VmapP)
		fRhoUCheck := utils.NewMatrix(2, el.K, []float64{1, 1, 0.55, 0.1, 1, 0.55, 0.1, 0.1})
		assert.Less(t, fRhoUCheck.Subtract(fRhoU).Apply(math.Abs).Max(), 0.0001)
		// Set face flux within global flux
		RhoF.AssignVector(el.VmapM, fRho)
		RhoUF.AssignVector(el.VmapM, fRhoU)
		EnerF.AssignVector(el.VmapM, fEner)
		fRhoUCheck = utils.NewMatrix(2, el.K, []float64{1, 1, 0.55, 0.1, 1, 0.55, 0.1, 0.1})
		assert.Less(t, fRhoUCheck.Subtract(fRhoU).Apply(math.Abs).Max(), 0.0001)
		rhoufCheck = utils.NewMatrix(el.Np, el.K, []float64{1, 1, 0.55, 0.1, 1, 1, 0.1, 0.1, 1, 1, 0.1, 0.1, 1, 0.55, 0.1, 0.1})
		assert.Less(t, rhoufCheck.Subtract(RhoUF).Apply(math.Abs).Max(), 0.0001)
	}
}
