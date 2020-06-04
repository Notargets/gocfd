package Euler1D

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

func TestFlux(t *testing.T) {
	/*
		Check face mapping
	*/
	{

		K := 4
		N := 1
		var c *Euler
		model := Galerkin_LF
		//model := Euler_DFR_Roe
		c = NewEuler(1, 20, 1, N, K, model, FREESTREAM)
		var (
			el                 = c.El
			s                  = c.State
			fRho, fRhoU, fEner utils.Matrix
			RhoF, RhoUF, EnerF utils.Matrix
		)
		c.State.Gamma = 1.4
		c.In = NewStateP(c.State.Gamma, 1, 0, 1)
		c.Out = NewStateP(c.State.Gamma, 0.125, 0, 0.1)
		c.InitializeSOD()
		c.MapSolutionSubset()
		RhoF, RhoUF, EnerF = s.Update(c.Rho, c.RhoU, c.Ener, c.FluxRanger, c.FluxSubset)
		fmt.Println(RhoF.Print("RhoF"))
		fmt.Println(RhoUF.Print("RhoUF"))
		fmt.Println(EnerF.Print("EnerF"))

		fRho, fRhoU, fEner = c.RoeFlux(c.Rho, c.RhoU, c.Ener, RhoF, RhoUF, EnerF, el.VmapM, el.VmapP)
		//fRho, fRhoU, fEner = c.LaxFlux(c.Rho, c.RhoU, c.Ener, RhoF, RhoUF, EnerF)

		// Set face flux within global flux
		RhoF.AssignVector(el.VmapM, fRho)
		RhoUF.AssignVector(el.VmapM, fRhoU)
		EnerF.AssignVector(el.VmapM, fEner)
		fmt.Println(RhoF.Print("RhoF"))
		fmt.Println(RhoUF.Print("RhoUF"))
		fmt.Println(EnerF.Print("EnerF"))
		rhofCheck := utils.NewMatrix(el.Np, el.K, []float64{0, 0, 0.5216, 0, 0, 0.5216, 0, 0})
		assert.Less(t, rhofCheck.Subtract(RhoF).Apply(math.Abs).Max(), 0.0001)
		enerfCheck := utils.NewMatrix(el.Np, el.K, []float64{0, 0, 4.0979, 0, 0, 4.0979, 0, 0})
		assert.Less(t, enerfCheck.Subtract(EnerF).Apply(math.Abs).Max(), 0.0001)
	}
}
