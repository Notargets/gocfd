package model_problems

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG1D"
)

func TestFlux(t *testing.T) {
	/*
		Check face mapping
	*/
	{

		K := 4
		N := 1
		VX, EToV := DG1D.SimpleMesh1D(0, 1, K)
		var c *EulerDFR
		c = &EulerDFR{
			CFL:       1,
			State:     NewFieldState(),
			FinalTime: 20,
			El:        DG1D.NewElements1D(N, VX, EToV),
		}
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
		s.Update(c.Rho, c.RhoU, c.Ener)
		RhoF = c.RhoU.Copy()
		RhoUF = s.Q.Copy().Scale(2.).Add(s.Pres)
		EnerF = c.Ener.Copy().Add(s.Pres).ElMul(s.U)
		fmt.Println(RhoF.Print("RhoF"))
		fmt.Println(RhoUF.Print("RhoUF"))
		fmt.Println(EnerF.Print("EnerF"))

		fRho, fRhoU, fEner = c.RoeFlux(c.Rho, c.RhoU, c.Ener, RhoF, RhoUF, EnerF)
		//fRho, fRhoU, fEner = c.LFFlux(c.Rho, c.RhoU, c.Ener, RhoF, RhoUF, EnerF)

		// Set face flux within global flux
		RhoF.AssignVector(el.VmapM, fRho)
		RhoUF.AssignVector(el.VmapM, fRhoU)
		EnerF.AssignVector(el.VmapM, fEner)
		fmt.Println(RhoF.Print("RhoF"))
		fmt.Println(RhoUF.Print("RhoUF"))
		fmt.Println(EnerF.Print("EnerF"))
	}
}
