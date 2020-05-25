package model_problems

import (
	"fmt"
	"testing"

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
		c.State.Gamma = 1.4
		c.In = NewStateP(c.State.Gamma, 1, 0, 1)
		c.Out = NewStateP(c.State.Gamma, 0.125, 0, 0.1)
		c.InitializeSOD()
		fmt.Println(c.Rho.Print("Rho"))
	}
}
