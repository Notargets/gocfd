package Euler2D

import (
	"fmt"
	"testing"
)

func TestEuler(t *testing.T) {
	c := NewEuler(1, 1, 1, "../../DG2D/test_tris_5.neu", FLUX_Average, FREESTREAM)
	fmt.Println(c.Rho.Print("Rho"))
	fmt.Println(c.RhoU.Print("RhoU"))
	fmt.Println(c.RhoV.Print("RhoV"))
	fmt.Println(c.RhoE.Print("RhoE"))
	fmt.Println(c.F1x.Print("F1x"))
	fmt.Println(c.F2x.Print("F2x"))
	fmt.Println(c.F3x.Print("F3x"))
	fmt.Println(c.F4x.Print("F4x"))
	fmt.Println(c.F1y.Print("F1y"))
	fmt.Println(c.F2y.Print("F2y"))
	fmt.Println(c.F3y.Print("F3y"))
	fmt.Println(c.F4y.Print("F4y"))
}
