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
}
