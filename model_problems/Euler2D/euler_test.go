package Euler2D

import (
	"fmt"
	"strconv"
	"testing"
)

func TestEuler(t *testing.T) {
	c := NewEuler(1, 1, 1, "../../DG2D/test_tris_5.neu", FLUX_Average, FREESTREAM)
	fmt.Println(c.Q[0].Print("Rho"))
	fmt.Println(c.Q[1].Print("RhoU"))
	fmt.Println(c.Q[2].Print("RhoV"))
	fmt.Println(c.Q[3].Print("RhoE"))
	for ii := 0; ii < 4; ii++ {
		label := strconv.Itoa(ii)
		fmt.Println(c.Fx[ii].Print("Fx" + "[" + label + "]"))
	}
	for ii := 0; ii < 4; ii++ {
		label := strconv.Itoa(ii)
		fmt.Println(c.Fy[ii].Print("Fy" + "[" + label + "]"))
	}
}
