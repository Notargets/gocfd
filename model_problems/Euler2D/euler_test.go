package Euler2D

import (
	"fmt"
	"strconv"
	"testing"

	"github.com/notargets/gocfd/utils"
)

func TestEuler(t *testing.T) {
	c := NewEuler(1, 1, 1, "../../DG2D/test_tris_5.neu", FLUX_Average, FREESTREAM)
	{ // Test face flux averaging
		el := c.dfr.FluxElement
		for ii := 0; ii < 4; ii++ {
			for i := 0; i < el.Np; i++ {
				c.Fx[ii].Data()[i] = float64(ii + 1)
				c.Fx[ii].Data()[i+el.Np] = 2 * float64(ii+1)
			}
		}
		PrintFlux(c.Fx)
		c.AverageFlux()
		fmt.Printf("After averaging")
		PrintFlux(c.Fx)
	}
}

func PrintQ(Q [4]utils.Matrix) {
	var (
		label string
	)
	for ii := 0; ii < 4; ii++ {
		switch ii {
		case 0:
			label = "Rho"
		case 1:
			label = "RhoU"
		case 2:
			label = "RhoV"
		case 3:
			label = "RhoE"
		}
		fmt.Println(Q[ii].Print(label))
	}
}
func PrintFlux(F [4]utils.Matrix) {
	for ii := 0; ii < 4; ii++ {
		label := strconv.Itoa(ii)
		fmt.Println(F[ii].Print("F" + "[" + label + "]"))
	}
}
