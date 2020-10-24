package Euler2D

import (
	"fmt"
	"math"
	"strconv"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

func TestEuler(t *testing.T) {
	c := NewEuler(1, 1, 1, "../../DG2D/test_tris_5.neu", FLUX_Average, FREESTREAM)
	{ // Test face flux averaging
		el := c.dfr.FluxElement
		for ii := 0; ii < 4; ii++ {
			for i := 0; i < el.Np; i++ {
				c.Fx[ii].Data()[i] = float64(i + 1)
				c.Fx[ii].Data()[i+el.Np] = float64(i + 1)
			}
		}
		/*
			PrintFlux(c.Fx[0:1])
			c.AverageFlux()
			fmt.Printf("After averaging\n")
			PrintFlux(c.Fx[0:1])
		*/
		c.AverageFlux()
		assert.True(t, nearVec(c.Fx[0].Data(),
			[]float64{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 11, 11, 11, 1, 2, 3, 4, 5, 6, 11, 11, 11, 10, 11, 12, 13, 14, 15},
			0.00001))
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
func PrintFlux(F []utils.Matrix) {
	for ii := 0; ii < len(F); ii++ {
		label := strconv.Itoa(ii)
		fmt.Println(F[ii].Print("F" + "[" + label + "]"))
	}
}

func nearVec(a, b []float64, tol float64) (l bool) {
	for i, val := range a {
		if !near(b[i], val, tol) {
			fmt.Printf("Diff = %v, Left[%d] = %v, Right[%d] = %v\n", math.Abs(val-b[i]), i, val, i, b[i])
			return false
		}
	}
	return true
}

func near(a, b float64, tolI ...float64) (l bool) {
	var (
		tol float64
	)
	if len(tolI) == 0 {
		tol = 1.e-08
	} else {
		tol = tolI[0]
	}
	bound := math.Max(tol, tol*math.Abs(a))
	if math.Abs(a-b) <= bound {
		l = true
	}
	return
}