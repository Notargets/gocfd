package sod_shock_tube

import (
	"fmt"
	"math"
	"sync"
	"testing"

	"github.com/stretchr/testify/assert"
)

var (
	o_test sync.Once
)

func TestSOD(t *testing.T) {
	// plotInteractive(SOD_calc(0.1))
	sod := NewSOD(0.1)
	X, Rho, _, _, _ := sod.Get()
	fmt.Printf("X = %v\n", X)
	fmt.Printf("Rho = %v\n", Rho)
	xCheck := []float64{0, 0.3815784043380077, 0.3817784043380077, 0.39280783577858336, 0.40393726721915907, 0.4150666986597348, 0.42619613010031043, 0.4373255615408861, 0.4484549929814618, 0.4595844244220375, 0.47071385586261316, 0.4818432873031888, 0.49277271874376455, 0.49287271874376454, 0.4930727187437645, 0.5926452620047974, 0.5928452620047974, 0.675115573202932, 0.675315573202932, 1}
	rhoCheck := []float64{1, 1, 0.9992959031724784, 0.9240353444481086, 0.852758969991083, 0.7859504402212434, 0.7233963393812908, 0.6648901587403833, 0.6102321829702019, 0.5592293765210307, 0.5116952699978237, 0.467449846536279, 0.4270320564069276, 0.42667562327066666, 0.4263194281781805, 0.4263194281781805, 0.26557371170513905, 0.26557371170513905, 0.125, 0.125}
	assert.True(t, isNear(xCheck, X, 0.001))
	assert.True(t, isNear(rhoCheck, Rho, 0.001))
	assert.True(t, math.Abs(sod.x4-0.6752) < 0.0001)
	sod = NewSOD(0.2)
	assert.True(t, math.Abs(sod.x4-0.8504) < 0.0001)
}

func isNear(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i, val := range a {
		if math.Abs(b[i]-val) > tol {
			return false
		}
	}
	return true
}
