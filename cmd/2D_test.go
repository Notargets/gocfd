package cmd

import (
	"testing"

	"github.com/magiconair/properties/assert"
)

func TestRun2D(t *testing.T) {
	var (
		err error
	)
	fileInput := []byte(`
Title: Test Case
CFL: 1.
FluxType: Roe
Case: IVortex # Can be "Freestream"
PolynomialOrder: 2
BCs: 
  Inflow:
      37:
         NPR: 4.0
  Outflow:
      22:
         P: 1.5
`)
	var input InputParameters
	if err = input.Parse(fileInput); err != nil {
		panic(err)
	}
	// Check Inflow BC number 37
	assert.Equal(t, input.BCs["Inflow"][37]["NPR"], float64(4))
	assert.Equal(t, input.BCs["Outflow"][22]["P"], float64(1.5))
	input.Print()
}
