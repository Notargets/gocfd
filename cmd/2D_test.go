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
InitType: Freestream # Can be IVortex or Freestream
FluxType: Roe
PolynomialOrder: 2
FinalTime: 4.
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
	assert.Equal(t, input.BCs["Inflow"][37]["NPR"], 4.)
	// Check Outflow BC number 22
	assert.Equal(t, input.BCs["Outflow"][22]["P"], 1.5)
	input.Print()
	assert.Equal(t, input.FinalTime, 4.)
}
