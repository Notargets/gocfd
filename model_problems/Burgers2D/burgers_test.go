package Burgers2D

import (
	"testing"

	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/InputParameters"
)

func TestNewBurgers2D(t *testing.T) {
	P := 2
	b2d := &Burgers2D{}
	pm := &InputParameters.PlotMeta{PlotMesh: false}
	mesh := "../../DG2D/test_data/test_tris_9.neu"
	b2d.DFR = DG2D.NewDFR2D(P, pm, true, mesh)
	b2d.DFR.VX.Transpose().Print("VX")
	b2d.DFR.VY.Transpose().Print("VY")
	return
}
