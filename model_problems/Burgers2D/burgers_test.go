package Burgers2D

import (
	"testing"

	"github.com/notargets/gocfd/DG2D"
)

func TestNewBurgers2D(t *testing.T) {
	P := 2
	b2d := &Burgers2D{}
	mesh := "../../DG2D/test_data/test_tris_9.neu"
	b2d.DFR = DG2D.NewDFR2D(P, true, mesh)
	b2d.DFR.VX.Transpose().Print("VX")
	b2d.DFR.VY.Transpose().Print("VY")
	// cc := chart2d.NewChart2D(1920, 1080, 0, 10, 0, 10)
	// cc.AddTriMesh()
	return
}
