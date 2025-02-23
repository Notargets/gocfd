package DG2D

import (
	"testing"
)

func TestPlotEquiTri(t *testing.T) {
	dfr := CreateEquiTriMesh(1, 50)
	_ = dfr
	if testing.Verbose() {
		PlotDFRElements(dfr)
	}
}
