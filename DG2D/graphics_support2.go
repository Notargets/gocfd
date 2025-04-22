package DG2D

import (
	"fmt"
	"image/color"
	"math"

	"github.com/notargets/avs/assets"
	"github.com/notargets/avs/chart2d"
	"github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/utils"
)

type RenderText struct {
	Color color.RGBA
	Text  string
	Pitch uint32
	X, Y  float32
}

func PlotLinesAndText(lines map[color.RGBA][]float32,
	text []RenderText) {
	var (
		xMin, xMax = float32(math.MaxFloat32), -float32(math.MaxFloat32)
		yMin, yMax = float32(math.MaxFloat32), -float32(math.MaxFloat32)
	)
	for _, line := range lines {
		xMin, xMax, yMin, yMax = getMinMax(line, xMin, xMax, yMin, yMax)
	}
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	for col, line := range lines {
		ch.AddLine(line, col)
	}
	for _, txt := range text {
		tf := assets.NewTextFormatter("NotoSans",
			"Regular", txt.Pitch,
			txt.Color, true, false)
		ch.Printf(tf, txt.X, txt.Y, "%s", txt.Text)
	}
	for {
	}
}

func AddLine(x1, y1, x2, y2 float64, col color.RGBA,
	lines map[color.RGBA][]float32) {
	lines[col] = append(lines[col],
		float32(x1), float32(y1),
		float32(x2), float32(y2),
	)
}

func AddCrossHairs(xy []float32, col color.RGBA,
	lines map[color.RGBA][]float32) {
	var (
		lenXY = len(xy) / 2
		size  = float32(0.02)
	)
	for i := 0; i < lenXY; i++ {
		lines[col] = append(lines[col],
			xy[2*i]-size, xy[2*i+1],
			xy[2*i]+size, xy[2*i+1],
			xy[2*i], xy[2*i+1]-size,
			xy[2*i], xy[2*i+1]+size,
		)
	}
}

func PlotTriMesh(gm geometry.TriMesh) {
	var (
		xMin, xMax = float32(math.MaxFloat32), float32(-math.MaxFloat32)
		yMin, yMax = float32(math.MaxFloat32), float32(-math.MaxFloat32)
	)
	xMin, xMax, yMin, yMax = getMinMax(gm.XY, xMin, xMax, yMin, yMax)
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	ch.AddTriMesh(gm)
	for {
	}
}

func getMinMax(XY []float32, xi, xa, yi, ya float32) (xMin, xMax, yMin, yMax float32) {
	var (
		x, y  float32
		lenXY = len(XY) / 2
	)
	for i := 0; i < lenXY; i++ {
		x, y = XY[i*2+0], XY[i*2+1]
		if i == 0 {
			xMin = xi
			xMax = xa
			yMin = yi
			yMax = ya
		} else {
			if x < xMin {
				xMin = x
			}
			if x > xMax {
				xMax = x
			}
			if y < yMin {
				yMin = y
			}
			if y > yMax {
				yMax = y
			}
		}
	}
	return
}

func getFieldMinMax(field []float64) (fMin, fMax float64) {
	for i, f := range field {
		if i == 0 {
			fMin = f
			fMax = f
		}
		if f < fMin {
			fMin = f
		}
		if f > fMax {
			fMax = f
		}
	}
	return
}
func getFieldMinMax32(field []float32) (fMin, fMax float32) {
	for i, f := range field {
		if i == 0 {
			fMin = f
			fMax = f
		}
		if f < fMin {
			fMin = f
		}
		if f > fMax {
			fMax = f
		}
	}
	return
}

func PlotField(field []float64, gm geometry.TriMesh, FMin, FMax float64,
	xMM ...float64) {
	var xMin, xMax, yMin, yMax float32

	if len(xMM) == 4 {
		xMin, xMax = float32(xMM[0]), float32(xMM[1])
		yMin, yMax = float32(xMM[2]), float32(xMM[3])
	} else {
		xMin, xMax = -1, 1
		yMin, yMax = -1, 1
	}
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	var pField []float32
	var fMin, fMax float32
	pField = make([]float32, len(field))
	for i, f := range field {
		pField[i] = float32(f)
	}
	fMin, fMax = getFieldMinMax32(pField)
	vs := geometry.VertexScalar{
		TMesh:       &gm,
		FieldValues: pField,
	}
	fmt.Printf("fMin: %f, fMax: %f\n", fMin, fMax)
	if FMin == 0 && FMax == 0 {
		FMin = float64(fMin)
		FMax = float64(fMax)
	}
	ch.AddShadedVertexScalar(&vs, float32(FMin), float32(FMax))
	ch.AddTriMesh(gm)
	line := []float32{0, -5, 0, 5, -5, 0, 5, 0}
	ch.AddLine(line, utils2.RED)
	for {
	}
}

func (dfr *DFR2D) AverageGraphFieldVertices(field utils.Matrix) {
	// Replace vertices for plotting
	NpEdge := dfr.FluxElement.NpEdge + 2
	for k := 0; k < dfr.K; k++ {
		var iVm int
		for nEdge := 0; nEdge < 3; nEdge++ {
			iV := nEdge * (NpEdge - 1)
			iVp := iV + 1
			if nEdge == 0 {
				iVm = 3*(NpEdge-1) - 1
			} else {
				iVm = nEdge*(NpEdge-1) - 1
			}
			field.Set(iV, k, 0.5*(field.At(iVp, k)+field.At(iVm, k)))
		}
	}
}
