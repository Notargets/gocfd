package Euler2D

import (
	"fmt"
	"image/color"
	"time"

	"github.com/notargets/avs/chart2d"
	"github.com/notargets/avs/functions"
	graphics2D "github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/utils"
)

type PlotMeta struct {
	Plot                   bool
	Scale                  float64
	TranslateX, TranslateY float64
	Field                  FlowFunction
	FieldMinP, FieldMaxP   *float64 // nil if no forced min, max
	FrameTime              time.Duration
	StepsBeforePlot        int
	LineType               chart2d.LineType
}

type ChartState struct {
	chart *chart2d.Chart2D
	fs    *functions.FSurface
	gm    *graphics2D.TriMesh
}

func (c *Euler) GetPlotField(Q [4]utils.Matrix, plotField FlowFunction) (field utils.Matrix) {
	var (
		Kmax = c.dfr.K
		Np   = c.dfr.SolutionElement.Np
		fld  utils.Matrix
	)
	switch {
	case plotField <= Energy:
		//field = c.dfr.FluxInterpMatrix.Mul(Q[int(plotField)])
		fld = Q[int(plotField)]
	case plotField < ShockFunction:
		fld = utils.NewMatrix(Np, Kmax)
		fldD := fld.DataP
		for ik := 0; ik < Kmax*Np; ik++ {
			fldD[ik] = c.FS.GetFlowFunction(Q, ik, plotField)
		}
		//field = c.dfr.FluxInterpMatrix.Mul(fld)
	case plotField == ShockFunction:
		fld = utils.NewMatrix(Np, Kmax)
		fldD := fld.DataP
		//for ik := 0; ik < Kmax*Np; ik++ {
		for k := 0; k < Kmax; k++ {
			qElement := Q[3].Col(k)
			m := c.ShockFinder.ShockIndicator(qElement.DataP)
			for i := 0; i < Np; i++ {
				ik := k + i*Kmax
				fldD[ik] = m
			}
		}
	}
	field = c.dfr.FluxInterpMatrix.Mul(fld)
	return
}

func (c *Euler) PlotQ(Q [4]utils.Matrix, pm *PlotMeta) {
	var (
		plotField = pm.Field
		delay     = pm.FrameTime
		lineType  = pm.LineType
		scale     = pm.Scale
		translate = [2]float32{float32(pm.TranslateX), float32(pm.TranslateY)}
		oField    = c.GetPlotField(Q, plotField)
		fI        = c.dfr.ConvertScalarToOutputMesh(oField)
	)

	if c.chart.gm == nil {
		c.chart.gm = c.dfr.OutputMesh()
	}
	c.chart.fs = functions.NewFSurface(c.chart.gm, [][]float32{fI}, 0)
	fmt.Printf(" Plot>%s min,max = %8.5f,%8.5f\n", pm.Field.String(), oField.Min(), oField.Max())
	c.PlotFS(pm.FieldMinP, pm.FieldMaxP, 0.99*oField.Min(), 1.01*oField.Max(), scale, translate, lineType)
	utils.SleepFor(int(delay.Milliseconds()))
	return
}

func (c *Euler) PlotFS(fminP, fmaxP *float64, fmin, fmax float64, scale float64, translate [2]float32, ltO ...chart2d.LineType) {
	var (
		fs             = c.chart.fs
		trimesh        = fs.Tris
		lt             = chart2d.NoLine
		specifiedScale = fminP != nil || fmaxP != nil
		autoScale      = !specifiedScale
	)
	if c.chart.chart == nil {
		box := graphics2D.NewBoundingBox(trimesh.GetGeometry())
		box = box.Scale(float32(scale))
		box = box.Translate(translate)
		//c.chart.chart = chart2d.NewChart2D(1900, 1080, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
		c.chart.chart = chart2d.NewChart2D(3840, 2160, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
		go c.chart.chart.Plot()
		if specifiedScale {
			// Scale field min/max to preset values
			switch {
			case fminP != nil && fmaxP != nil:
				fmin, fmax = *fminP, *fmaxP
			case fminP != nil:
				fmin = *fminP
			case fmaxP != nil:
				fmax = *fmaxP
			}
			colorMap := utils2.NewColorMap(float32(fmin), float32(fmax), 1.)
			c.chart.chart.AddColorMap(colorMap)
		}
	}
	if autoScale {
		// Autoscale field min/max every time
		colorMap := utils2.NewColorMap(float32(fmin), float32(fmax), 1.)
		c.chart.chart.AddColorMap(colorMap)
	}
	white := color.RGBA{R: 255, G: 255, B: 255, A: 1}
	black := color.RGBA{R: 0, G: 0, B: 0, A: 1}
	_, _ = white, black
	if len(ltO) != 0 {
		lt = ltO[0]
	}
	if err := c.chart.chart.AddFunctionSurface("FSurface", *fs, lt, white); err != nil {
		panic("unable to add function surface series")
	}
}
