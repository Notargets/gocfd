package Euler2D

import (
	"fmt"
	"image/color"

	"github.com/notargets/gocfd/InputParameters"

	"github.com/notargets/avs/chart2d"
	"github.com/notargets/avs/functions"
	graphics2D "github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/utils"
)

type ChartState struct {
	chart *chart2d.Chart2D
	fs    *functions.FSurface
	gm    *graphics2D.TriMesh
}

func (c *Euler) GetPlotField(Q [4]utils.Matrix, plotField FlowFunction) (field utils.Matrix) {
	var (
		Kmax       = c.dfr.K
		Np         = c.dfr.SolutionElement.Np
		NpFlux     = c.dfr.FluxElement.Np
		fld        utils.Matrix
		skipInterp bool
	)
	switch {
	case plotField <= Energy:
		fld = Q[int(plotField)]
	case plotField < ShockFunction:
		fld = utils.NewMatrix(Np, Kmax)
		fldD := fld.DataP
		for ik := 0; ik < Kmax*Np; ik++ {
			fldD[ik] = c.FSFar.GetFlowFunction(Q, ik, plotField)
		}
	case plotField == ShockFunction:
		var sf *ModeAliasShockFinder
		if c.Limiter != nil {
			sf = c.Limiter.ShockFinder[0]
		} else {
			sf = NewAliasShockFinder(c.dfr.SolutionElement, 2)
		}
		fld = utils.NewMatrix(Np, Kmax)
		fldD := fld.DataP
		for k := 0; k < Kmax; k++ {
			qElement := Q[0].Col(k)
			m := sf.ShockIndicator(qElement.DataP)
			for i := 0; i < Np; i++ {
				ik := k + i*Kmax
				fldD[ik] = m
			}
		}
	case plotField == EpsilonDissipation:
		field = c.Dissipation.GetScalarEpsilonPlotField(c)
		skipInterp = true
	case plotField == EpsilonDissipationC0:
		field = c.Dissipation.GetC0EpsilonPlotField(c)
		skipInterp = true
	case plotField == XGradientDensity || plotField == YGradientDensity ||
		plotField == XGradientXMomentum || plotField == YGradientXMomentum ||
		plotField == XGradientYMomentum || plotField == YGradientYMomentum ||
		plotField == XGradientEnergy || plotField == YGradientEnergy:
		GradX, GradY := utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
		DOFX, DOFY := utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
		var varNum int
		switch plotField {
		case XGradientDensity, YGradientDensity:
			varNum = 0
		case XGradientXMomentum, YGradientXMomentum:
			varNum = 1
		case XGradientYMomentum, YGradientYMomentum:
			varNum = 2
		case XGradientEnergy, YGradientEnergy:
			varNum = 3
		}
		c.GetSolutionGradientUsingRTElement(-1, varNum, Q, GradX, GradY, DOFX, DOFY)
		if plotField < 300 {
			field = GradX
		} else {
			field = GradY
		}
		skipInterp = true
	}
	if !skipInterp {
		field = c.dfr.FluxInterp.Mul(fld)
	}
	return
}

func (c *Euler) PlotQ(pm *InputParameters.PlotMeta, width, height int) {
	var (
		Q         = c.RecombineShardsKBy4(c.Q)
		plotField = FlowFunction(pm.Field)
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
	fmt.Printf(" Plot>%s min,max = %8.5f,%8.5f\n", plotField.String(), oField.Min(), oField.Max())
	c.PlotFS(width, height,
		pm.FieldMinP, pm.FieldMaxP, 0.99*oField.Min(), 1.01*oField.Max(),
		scale, translate, lineType)
	utils.SleepFor(int(delay.Milliseconds()))
	return
}

func (c *Euler) PlotFS(width, height int,
	fminP, fmaxP *float64, fmin, fmax float64,
	scale float64, translate [2]float32, ltO ...chart2d.LineType) {
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
		//c.chart.chart = chart2d.NewChart2D(3840, 2160, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
		c.chart.chart = chart2d.NewChart2D(width, height, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
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
