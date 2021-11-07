package utils

import (
	"image/color"
	"time"

	"github.com/notargets/avs/functions"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"

	graphics2D "github.com/notargets/avs/geometry"
)

type ColorName uint8

const (
	White ColorName = iota
	Blue
	Red
	Green
	Black
)

func GetColor(name ColorName) (c color.RGBA) {
	switch name {
	case White:
		c = color.RGBA{
			R: 255,
			G: 255,
			B: 255,
			A: 0,
		}
	case Blue:
		c = color.RGBA{
			R: 50,
			G: 0,
			B: 255,
			A: 0,
		}
	case Red:
		c = color.RGBA{
			R: 255,
			G: 0,
			B: 50,
			A: 0,
		}
	case Green:
		c = color.RGBA{
			R: 25,
			G: 255,
			B: 25,
			A: 0,
		}
	case Black:
		c = color.RGBA{
			R: 0,
			G: 0,
			B: 0,
			A: 0,
		}
	}
	return
}

func SleepFor(milliseconds int) {
	time.Sleep(time.Duration(milliseconds) * time.Millisecond)
}

func ArraysTo2Vector(r1, r2 []float64, scaleO ...float64) (g [][2]float64) {
	var (
		scale float64 = 1
	)
	g = make([][2]float64, len(r1))
	if len(scaleO) > 0 {
		scale = scaleO[0]
	}
	for i := range r1 {
		g[i][0] = r1[i] * scale
		g[i][1] = r2[i] * scale
	}
	return
}

func ArraysToPoints(r1, r2 []float64) (points []graphics2D.Point) {
	points = make([]graphics2D.Point, len(r1))
	for i := range r1 {
		points[i].X[0] = float32(r1[i])
		points[i].X[1] = float32(r2[i])
	}
	return
}

type LineChart struct {
	Chart    *chart2d.Chart2D
	ColorMap *utils2.ColorMap
}

func NewLineChart(width, height int, xmin, xmax, fmin, fmax float64) (lc *LineChart) {
	lc = &LineChart{
		Chart:    chart2d.NewChart2D(width, height, float32(xmin), float32(xmax), float32(fmin), float32(fmax)),
		ColorMap: utils2.NewColorMap(-1, 1, 1),
	}
	go lc.Chart.Plot()
	return
}

func (lc *LineChart) Plot(graphDelay time.Duration, x, f []float64, lineColor float64, lineName string) {
	/*
		lineColor goes from -1 (red) to 1 (blue)
	*/
	pSeries := func(field []float64, name string, color float32, gl chart2d.GlyphType) {
		if err := lc.Chart.AddSeries(name, x, f,
			gl, chart2d.Solid, lc.ColorMap.GetRGB(color)); err != nil {
			panic("unable to add graph series")
		}
	}
	pSeries(f, lineName, float32(lineColor), chart2d.NoGlyph)
	time.Sleep(graphDelay)
	return
}

type SurfacePlot struct {
	Chart        *chart2d.Chart2D
	ColorMap     *utils2.ColorMap
	GraphicsMesh *graphics2D.TriMesh
}

func NewSurfacePlot(width, height int, xmin, xmax, ymin, ymax float64,
	gm *graphics2D.TriMesh) (sp *SurfacePlot) {
	sp = &SurfacePlot{
		Chart:        chart2d.NewChart2D(width, height, float32(xmin), float32(xmax), float32(ymin), float32(ymax)),
		GraphicsMesh: gm,
	}
	go sp.Chart.Plot()
	return
}

func (sp *SurfacePlot) AddColorMap(fmin, fmax float64) {
	sp.ColorMap = utils2.NewColorMap(float32(fmin), float32(fmax), 1.)
}

func (sp *SurfacePlot) AddFunctionSurface(field []float32) {
	/*
		oField should be [NpFlux, K]
	*/
	var (
		noLine = chart2d.NoLine
		white  = color.RGBA{R: 255, G: 255, B: 255, A: 1}
		black  = color.RGBA{R: 0, G: 0, B: 0, A: 1}
	)
	_, _ = white, black
	fs := functions.NewFSurface(sp.GraphicsMesh, [][]float32{field}, 0)
	if err := sp.Chart.AddFunctionSurface("FSurface", *fs, noLine, white); err != nil {
		panic("unable to add function surface series")
	}
}
