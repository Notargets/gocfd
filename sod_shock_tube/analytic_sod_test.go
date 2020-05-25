package sod_shock_tube

import (
	"fmt"
	"sync"
	"testing"
	"time"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"
)

var (
	o_test        sync.Once
	chart_test    *chart2d.Chart2D
	colormap_test *utils2.ColorMap
)

func TestSOD(t *testing.T) {
	X, Rho, P, U, E := SOD_calc(0.1)
	fmt.Printf("Rho = %v\n", Rho)
	Plot(0, 1, X, Rho, P, U, E)
	for {
		time.Sleep(10 * time.Second)
	}
}

func Plot(xmin, xmax float32, X, Rho, P, U, E []float64) {
	var (
		fmin, fmax = float32(-1.5), float32(5)
	)
	o_test.Do(func() {
		chart_test = chart2d.NewChart2D(1920, 1280, xmin, xmax, fmin, fmax)
		colormap_test = utils2.NewColorMap(-1, 1, 1)
		go chart_test.Plot()
	})
	pSeries := func(field []float64, name string, color float32, gl chart2d.GlyphType) {
		if err := chart_test.AddSeries(name, X, field,
			gl, chart2d.Solid, colormap_test.GetRGB(color)); err != nil {
			panic("unable to add graph series")
		}
	}
	pSeries(Rho, "Rho", -0.7, chart2d.NoGlyph)
	pSeries(P, "P", -0.2, chart2d.NoGlyph)
	pSeries(U, "U", 0.2, chart2d.NoGlyph)
	pSeries(E, "E", 0.7, chart2d.NoGlyph)
}
