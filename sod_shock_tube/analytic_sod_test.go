package sod_shock_tube

import (
	"fmt"
	"math"
	"sync"
	"testing"
	"time"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"
)

var (
	o_test        sync.Once
	chart_test    *chart2d.Chart2D
	colormap_test *utils2.ColorMap
)

func TestSOD(t *testing.T) {
	// plotInteractive(SOD_calc(0.1))
	X, Rho, P, U, E := SOD_calc(0.1)
	fmt.Printf("X = %v\n", X)
	fmt.Printf("Rho = %v\n", Rho)
	xCheck := []float64{0, 0.3816783943380077, 0.3816784143380077, 0.392369832011622, 0.4030612596852364, 0.4137526873588507, 0.424444115032465, 0.43513554270607935, 0.44582697037969365, 0.45651839805330796, 0.46720982572692227, 0.47790125340053663, 0.48859266107415095, 0.48859267107415094, 0.48859269107415093, 0.5890952206134527, 0.5890952406134528, 0.6638699795418159, 0.663869999541816, 1}
	rhoCheck := []float64{1, 1, 0.9999999295704803, 0.9269348762386385, 0.8582043737798012, 0.7936127019079521, 0.7329700565197678, 0.6760925271668523, 0.6228020040984605, 0.5729260853042116, 0.52629798355681, 0.4827564334547598, 0.44214567178754377, 0.4421456351263143, 0.4421455984650873, 0.4421455984650873, 0.27393934779985196, 0.27393934779985196, 0.125, 0.125}
	assert.True(t, isNear(xCheck, X, 0.001))
	assert.True(t, isNear(rhoCheck, Rho, 0.001))
	_, _, _ = P, U, E
}

func isNear(a, b []float64, tol float64) bool {
	if len(a) != len(b) {
		return false
	}
	for i, val := range a {
		if math.Abs(b[i]-val) > tol {
			return false
		}
	}
	return true
}

func plotInteractive(X, Rho, P, U, E []float64) {
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
	pSeries(Rho, "Rho", -0.7, chart2d.CrossGlyph)
	pSeries(E, "E", 0.7, chart2d.BoxGlyph)
	//pSeries(P, "P", 0.7, chart2d.BoxGlyph)
	//pSeries(U, "U", 0.2, chart2d.NoGlyph)
	//pSeries(E, "E", 0.7, chart2d.NoGlyph)
}
