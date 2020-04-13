package model_problems

import (
	"fmt"
	"math"
	"sync"
	"time"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
)

type Maxwell1D struct {
	// Input parameters
	CFL, FinalTime float64
	El             *DG1D.Elements1D
	RHSOnce        sync.Once
	E, H           utils.Matrix
	Epsilon, Mu    utils.Matrix
}

func NewMaxwell1D(CFL, FinalTime float64, Elements *DG1D.Elements1D) *Maxwell1D {
	var (
		epsData = utils.ConstArray(Elements.K, 1)
	)
	//Mu1 := utils.NewVectorConstant(Elements.K, 1)
	for i := Elements.K / 2; i < Elements.K; i++ {
		epsData[i] = 2
	}
	Eps1 := utils.NewVector(Elements.K, epsData)
	ones := utils.NewVectorConstant(Elements.Np, 1)
	Epsilon := Eps1.Outer(ones)
	fmt.Printf("Epsilon = \n%v\n", mat.Formatted(Epsilon, mat.Squeeze()))
	return &Maxwell1D{
		CFL:       CFL,
		FinalTime: FinalTime,
		El:        Elements,
	}
}

func (c *Maxwell1D) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
		el        = c.El
		chart     *chart2d.Chart2D
		colorMap  *utils2.ColorMap
		chartName string
	)
	return
	xmin := el.X.Row(1).Subtract(el.X.Row(0)).Apply(math.Abs).Min()
	dt := xmin * c.CFL
	Ns := math.Ceil(c.FinalTime / dt)
	dt = c.FinalTime / Ns
	Nsteps := int(Ns)
	fmt.Printf("Min Dist = %8.6f, dt = %8.6f, Nsteps = %d\n\n", xmin, dt, Nsteps)
	U := el.X.Copy().Apply(math.Sin)
	//fmt.Printf("U = \n%v\n", mat.Formatted(U, mat.Squeeze()))
	resid := utils.NewMatrix(el.Np, el.K)
	if showGraph {
		chart = chart2d.NewChart2D(1024, 768, float32(el.X.Min()), float32(el.X.Max()), -1, 1)
		colorMap = utils2.NewColorMap(-1, 1, 1)
		chartName = "Maxwell1D"
		go chart.Plot()
	}
	var Time, timelocal float64
	for tstep := 0; tstep < Nsteps; tstep++ {
		for INTRK := 0; INTRK < 5; INTRK++ {
			timelocal = Time + dt*utils.RK4c[INTRK]
			RHSU := c.RHS(U, timelocal)
			// resid = rk4a(INTRK) * resid + dt * rhsu;
			resid.Scale(utils.RK4a[INTRK]).Add(RHSU.Scale(dt))
			// u += rk4b(INTRK) * resid;
			U.Add(resid.Copy().Scale(utils.RK4b[INTRK]))
		}
		Time += dt
		if showGraph {
			if len(graphDelay) != 0 {
				time.Sleep(graphDelay[0])
			}
			if err := chart.AddSeries(chartName,
				el.X.Transpose().RawMatrix().Data,
				U.Transpose().RawMatrix().Data,
				chart2d.CrossGlyph, chart2d.Dashed,
				colorMap.GetRGB(0)); err != nil {
				panic("unable to add graph series")
			}
		}
		if tstep%50 == 0 {
			fmt.Printf("Time = %8.4f, max_resid[%d] = %8.4f, umin = %8.4f, umax = %8.4f\n", Time, tstep, resid.Max(), U.Col(0).Min(), U.Col(0).Max())
		}
	}
	fmt.Printf("U = \n%v\n", mat.Formatted(U, mat.Squeeze()))
}

func (c *Maxwell1D) RHS(U utils.Matrix, time float64) (RHSU utils.Matrix) {
	return
}
