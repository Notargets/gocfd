package model_problems

import (
	"fmt"
	"math"
	"sync"
	"time"

	"gonum.org/v1/gonum/mat"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/avs/chart2d"

	utils2 "github.com/notargets/avs/utils"
)

type Convection1D struct {
	// Input parameters
	a, CFL, FinalTime float64
	El                *DG1D.Elements1D
	RHSOnce           sync.Once
	UFlux             utils.Matrix
}

func NewConvection(a, CFL, FinalTime float64, Elements *DG1D.Elements1D) *Convection1D {
	return &Convection1D{
		a:         a,
		CFL:       CFL,
		FinalTime: FinalTime,
		El:        Elements,
	}
}

func (c *Convection1D) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
		el        = c.El
		chart     *chart2d.Chart2D
		colorMap  *utils2.ColorMap
		chartName string
	)
	xmin := el.X.Row(1).Subtract(el.X.Row(0)).Apply(math.Abs).Min()
	dt := 0.5 * xmin * (c.CFL / c.a)
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
		chartName = "Advect1D"
		if err := chart.AddSeries(chartName,
			ToFloat32Slice(el.X.Transpose().RawMatrix().Data),
			ToFloat32Slice(U.Transpose().RawMatrix().Data),
			chart2d.NoGlyph, chart2d.Solid,
			colorMap.GetRGB(0)); err != nil {
			panic("unable to add graph series")
		}
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
				ToFloat32Slice(el.X.Transpose().RawMatrix().Data),
				ToFloat32Slice(U.Transpose().RawMatrix().Data),
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

func (c *Convection1D) RHS(U utils.Matrix, time float64) (RHSU utils.Matrix) {
	var (
		uin   float64
		alpha = 0.0 // flux splitting parameter, 0 is full upwinding
		el    = c.El
	)
	c.RHSOnce.Do(func() {
		aNX := el.NX.Copy().Scale(c.a)
		aNXabs := aNX.Copy().Apply(math.Abs).Scale(1. - alpha)
		c.UFlux = aNX.Subtract(aNXabs)
	})
	// Face fluxes
	// du = (u(vmapM)-u(vmapP)).dm(a*nx-(1.-alpha)*abs(a*nx))/2.;
	duNr := el.Nfp * el.NFaces
	duNc := el.K
	dU := U.Subset(el.VmapM, duNr, duNc).Subtract(U.Subset(el.VmapP, duNr, duNc)).ElementMultiply(c.UFlux).Scale(0.5)

	// Boundaries
	// Inflow boundary
	// du(mapI) = (u(vmapI)-uin).dm(a*nx(mapI)-(1.-alpha)*abs(a*nx(mapI)))/2.;
	uin = -math.Sin(c.a * time)
	dU.Assign(el.MapI, U.Subset(el.VmapI, duNr, duNc).AddScalar(-uin).ElementMultiply(c.UFlux.Subset(el.MapI, duNr, duNc)).Scale(0.5))
	dU.AssignScalar(el.MapO, 0)

	// rhsu = -a*rx.dm(Dr*u) + LIFT*(Fscale.dm(du));
	// Important: must change the order from Fscale.dm(du) to du.dm(Fscale) here because the dm overwrites the target
	RHSU = el.Rx.Copy().Scale(-c.a).ElementMultiply(el.Dr.Mul(U)).Add(el.LIFT.Mul(dU.ElementMultiply(el.FScale)))
	return
}

func ToFloat32Slice(A []float64) (R []float32) {
	R = make([]float32, len(A))
	for i, val := range A {
		R[i] = float32(val)
	}
	return R
}
