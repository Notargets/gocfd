package model_problems

import (
	"fmt"
	"math"
	"sync"
	"time"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/avs/chart2d"

	utils2 "github.com/notargets/avs/utils"
)

type Advection1D struct {
	// Input parameters
	a, CFL, FinalTime float64
	El                *DG1D.Elements1D
	RHSOnce           sync.Once
	UFlux             utils.Matrix
}

func NewAdvection1D(a, CFL, FinalTime float64, N, K int) *Advection1D {
	VX, EToV := DG1D.SimpleMesh1D(0, 2*math.Pi, K)
	return &Advection1D{
		a:         a,
		CFL:       CFL,
		FinalTime: FinalTime,
		El:        DG1D.NewElements1D(N, VX, EToV),
	}
}

func (c *Advection1D) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
		el           = c.El
		chart        *chart2d.Chart2D
		colorMap     *utils2.ColorMap
		chartName    string
		logFrequency = 50
	)
	xmin := el.X.Row(1).Subtract(el.X.Row(0)).Apply(math.Abs).Min()
	dt := 0.5 * xmin * (c.CFL / c.a)
	Ns := math.Ceil(c.FinalTime / dt)
	dt = c.FinalTime / Ns
	Nsteps := int(Ns)
	U := el.X.Copy().Apply(math.Sin)
	resid := utils.NewMatrix(el.Np, el.K)
	if showGraph {
		chart = chart2d.NewChart2D(1024, 768, float32(el.X.Min()), float32(el.X.Max()), -1, 1)
		colorMap = utils2.NewColorMap(-1, 1, 1)
		chartName = "Advect1D"
		go chart.Plot()
	}
	var Time, timelocal float64
	for tstep := 0; tstep < Nsteps; tstep++ {
		if showGraph {
			if err := chart.AddSeries(chartName, el.X.Transpose().RawMatrix().Data, U.Transpose().RawMatrix().Data,
				chart2d.CrossGlyph, chart2d.Dashed, colorMap.GetRGB(0)); err != nil {
				panic("unable to add graph series")
			}
			if len(graphDelay) != 0 {
				time.Sleep(graphDelay[0])
			}
		}
		for INTRK := 0; INTRK < 5; INTRK++ {
			//U = el.SlopeLimitN(U.Copy())
			timelocal = Time + dt*utils.RK4c[INTRK]
			RHSU := c.RHS(U, timelocal)
			// resid = rk4a(INTRK) * resid + dt * rhsu;
			resid.Scale(utils.RK4a[INTRK]).Add(RHSU.Scale(dt))
			// u += rk4b(INTRK) * resid;
			U.Add(resid.Copy().Scale(utils.RK4b[INTRK]))
		}
		Time += dt
		if tstep%logFrequency == 0 {
			fmt.Printf("Time = %8.4f, max_resid[%d] = %8.4f, umin = %8.4f, umax = %8.4f\n", Time, tstep, resid.Max(), U.Col(0).Min(), U.Col(0).Max())
		}
	}
}

func (c *Advection1D) RHS(U utils.Matrix, time float64) (RHSU utils.Matrix) {
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
	dU := U.Subset(el.VmapM, duNr, duNc).Subtract(U.Subset(el.VmapP, duNr, duNc)).ElMul(c.UFlux).Scale(0.5)

	// Boundaries
	// Inflow boundary
	// du(mapI) = (u(vmapI)-uin).dm(a*nx(mapI)-(1.-alpha)*abs(a*nx(mapI)))/2.;
	uin = -math.Sin(c.a * time)
	dU.Assign(el.MapI, U.Subset(el.VmapI, duNr, duNc).AddScalar(-uin).ElMul(c.UFlux.Subset(el.MapI, duNr, duNc)).Scale(0.5))
	dU.AssignScalar(el.MapO, 0)

	// rhsu = -a*rx.dm(Dr*u) + LIFT*(Fscale.dm(du));
	// Important: must change the order from Fscale.dm(du) to du.dm(Fscale) here because the dm overwrites the target
	RHSU = el.Rx.Copy().Scale(-c.a).ElMul(el.Dr.Mul(U)).Add(el.LIFT.Mul(dU.ElMul(el.FScale)))
	return
}
