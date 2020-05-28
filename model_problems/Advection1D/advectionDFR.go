package Advection1D

import (
	"fmt"
	"math"
	"sync"
	"time"

	"github.com/notargets/avs/chart2d"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type AdvectionDFR struct {
	// Input parameters
	a, CFL, FinalTime float64
	El                *DG1D.Elements1D
	UFlux, F          utils.Matrix
	RHSOnce, PlotOnce sync.Once
	chart             *chart2d.Chart2D
	colorMap          *utils2.ColorMap
}

func NewAdvectionDFR(a, CFL, FinalTime, XMax float64, N, K int) *AdvectionDFR {
	if XMax == 0 {
		XMax = 2 * math.Pi
	}
	VX, EToV := DG1D.SimpleMesh1D(0, XMax, K)
	return &AdvectionDFR{
		a:         a,
		CFL:       CFL,
		FinalTime: FinalTime,
		El:        DG1D.NewElements1D(N, VX, EToV),
	}
}

func (c *AdvectionDFR) Run(showGraph bool, graphDelay ...time.Duration) {
	var (
		el           = c.El
		logFrequency = 50
	)
	xmin := el.X.Row(1).Subtract(el.X.Row(0)).Apply(math.Abs).Min()
	dt := 0.5 * xmin * (c.CFL / c.a)
	Ns := math.Ceil(c.FinalTime / dt)
	dt = c.FinalTime / Ns
	Nsteps := int(Ns)
	U := el.X.Copy().Apply(math.Sin)
	fmt.Printf("Umin, Umax = %8.5f, %8.5f\n", U.Min(), U.Max())
	resid := utils.NewMatrix(el.Np, el.K)

	var Time, timelocal float64
	for tstep := 0; tstep < Nsteps; tstep++ {
		c.Plot(showGraph, graphDelay, U)
		for INTRK := 0; INTRK < 5; INTRK++ {
			timelocal = Time + dt*utils.RK4c[INTRK]
			RHSU := c.RHS(U, timelocal)
			// resid = rk4a(INTRK) * resid + dt * rhsu;
			resid.Scale(utils.RK4a[INTRK]).Add(RHSU.Scale(dt))
			// u += rk4b(INTRK) * resid;
			U.Add(resid.Copy().Scale(utils.RK4b[INTRK]))
		}
		//c.Plot(showGraph, graphDelay, c.F)
		Time += dt
		if tstep%logFrequency == 0 {
			fmt.Printf("Time = %8.4f, max_resid[%d] = %8.4f, umin = %8.4f, umax = %8.4f\n", Time, tstep, resid.Max(), U.Col(0).Min(), U.Col(0).Max())
		}
	}
}

func (c *AdvectionDFR) RHS(U utils.Matrix, Time float64) (RHSU utils.Matrix) {
	var (
		uin   float64
		alpha = 0.0 // flux splitting parameter, 0 is full upwinding
		el    = c.El
		limit = true
	)
	c.RHSOnce.Do(func() {
		aNX := el.NX.Copy().Scale(c.a)
		aNXabs := aNX.Copy().Apply(math.Abs).Scale(1. - alpha)
		c.UFlux = aNX.Subtract(aNXabs)
	})
	// Global Flux
	uin = -math.Sin(c.a * Time)
	//U.AssignScalar(el.MapI, uin)
	c.F = U.Copy().Scale(c.a)

	// Face fluxes
	// du = (u(vmapM)-u(vmapP)).dm(a*nx-(1.-alpha)*abs(a*nx))/2.;
	duNr := el.Nfp * el.NFaces
	duNc := el.K
	dU := U.Subset(el.VmapM, duNr, duNc).Subtract(U.Subset(el.VmapP, duNr, duNc)).ElMul(c.UFlux).Scale(0.5)

	// Boundaries
	// Inflow boundary
	// du(mapI) = (u(vmapI)-uin).dm(a*nx(mapI)-(1.-alpha)*abs(a*nx(mapI)))/2.;
	dU.Assign(el.MapI, U.Subset(el.VmapI, duNr, duNc).AddScalar(-uin).ElMul(c.UFlux.Subset(el.MapI, duNr, duNc)).Scale(0.5))
	dU.AssignScalar(el.MapO, 0)

	// Add the average flux, Avg(Fl, Fr)
	Fface := dU.Add(c.F.Subset(el.VmapM, duNr, duNc).Add(c.F.Subset(el.VmapP, duNr, duNc)).Scale(0.5))

	/*
			The flux value on each side of the face is different - average the two sides after each is calculated from the element
		Example: K = 5 elements, 2 faces per element, one face point per face (assumed)
		Fface = Fface(Nfp*NFaces, K)
		Note that Fface[face1,element0] ~= Fface[face0,element1], we need to average the two to make them continuous
		Before --> Fface =
		-0.2658  5.8548   3.9147  -3.4354  -6.0378
		 5.8841  3.9328  -3.4535  -6.0672  -0.2962
		After --> Fface =
		-0.2658  5.8548   3.9147  -3.4354  -6.0378
		 5.8841  3.9328  -3.4535  -6.0672  -0.2962
	*/
	for k := 0; k < el.K-1; k++ {
		avg := 0.5 * (Fface.At(1, k) + Fface.M.At(0, k+1))
		Fface.M.Set(1, k, avg)
		Fface.M.Set(0, k+1, avg)
	}

	// Set the global flux values at the face to the numerical flux
	c.F.AssignVector(el.VmapM, Fface)

	if limit {
		el.SlopeLimitN(c.F, 20)
	}

	RHSU = el.Dr.Mul(c.F).ElMul(el.Rx).Scale(-1)
	return
}

func (c *AdvectionDFR) Plot(showGraph bool, graphDelay []time.Duration, U utils.Matrix) {
	var (
		el         = c.El
		pMin, pMax = float32(-1), float32(1)
	)
	if !showGraph {
		return
	}
	c.PlotOnce.Do(func() {
		c.chart = chart2d.NewChart2D(1280, 1024, float32(el.X.Min()), float32(el.X.Max()), pMin, pMax)
		c.colorMap = utils2.NewColorMap(-1, 1, 1)
		go c.chart.Plot()
	})

	if err := c.chart.AddSeries("U", el.X.Transpose().RawMatrix().Data, U.Transpose().RawMatrix().Data,
		chart2d.NoGlyph, chart2d.Solid, c.colorMap.GetRGB(0)); err != nil {
		panic("unable to add graph series")
	}
	if len(graphDelay) != 0 {
		time.Sleep(graphDelay[0])
	}
}
