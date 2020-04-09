package model_problems

import (
	"fmt"
	"math"
	"os"
	"sync"

	"github.com/notargets/gocfd/utils"

	"gonum.org/v1/gonum/mat"

	"github.com/notargets/gocfd/DG1D"
)

var (
	rk4a = []float64{
		0.0,
		-567301805773.0 / 1357537059087.0,
		-2404267990393.0 / 2016746695238.0,
		-3550918686646.0 / 2091501179385.0,
		-1275806237668.0 / 842570457699.0,
	}
	rk4b = []float64{
		1432997174477.0 / 9575080441755.0,
		5161836677717.0 / 13612068292357.0,
		1720146321549.0 / 2090206949498.0,
		3134564353537.0 / 4481467310338.0,
		2277821191437.0 / 14882151754819.0,
	}
	rk4c = []float64{
		0.0,
		1432997174477.0 / 9575080441755.0,
		2526269341429.0 / 6820363962896.0,
		2006345519317.0 / 3224310063776.0,
		2802321613138.0 / 2924317926251.0,
	}
)

type Convection1D struct {
	// Input parameters
	a, CFL, FinalTime float64
	El                *DG1D.Elements1D
	RHSOnce           sync.Once
	Uflux             utils.Matrix
}

func NewConvection(a, CFL, FinalTime float64, Elements *DG1D.Elements1D) *Convection1D {
	return &Convection1D{
		a:         a,
		CFL:       CFL,
		FinalTime: FinalTime,
		El:        Elements,
	}
}

func (c *Convection1D) Run() {
	var (
		el = c.El
	)
	xmin := el.X.Row(1).Subtract(el.X.Row(0)).Apply(math.Abs).Min()
	dt := 0.5 * xmin * (c.CFL / c.a)
	Ns := math.Ceil(c.FinalTime / dt)
	dt = c.FinalTime / Ns
	Nsteps := int(Ns)
	fmt.Printf("Min Dist = %8.6f, dt = %8.6f, Nsteps = %d\n\n", xmin, dt, Nsteps)
	U := el.X.Copy().Apply(math.Sin)
	fmt.Printf("U = \n%v\n", mat.Formatted(U, mat.Squeeze()))
	resid := utils.NewMatrix(el.Np, el.K)
	var time, timelocal float64
	for tstep := 0; tstep < Nsteps; tstep++ {
		for INTRK := 0; INTRK < 5; INTRK++ {
			timelocal = time + dt*rk4c[INTRK]
			RHSU := c.RHS(U, timelocal)
			resid.Scale(rk4a[INTRK]).Add(RHSU.Scale(dt))
			U.Add(resid.Copy().Scale(rk4b[INTRK]))
		}
		time += dt
		fmt.Printf("max_resid[%d] = %5.2f, time = %8.4f, dt = %8.4f\n", tstep, resid.Max(), time, dt)
		if tstep == 10 {
			fmt.Printf("U = \n%v\n", mat.Formatted(U, mat.Squeeze()))
			os.Exit(1)
		}
	}
	/*
	   // outer time step loop
	   resid = zeros(Np,K); // Runge-Kutta residual storage
	   time = 0;
	   for (tstep=1; tstep<=Nsteps; tstep++) {
	       for (int INTRK=1; INTRK<=5; INTRK++) {
	           timelocal = time + rk4c(INTRK) * dt;
	           this->RHS(u, timelocal, a);
	           resid = rk4a(INTRK) * resid + dt * rhsu;
	           u += rk4b(INTRK) * resid;
	       }
	       time = time+dt;
	       //umLOG(1, "max_resid[%d] = %g, time = %g, dt = %g\n", tstep, resid.max_val(), time, dt);
	       this->Report(false);
	   }
	   Summary();
	   u.print(stdout, "Solution U");
	*/
}

func (c *Convection1D) RHS(U utils.Matrix, time float64) (RHSU utils.Matrix) {
	var (
		uin   float64
		alpha = 1.0 // flux splitting parameter, 0 is full upwinding
		el    = c.El
	)
	c.RHSOnce.Do(func() {
		aNX := el.NX.Copy().Scale(0.5 * c.a)
		aNXabs := aNX.Copy().Apply(math.Abs).Scale(1. - alpha)
		c.Uflux = aNX.Subtract(aNXabs)
	})
	// Face fluxes
	// du = (u(vmapM)-u(vmapP)).dm(a*nx-(1.-alpha)*abs(a*nx))/2.;
	duNr := el.Nfp * el.NFaces
	duNc := el.K
	dU := U.Subset(el.VmapM, duNr, duNc).Subtract(U.Subset(el.VmapP, duNr, duNc)).ElementMultiply(c.Uflux)

	// Boundaries
	// Inflow boundary
	uin = -math.Sin(c.a * time)
	dU.Assign(el.MapI, U.Subset(el.VmapI, duNr, duNc)).AddScalar(-uin).ElementMultiply(c.Uflux.Subset(el.MapI, duNr, duNc))
	dU.AssignScalar(el.MapO, 0)

	RHSU = el.Rx.Copy().Scale(-c.a).ElementMultiply(el.Dr.Mul(U)).Add(el.LIFT.Mul(el.FScale.ElementMultiply(dU)))
	return
}
