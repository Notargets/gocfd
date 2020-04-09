package model_problems

import (
	"fmt"
	"math"
	"sync"

	"github.com/notargets/gocfd/utils"

	"gonum.org/v1/gonum/mat"

	"github.com/notargets/gocfd/DG1D"
)

type Convection1D struct {
	// Input parameters
	a, CFL, FinalTime float64
	El                *DG1D.Elements1D
	RHSOnce           sync.Once
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
	timelocal := 0.00
	c.RHS(U, timelocal)
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

func (c *Convection1D) RHS(U utils.Matrix, time float64) (Resid utils.Matrix) {
	var (
		uin   float64
		alpha = 1. // flux splitting parameter, 0 is full upwinding
		el    = c.El
		Uflux utils.Matrix
	)
	c.RHSOnce.Do(func() {
		aNX := el.NX.Copy().Scale(0.5 * c.a)
		aNXabs := aNX.Copy().Apply(math.Abs).Scale(1. - alpha)
		Uflux = aNX.Subtract(aNXabs)
	})
	// Face fluxes
	dU := U.Subset(el.VmapM).Subtract(U.Subset(el.VmapP)).ElementMultiply(Uflux)
	// Boundaries
	// Inflow boundary
	uin = -math.Sin(c.a * time)
	dU.SubAssign(el.MapI, U.Copy().Subset(el.VmapI).AddScalar(-uin).ElementMultiply(Uflux.Subset(el.MapI)))
	dU.SubAssignScalar(el.MapO, 0)
	fmt.Printf("dU = \n%v\n", mat.Formatted(dU, mat.Squeeze()))
	/*
	   double alpha = 1.;
	   double uin;
	   DMat du = zeros(Nfp*Nfaces, K);

	   // Face fluxes
	   du = (u(vmapM)-u(vmapP)).dm(a*nx-(1.-alpha)*abs(a*nx))/2.;

	   // Boundaries
	   // Inflow boundary
	   uin = -sin(a*time);
	   du(mapI) = (u(vmapI)-uin).dm(a*nx(mapI)-(1.-alpha)*abs(a*nx(mapI)))/2.;
	   du(mapO) = 0.;

	   rhsu = -a*rx.dm(Dr*u) + LIFT*(Fscale.dm(du));
	*/
	return
}
