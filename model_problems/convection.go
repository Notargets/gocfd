package model_problems

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"

	"github.com/notargets/gocfd/DG1D"
)

type Convection1D struct {
	// Input parameters
	a, CFL, FinalTime float64
	El                *DG1D.Elements1D
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
