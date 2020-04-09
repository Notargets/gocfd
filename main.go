package main

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

const (
	K      = 10
	N      = 8
	NFaces = 2
	Nfp    = 1
)

func main() {
	X := DG1D.Startup1D(K, N, NFaces, Nfp)
	U := X.Copy().Apply(math.Sin)
	_ = U
	//fmt.Printf("U = \n%v\n", mat.Formatted(U, mat.Squeeze()))
	run(X)
}

func run(X utils.Matrix) {
	var (
		a         = 2 * math.Pi
		FinalTime = 5.
		CFL       = 0.75
	)
	_, _, _ = a, FinalTime, CFL
	xmin := X.Row(1).Subtract(X.Row(0)).Apply(math.Abs).Min()
	dt := 0.5 * xmin * (CFL / a)
	Nsteps := math.Ceil(FinalTime / dt)
	fmt.Printf("Min Dist = %8.6f, dt = %8.6f, Nsteps = %5.2f\n\n", xmin, dt, Nsteps)
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
