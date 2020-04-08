package main

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
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
	fmt.Printf("U = \n%v\n", mat.Formatted(U, mat.Squeeze()))
	run(X)
}

func run(X utils.Matrix) {
	var (
		a         = 2 * math.Pi
		FinalTime = 5
	)
	_, _ = a, FinalTime
	_, nc := X.Dims()
	xmin := X.Slice(0, 1, 0, nc)
	ymin := X.Slice(1, 2, 0, nc)
	fmt.Printf("xmin = \n%v\n", mat.Formatted(xmin, mat.Squeeze()))
	fmt.Printf("ymin = \n%v\n", mat.Formatted(ymin, mat.Squeeze()))
	fmt.Printf("X = \n%v\n", mat.Formatted(X, mat.Squeeze()))
	/*
	   double a; // advection speed
	   double xmin, CFL;
	   double timelocal;

	   a = 2.*pi;

	   InitRun();

	   // compute time step size
	   xmin = (abs(x(1,All)-x(2,All))).min_val();
	   CFL=0.75;
	   dt   = .5 * (CFL/(2*pi)*xmin);
	   Nsteps = ceil(FinalTime/dt);
	   dt = FinalTime/Nsteps;

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
