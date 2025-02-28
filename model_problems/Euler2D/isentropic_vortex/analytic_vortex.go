package isentropic_vortex

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type IVortex struct {
	Beta, X0, Y0, Gamma float64
	Ufs                 float64
}

func NewIVortex(Beta, X0, Y0, Gamma float64, UfsO ...float64) (iv *IVortex) {
	var (
		Ufs = 1.0
	)
	if len(UfsO) > 0 {
		Ufs = UfsO[0]
	}
	iv = &IVortex{
		Beta:  Beta,
		X0:    X0,
		Y0:    Y0,
		Gamma: Gamma,
		Ufs:   Ufs,
	}
	return
}

func (iv *IVortex) GetState(t, x, y float64) (u, v, rho, p float64) {
	/*
	  // reset static arrays to constant values
	  u=1.0;  v=0.0;

	  // base flow parameters
	  double xo=5.0, yo=0.0, beta=5.0;
	  double fac = 16.0*gamma*pi*pi;

	  xmut = xi - u*ti;   ymvt = yi - v*ti;
	  rsqr = sqr(xmut-xo)+sqr(ymvt-yo);
	  ex1r = exp(1.0-rsqr);

	  // perturbed density
	  u -= beta * ex1r.dm(ymvt-yo)/(2.0*pi);
	  v += beta * ex1r.dm(xmut-xo)/(2.0*pi);

	  tv1  = (1.0-(gm1*SQ(beta)*exp(2.0*(1.0-rsqr))/fac));
	  rho1 = pow(tv1, 1.0/gm1);
	  p1   = pow(rho1, gamma);
	*/
	var (
		oo2pi = 0.5 * (1. / math.Pi)
		Gamma = iv.Gamma
		GM1   = Gamma - 1
		OOGM1 = 1. / GM1
		pi2   = math.Pi * math.Pi
		beta  = iv.Beta
		beta2 = beta * beta
		fac   = 16 * Gamma * pi2
	)
	u, v = iv.Ufs, 0.          // start with freestream values, perturb them later
	xmut, ymvt := x-u*t, y-v*t // vortex center location at time t
	r2 := utils.POW(xmut-iv.X0, 2) + utils.POW(ymvt-iv.Y0, 2)
	ex1r := math.Exp(1 - r2)
	tv1 := 1.0 - (GM1 * beta2 * math.Exp(2.0*(1.0-r2)) / fac)
	u -= beta * ex1r * (ymvt - iv.Y0) * oo2pi
	v += beta * ex1r * (xmut - iv.X0) * oo2pi
	rho = math.Pow(tv1, OOGM1)
	p = math.Pow(rho, Gamma)
	return
}

func (iv *IVortex) GetStateC(t, x, y float64) (Rho, RhoU, RhoV, E float64) {
	var (
		ooGM1 = 1. / (iv.Gamma - 1.)
	)
	u, v, rho, p := iv.GetState(t, x, y)
	q := 0.5 * rho * (u*u + v*v)
	Rho, RhoU, RhoV, E = rho, rho*u, rho*v, p*ooGM1+q
	return
}

func (iv *IVortex) GetDivergence(t, x, y float64) (div [4]float64) {
	// From a matlab calculation of the divergence of the analytic solution
	// Matlab script is in research/test_cases/isentropic_vortex.m
	var (
		beta   = iv.Beta
		gamma  = iv.Gamma
		x0, y0 = iv.X0, iv.Y0
		exp    = math.Exp
		pow    = math.Pow
		pi     = math.Pi
		pi2    = pi * pi
		twopi  = 2. * pi
		beta2  = beta * beta
		GM1    = gamma - 1.
	)
	div[0] = (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1)*(t*2.0-x*2.0+x0*2.0)*(-1.0/2.0))/3.141592653589793 - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*(x-x0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1))/(twopi) + (beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*1.6e+1) + ((beta2*beta)*1.0/(pi2*3.141592653589793)*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*(x-x0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0))/(gamma*3.2e+1)
	div[1] = beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma-1.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0)*(t*4.0-x*4.0+x0*4.0)*(-1.0/1.6e+1) - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0))/(twopi)-(beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*(y-y0))/(twopi))*(x-x0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1))/(twopi) + (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(y-y0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1)*(t*2.0-x*2.0+x0*2.0))/3.141592653589793 + (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(x-x0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1))/(twopi) - (beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*1.6e+1) - ((beta2*beta)*1.0/(pi2*3.141592653589793)*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(x-x0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0))/(gamma*3.2e+1)
	div[2] = (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1)*(-1.0/2.0))/3.141592653589793 + (beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma-1.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0))/1.6e+1 - (beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(x-x0, 2.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1))/4.0 - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(x-x0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1)*(t*2.0-x*2.0+x0*2.0))/(twopi) - (beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(x-x0)*(y-y0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1)*(t*2.0-x*2.0+x0*2.0))/4.0 + ((beta2*beta2)*1.0/(pi2*pi2)*exp(pow(y-y0, 2.0)*-4.0-pow(t-x+x0, 2.0)*4.0+4.0)*(y*4.0-y0*4.0)*pow(x-x0, 2.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0))/(gamma*6.4e+1) + ((beta2*beta)*1.0/(pi2*3.141592653589793)*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(x-x0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*3.2e+1)
	div[3] = ((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1)*((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(x*2.0-x0*2.0))/4.0+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0)*(t*4.0-x*4.0+x0*4.0))/4.0+(beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(y-y0)*(t*2.0-x*2.0+x0*2.0))/3.141592653589793)*(-1.0/2.0)+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma-1.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0)*(t*4.0-x*4.0+x0*4.0))/1.6e+1+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma-1.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*1.6e+1-1.6e+1)+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*3.2e+1)) + (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(x-x0)*(((((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0))/(twopi)-(beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*(y-y0))/(twopi))*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*2.0-(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(x-x0, 2.0))/4.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1))/2.0+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma-1.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0))/1.6e+1+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma-1.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0))/(gamma*1.6e+1-1.6e+1)+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*(pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1-1.0))/(gamma*3.2e+1)))/(twopi) - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0)*(((pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1))/2.0+pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma)+pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma)/GM1)*(t*2.0-x*2.0+x0*2.0))/(twopi) - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*(x-x0)*(((pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+(beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1))/2.0+pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma)+pow(pow((beta2*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*GM1*(-1.0/1.6e+1))/gamma+1.0, 1.0/GM1), gamma)/GM1))/(twopi)
	return
}

func (iv *IVortex) GetFlux(t, x, y float64) (Fx, Fy [4]float64) {
	rho, rhoU, rhoV, E := iv.GetStateC(t, x, y)
	Fx, Fy = FluxCalc(iv.Gamma, rho, rhoU, rhoV, E)
	return
}

func FluxCalc(Gamma, rho, rhoU, rhoV, E float64) (Fx, Fy [4]float64) {
	var (
		GM1 = Gamma - 1.
	)
	u := rhoU / rho
	v := rhoV / rho
	u2 := u*u + v*v
	q := 0.5 * rho * u2
	p := GM1 * (E - q)
	Fx, Fy =
		[4]float64{rhoU, rhoU*u + p, rhoU * v, u * (E + p)},
		[4]float64{rhoV, rhoV * u, rhoV*v + p, v * (E + p)}
	return
}
