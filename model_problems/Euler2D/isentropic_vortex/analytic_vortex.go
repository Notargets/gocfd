package isentropic_vortex

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type IVortex struct {
	Beta, X0, Y0, Gamma float64
}

func NewIVortex(Beta, X0, Y0, Gamma float64) (iv *IVortex) {
	iv = &IVortex{
		Beta:  Beta,
		X0:    X0,
		Y0:    Y0,
		Gamma: Gamma,
	}
	return
}

func (iv *IVortex) GetState(t, x, y float64) (u, v, rho, p float64) {
	var (
		oo2pi  = 0.5 * (1. / math.Pi)
		r2     = utils.POW(x-t-iv.X0, 2) + utils.POW(y-iv.Y0, 2)
		omr2   = 1. - r2
		bt     = oo2pi * iv.Beta * math.Exp(omr2)
		bt2    = bt * bt
		Gamma  = iv.Gamma
		GM1    = Gamma - 1
		GM1OGM = GM1 / Gamma
	)
	u = 1. - bt*(y-iv.Y0)
	v = bt * (x - iv.X0)
	rho = math.Pow(1.-bt2*GM1OGM/4., 1./GM1)
	p = math.Pow(rho, Gamma)
	return
}

func (iv *IVortex) GetStateC(t, x, y float64) (Rho, RhoU, RhoV, RhoE float64) {
	var (
		ooGM1 = 1. / (iv.Gamma - 1.)
	)
	u, v, rho, p := iv.GetState(t, x, y)
	q := 0.5 * rho * (u*u + v*v)
	Rho, RhoU, RhoV = rho, rho*u, rho*v
	RhoE = p*ooGM1 + q
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
	div[0] = (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1))*(t*2.0-x*2.0+x0*2.0)*(-1.0/2.0))/pi - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*(x-x0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/(twopi) + ((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*1.6e+1) + ((beta2*beta)*1.0/(pi2*pi)*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*(x-x0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0))/(gamma*3.2e+1)
	div[1] = (beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), GM1)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0)*(t*4.0-x*4.0+x0*4.0)*(-1.0/1.6e+1) - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0))/(twopi)-(beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*(y-y0))/(twopi))*(x-x0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/(twopi) + (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(y-y0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1))*(t*2.0-x*2.0+x0*2.0))/pi + (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(x-x0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/(twopi) - ((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*1.6e+1) - ((beta2*beta)*1.0/(pi2*pi)*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(x-x0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0))/(gamma*3.2e+1)
	div[2] = (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1))*(-1.0/2.0))/pi + ((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), GM1)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0))/1.6e+1 - ((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(x-x0, 2.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/4.0 - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(x-x0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1))*(t*2.0-x*2.0+x0*2.0))/(twopi) - ((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(x-x0)*(y-y0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1))*(t*2.0-x*2.0+x0*2.0))/4.0 + ((beta2*beta2)*1.0/(pi2*pi2)*exp(pow(y-y0, 2.0)*-4.0-pow(t-x+x0, 2.0)*4.0+4.0)*(y*4.0-y0*4.0)*pow(x-x0, 2.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0))/(gamma*6.4e+1) + ((beta2*beta)*1.0/(pi2*pi)*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(x-x0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*3.2e+1)
	div[3] = ((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, -1.0/(GM1))*(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1))*(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(x*2.0-x0*2.0))/4.0+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0)*(t*4.0-x*4.0+x0*4.0))/4.0+(beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*(y-y0)*(t*2.0-x*2.0+x0*2.0))/pi)*(-1.0/2.0)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), GM1)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*1.6e+1-1.6e+1)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*3.2e+1))+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), GM1)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0)*(t*4.0-x*4.0+x0*4.0))/1.6e+1-((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(((pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/2.0+pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), gamma)/(GM1))*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, -1.0/(GM1)-1.0)*(t*4.0-x*4.0+x0*4.0))/(gamma*1.6e+1)) + (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(x-x0)*(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, -1.0/(GM1))*(((((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0))/(twopi)-(beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*(y-y0))/(twopi))*((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0)*2.0-((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(x-x0, 2.0))/4.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/2.0+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), GM1)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0))/(gamma*1.6e+1-1.6e+1)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*(pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0))/(gamma*3.2e+1))+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(y*4.0-y0*4.0)*pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), GM1)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)-1.0))/1.6e+1-((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(((pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/2.0+pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), gamma)/(GM1))*(y*4.0-y0*4.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, -1.0/(GM1)-1.0))/(gamma*1.6e+1)))/(twopi) - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), gamma)+(((pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/2.0+pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), gamma)/(GM1))*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, -1.0/(GM1)))*(y-y0)*(t*2.0-x*2.0+x0*2.0))/(twopi) - (beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y*2.0-y0*2.0)*(pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), gamma)+(((pow((beta*exp(-pow(y-y0, 2.0)-pow(t-x+x0, 2.0)+1.0)*(y-y0))/(twopi)-1.0, 2.0)+((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*pow(x-x0, 2.0))/4.0)*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)))/2.0+pow(pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, 1.0/(GM1)), gamma)/(GM1))*pow(((beta2)*1.0/(pi2)*exp(pow(y-y0, 2.0)*-2.0-pow(t-x+x0, 2.0)*2.0+2.0)*(GM1)*(-1.0/1.6e+1))/gamma+1.0, -1.0/(GM1)))*(x-x0))/(twopi)
	return
}

func (iv *IVortex) GetFlux(t, x, y float64) (Fx, Fy [4]float64) {
	rho, rhoU, rhoV, rhoE := iv.GetStateC(t, x, y)
	Fx, Fy = FluxCalc(iv.Gamma, rho, rhoU, rhoV, rhoE)
	return
}

func FluxCalc(Gamma, rho, rhoU, rhoV, rhoE float64) (Fx, Fy [4]float64) {
	var (
		GM1 = Gamma - 1.
	)
	u := rhoU / rho
	v := rhoV / rho
	u2 := u*u + v*v
	q := 0.5 * rho * u2
	p := GM1 * (rhoE - q)
	E := rhoE / rho
	Fx, Fy =
		[4]float64{rhoU, rhoU*u + p, rhoU * v, u * (E + p)},
		[4]float64{rhoV, rhoV * u, rhoV*v + p, v * (E + p)}
	return
}
