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
