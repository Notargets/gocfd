package Euler2D

import (
	"fmt"
	"strings"
)

func (it InitType) Print() (txt string) {
	txt = InitPrintNames[it]
	return
}

type FluxType uint

const (
	FLUX_Average FluxType = iota
	FLUX_LaxFriedrichs
	FLUX_Roe
)

var (
	FluxNames = map[string]FluxType{
		"average": FLUX_Average,
		"lax":     FLUX_LaxFriedrichs,
		"roe":     FLUX_Roe,
	}
	FluxPrintNames = []string{"Average", "Lax Friedrichs", "Roe"}
)

func (ft FluxType) Print() (txt string) {
	txt = FluxPrintNames[ft]
	return
}

func NewFluxType(label string) (ft FluxType) {
	var (
		ok  bool
		err error
	)
	label = strings.ToLower(label)
	if ft, ok = FluxNames[label]; !ok {
		err = fmt.Errorf("unable to use flux named %s", label)
		panic(err)
	}
	return
}

func (c *Euler) CalculateFluxTransformed(k, i int, QQ [4][]float64) (Fr, Fs [4]float64) {
	var (
		Jdet = c.dfr.Jdet.Data()[k]
		Jinv = c.dfr.Jinv.Data()[4*k : 4*(k+1)]
	)
	Fx, Fy := c.CalculateFlux(k, i, QQ)
	for n := 0; n < 4; n++ {
		Fr[n] = Jdet * (Jinv[0]*Fx[n] + Jinv[1]*Fy[n])
		Fs[n] = Jdet * (Jinv[2]*Fx[n] + Jinv[3]*Fy[n])
	}
	return
}

func (c *Euler) CalculateFlux(k, i int, QQ [4][]float64) (Fx, Fy [4]float64) {
	// From https://www.theoretical-physics.net/dev/fluid-dynamics/euler.html
	var (
		Kmax = c.dfr.K
		ind  = k + i*Kmax
		qq   = [4]float64{QQ[0][ind], QQ[1][ind], QQ[2][ind], QQ[3][ind]}
	)
	Fx, Fy = c.FluxCalcMock(qq)
	return
}

func (c *Euler) FluxCalc(q [4]float64) (Fx, Fy [4]float64) {
	var (
		rho, rhoU, rhoV, E = q[0], q[1], q[2], q[3]
		oorho              = 1. / rho
		u                  = rhoU * oorho
		v                  = rhoV * oorho
		p                  = c.FS.GetFlowFunction(q, StaticPressure)
	)
	Fx, Fy =
		[4]float64{rhoU, rhoU*u + p, rhoU * v, u * (E + p)},
		[4]float64{rhoV, rhoV * u, rhoV*v + p, v * (E + p)}
	return
}
