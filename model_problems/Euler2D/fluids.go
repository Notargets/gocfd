package Euler2D

import "math"

type FlowFunction uint8

func (pm FlowFunction) String() string {
	strings := []string{
		"Density",
		"XMomentum",
		"YMomentum",
		"Energy",
		"Mach",
		"Static Pressure",
		"Dynamic Pressure",
		"Pressure Coefficient",
		"Sound Speed",
		"Velocity",
		"XVelocity",
		"YVelocity",
		"Enthalpy",
	}
	return strings[int(pm)]
}

const (
	Density FlowFunction = iota
	XMomentum
	YMomentum
	Energy
	Mach                // 4
	StaticPressure      // 5
	DynamicPressure     // 6
	PressureCoefficient // 7
	SoundSpeed          // 8
	Velocity            // 9
	XVelocity           // 10
	YVelocity           // 11
	Enthalpy            // 12
)

func (c *Euler) GetFlowFunction(Q [4]float64, pf FlowFunction) (f float64) {
	var (
		Gamma              = c.Gamma
		GM1                = Gamma - 1.
		rho, rhou, rhov, E = Q[0], Q[1], Q[2], Q[3]
		oorho              = 1. / rho
		q, p               float64
	)
	// Calculate q if needed
	switch pf {
	case StaticPressure, PressureCoefficient, SoundSpeed:
		q = 0.5 * (rhou*rhou + rhov*rhov) * oorho
	}
	// Calculate p if needed
	switch pf {
	case PressureCoefficient, SoundSpeed, Enthalpy, Mach:
		p = GM1 * (E - q)
	}

	switch pf {
	case Density:
		f = rho
	case XMomentum:
		f = rhou
	case YMomentum:
		f = rhov
	case Energy:
		f = E
	case StaticPressure:
		f = GM1 * (E - q)
	case DynamicPressure:
		f = 0.5 * (rhou*rhou + rhov*rhov) * oorho
	case PressureCoefficient:
		f = (p - c.Pinf) / c.QQinf
	case SoundSpeed:
		f = math.Sqrt(math.Abs(Gamma * p * oorho))
	case Velocity:
		f = math.Sqrt((rhou*rhou + rhov*rhov) * oorho)
	case XVelocity:
		f = rhou * oorho
	case YVelocity:
		f = rhov * oorho
	case Mach:
		C := math.Sqrt(math.Abs(Gamma * p * oorho))
		U := math.Sqrt((rhou*rhou + rhov*rhov)) * oorho
		f = U / C
	case Enthalpy:
		f = (E + p) / rho
	}
	return
}

func (c *Euler) GetFlowFunctionAtIndex(ind int, QQ [4][]float64, pf FlowFunction) (f float64) {
	var (
		QI = [4]float64{QQ[0][ind], QQ[1][ind], QQ[2][ind], QQ[3][ind]}
	)
	f = c.GetFlowFunction(QI, pf)
	return
}
