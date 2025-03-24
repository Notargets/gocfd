package Euler2D

import (
	"fmt"
	"math"
	"strings"

	"github.com/notargets/gocfd/utils"
)

type FlowFunction int16

// Map to store pre-tokenized names
var flowFunctionTokens map[FlowFunction][]string

func init() {
	flowFunctionTokens = make(map[FlowFunction][]string)
	allFunctions := []FlowFunction{
		Density, XMomentum, YMomentum, Energy, Mach, StaticPressure,
		DynamicPressure, PressureCoefficient, SoundSpeed, Velocity,
		XVelocity, YVelocity, Enthalpy, Entropy,
		ShockFunction, EpsilonDissipation, EpsilonDissipationC0,
		XGradientDensity, XGradientXMomentum, XGradientYMomentum, XGradientEnergy,
		YGradientDensity, YGradientXMomentum, YGradientYMomentum, YGradientEnergy,
	}

	for _, ff := range allFunctions {
		tokens := tokenize(ff.String())
		flowFunctionTokens[ff] = tokens
	}
}

// Normalize and tokenize a string into lowercase words
func tokenize(s string) []string {
	s = strings.ToLower(s)
	s = strings.ReplaceAll(s, "-", " ")
	s = strings.ReplaceAll(s, "_", " ")
	return strings.Fields(s)
}

// Match user input to best matching FlowFunction
func BestMatchFlowFunction(input string) (FlowFunction, bool) {
	inputTokens := tokenize(input)

	bestMatch := FlowFunction(-1)
	bestScore := -1

	for ff, tokens := range flowFunctionTokens {
		score := 0
		for _, inputToken := range inputTokens {
			for _, token := range tokens {
				if strings.Contains(token, inputToken) {
					score++
					break
				}
			}
		}
		if score > bestScore {
			bestScore = score
			bestMatch = ff
		}
	}

	return bestMatch, bestScore > 0
}

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
		"Entropy",
	}
	switch {
	case int(pm) < len(strings):
		return strings[int(pm)]
	case pm == ShockFunction:
		return "ShockFunction"
	case pm == EpsilonDissipation:
		return "Artificial Dissipation Epsilon"
	case pm == EpsilonDissipationC0:
		return "Artificial Dissipation Epsilon C0"
	case pm == XGradientDensity:
		return "X Direction Gradient of Density"
	case pm == XGradientXMomentum:
		return "X Direction Gradient of X Momentum"
	case pm == XGradientYMomentum:
		return "X Direction Gradient of Y Momentum"
	case pm == XGradientEnergy:
		return "X Direction Gradient of Energy"
	case pm == YGradientDensity:
		return "Y Direction Gradient of Density"
	case pm == YGradientXMomentum:
		return "Y Direction Gradient of X Momentum"
	case pm == YGradientYMomentum:
		return "Y Direction Gradient of Y Momentum"
	case pm == YGradientEnergy:
		return "Y Direction Gradient of Energy"
	default:
		return "Unknown"
	}
}

const (
	Density FlowFunction = iota
	XMomentum
	YMomentum
	Energy
	Mach                 // 4
	StaticPressure       // 5
	DynamicPressure      // 6
	PressureCoefficient  // 7
	SoundSpeed           // 8
	Velocity             // 9
	XVelocity            // 10
	YVelocity            // 11
	Enthalpy             // 12
	Entropy              //13
	ShockFunction        = 100
	EpsilonDissipation   = 101
	EpsilonDissipationC0 = 102
	XGradientDensity     = 200
	XGradientXMomentum   = 201
	XGradientYMomentum   = 202
	XGradientEnergy      = 203
	YGradientDensity     = 300
	YGradientXMomentum   = 301
	YGradientYMomentum   = 302
	YGradientEnergy      = 303
)

type FreeStream struct {
	Gamma             float64
	Qinf              [4]float64
	Pinf, QQinf, Cinf float64
	Alpha             float64
	Minf              float64
}

func NewFreeStream(Minf, Gamma, Alpha float64) (fs *FreeStream) {
	var (
		ooggm1 = 1. / (Gamma * (Gamma - 1.))
		uinf   = Minf * math.Cos(Alpha*math.Pi/180.)
		vinf   = Minf * math.Sin(Alpha*math.Pi/180.)
	)
	qq := [4]float64{1, uinf, vinf, ooggm1 + 0.5*Minf*Minf}
	fs = &FreeStream{
		Gamma: Gamma,
		Qinf:  qq,
		Alpha: Alpha,
		Minf:  Minf,
	}
	fs.Pinf = fs.GetFlowFunctionQQ(qq, StaticPressure)
	fs.QQinf = fs.GetFlowFunctionQQ(qq, DynamicPressure)
	fs.Cinf = fs.GetFlowFunctionQQ(qq, SoundSpeed)
	return
}

func NewFreestreamFromQinf(gamma float64, qq [4]float64) (fs *FreeStream) {
	fs = &FreeStream{
		Gamma: gamma,
		Qinf:  qq,
	}
	fs.Pinf = fs.GetFlowFunctionQQ(qq, StaticPressure)
	fs.QQinf = fs.GetFlowFunctionQQ(qq, DynamicPressure)
	fs.Cinf = fs.GetFlowFunctionQQ(qq, SoundSpeed)
	return
}

func (fs *FreeStream) Print() string {
	return fmt.Sprintf("Minf[%5.2f] Gamma[%5.2f] Alpha[%5.2f] Q[%8.5f,%8.5f,%8.5f,%8.5f]\n",
		fs.Minf, fs.Gamma, fs.Alpha,
		fs.Qinf[0], fs.Qinf[1], fs.Qinf[2], fs.Qinf[3])
}

func (fs *FreeStream) GetFlowFunction(Q [4]utils.Matrix, ind int, pf FlowFunction) (f float64) {
	return fs.GetFlowFunctionBase(Q[0].DataP[ind], Q[1].DataP[ind], Q[2].DataP[ind], Q[3].DataP[ind], pf)
}

func (fs *FreeStream) GetFlowFunctionQQ(Q [4]float64, pf FlowFunction) (f float64) {
	return fs.GetFlowFunctionBase(Q[0], Q[1], Q[2], Q[3], pf)
}

func (fs *FreeStream) GetFlowFunctionBase(rho, rhoU, rhoV, E float64, pf FlowFunction) (f float64) {
	var (
		Gamma = fs.Gamma
		GM1   = Gamma - 1.
		oorho = 1. / rho
		q, p  float64
	)
	switch pf {
	case Density:
		f = rho
	case XMomentum:
		f = rhoU
	case YMomentum:
		f = rhoV
	case Energy:
		f = E
	case XVelocity:
		f = rhoU * oorho
	case YVelocity:
		f = rhoV * oorho
	case Velocity, DynamicPressure, StaticPressure, PressureCoefficient, SoundSpeed, Enthalpy, Entropy, Mach:
		u, v := rhoU*oorho, rhoV*oorho
		U2 := u*u + v*v
		q = 0.5 * rho * U2
		p = GM1 * (E - q)
		switch pf {
		case Velocity:
			f = math.Sqrt(U2)
		case DynamicPressure:
			f = q
		case StaticPressure:
			f = p
		case PressureCoefficient:
			f = (p - fs.Pinf) / fs.QQinf
		case SoundSpeed:
			f = math.Sqrt(math.Abs(Gamma * p * oorho))
		case Enthalpy:
			f = (E + p) / rho
		case Entropy:
			f = math.Log(p) - Gamma*math.Log(rho)
		case Mach:
			C := math.Sqrt(math.Abs(Gamma * p * oorho))
			U := math.Sqrt(U2)
			f = U / C
		}
	}
	return
}
