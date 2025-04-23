package Euler2D

import (
	"fmt"
	"math"
	"strings"

	"github.com/notargets/gocfd/utils"
)

func SwitchToEntropyVariables(Q [4]utils.Matrix, gamma float64) {
	/*
		Let:
		  ρ   = density
		  v₁, v₂ = velocity components
		  v²  = v₁² + v₂²
		  E   = total energy
		  p   = pressure
		  T   = temperature
		  s   = entropy = ln(p / ρ^γ)
		  γ   = ratio of specific heats

		Then the entropy variables 𝒘 ∈ ℝ⁴ are:
		          ┌                                       ┐
		          │ (γ - s)/(γ - 1) - v²/(2T)             │
		𝒘 =       │ v₁ / T                                │
		          │ v₂ / T                                │
		          │ -1 / T                                │
		          └                                       ┘
	*/
	var (
		Qp0   = Q[0].DataP
		Qp1   = Q[1].DataP
		Qp2   = Q[2].DataP
		Qp3   = Q[3].DataP
		gm1   = gamma - 1.
		oogm1 = 1. / gm1
	)
	for i := range Qp0 {
		rho, rhoU, rhoV, rhoE := Qp0[i], Qp1[i], Qp2[i], Qp3[i]
		u, v := rhoU/rho, rhoV/rho
		v2 := u*u + v*v
		p := gm1 * (rhoE - 0.5*rho*v2)
		ooT := rho / p
		s := math.Log(p / (math.Pow(rho, gamma)))
		Qp0[i] = (gamma-s)*oogm1 - 0.5*v2*ooT
		Qp1[i] = u * ooT
		Qp2[i] = v * ooT
		Qp3[i] = -ooT
	}
}

func SwitchToConservedVariables(W [4]utils.Matrix, gamma float64) {
	/*
		Given entropy variables:

		          ┌                    ┐
		          │ w₁                 │
		𝒘 =       │ w₂                 │
		          │ w₃                 │
		          │ w₄                 │
		          └                    ┘

		Then compute:
		  T   = -1 / w₄                         (Temperature)
		  v₁  = w₂ / T                          (Velocity component 1)
		  v₂  = w₃ / T                          (Velocity component 2)
		  v²  = v₁² + v₂²                       (Velocity magnitude squared)

		  ρ   = exp( (γ - 1)/γ ⋅ (w₁ - v²/(2T)) )
		  p   = ρ ⋅ T                           (Ideal gas law)
		  E   = p/(γ - 1) + ½ ⋅ ρ ⋅ v²          (Total energy)

		Conserved variables:
		          ┌                     ┐
		          │ ρ                   │
		𝒖 =       │ ρ ⋅ v₁              │
		          │ ρ ⋅ v₂              │
		          │ E                   │
		          └                     ┘
	*/
	var (
		Wp0    = W[0].DataP
		Wp1    = W[1].DataP
		Wp2    = W[2].DataP
		Wp3    = W[3].DataP
		gm1    = gamma - 1.
		oogm1  = 1. / gm1
		gm1ogm = gm1 / gamma
	)
	for i := range Wp0 {
		w1, w2, w3, w4 := Wp0[i], Wp1[i], Wp2[i], Wp3[i]
		T := -1. / w4
		ooT := 1. / T
		u, v := w2*ooT, w3*ooT
		v2 := u*u + v*v
		rho := math.Exp(gm1ogm * (w1 - 0.5*v2*ooT))
		p := rho * T
		Wp0[i] = rho
		Wp1[i] = rho * u
		Wp2[i] = rho * v
		Wp3[i] = p*oogm1 + 0.5*rho*v2
	}
}

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
		return "R Direction Gradient of Density"
	case pm == XGradientXMomentum:
		return "R Direction Gradient of R Momentum"
	case pm == XGradientYMomentum:
		return "R Direction Gradient of S Momentum"
	case pm == XGradientEnergy:
		return "R Direction Gradient of Energy"
	case pm == YGradientDensity:
		return "S Direction Gradient of Density"
	case pm == YGradientXMomentum:
		return "S Direction Gradient of R Momentum"
	case pm == YGradientYMomentum:
		return "S Direction Gradient of S Momentum"
	case pm == YGradientEnergy:
		return "S Direction Gradient of Energy"
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
	Entropy              // 13
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
