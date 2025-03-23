package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/gocfd/model_problems/Euler2D/isentropic_vortex"
)

func PrintU(label string, U [4]float64) {
	fmt.Printf("%s:", label)
	for i := 0; i < 4; i++ {
		fmt.Printf("%0.2f ", U[i])
	}
	fmt.Printf("\n")
}

type TestField uint8

const (
	NORMALSHOCKTESTM2 = TestField(iota)
	NORMALSHOCKTESTM5 = TestField(iota)
	NORMALSHOCKTESTM18
	NORMALSHOCKTESTM12
	FIXEDVORTEXTEST
	MOVINGVORTEXTEST
	RADIAL1TEST
	RADIAL2TEST
	RADIAL3TEST
	RADIAL4TEST
	INTEGERTEST
)

func (tf TestField) String() string {
	switch tf {
	case NORMALSHOCKTESTM2:
		return "NORMALSHOCKTESTM2"
	case NORMALSHOCKTESTM5:
		return "NORMALSHOCKTESTM5"
	case NORMALSHOCKTESTM18:
		return "NORMALSHOCKTESTM18"
	case NORMALSHOCKTESTM12:
		return "NORMALSHOCKTESTM12"
	case FIXEDVORTEXTEST:
		return "FIXEDVORTEXTEST"
	case MOVINGVORTEXTEST:
		return "MOVINGVORTEXTEST"
	case RADIAL1TEST:
		return "RADIAL1TEST"
	case RADIAL2TEST:
		return "RADIAL2TEST"
	case RADIAL3TEST:
		return "RADIAL3TEST"
	case RADIAL4TEST:
		return "RADIAL4TEST"
	}
	return ""
}

func SetTestField(X, Y []float64, tf TestField) (field []float64) {
	switch tf {
	case NORMALSHOCKTESTM12:
		field = setShockConditions(X, 1.2, 0)
	case NORMALSHOCKTESTM18:
		field = setShockConditions(X, 1.8, 0)
	case NORMALSHOCKTESTM2:
		field = setShockConditions(X, 2, 0)
	case NORMALSHOCKTESTM5:
		field = setShockConditions(X, 5, 0)
	case FIXEDVORTEXTEST:
		iv := isentropic_vortex.NewIVortex(5, 0, 0, 1.4, 0)
		field = setIsoVortexConditions(X, Y, iv)
	case MOVINGVORTEXTEST:
		iv := isentropic_vortex.NewIVortex(5, 0, 0, 1.4)
		field = setIsoVortexConditions(X, Y, iv)
	case RADIAL1TEST:
		field = setRadial(X, Y, 1)
	case RADIAL2TEST:
		field = setRadial(X, Y, 2)
	case RADIAL3TEST:
		field = setRadial(X, Y, 3)
	case RADIAL4TEST:
		field = setRadial(X, Y, 4)
	case INTEGERTEST:
		field = setInteger(len(X), 4)
	}
	return
}

func SetTestFieldQ(dfr *DFR2D, tf TestField, Q [4]utils.Matrix) {
	var (
		Np, Kmax = Q[0].Dims()
	)
	X, Y := dfr.SolutionX.Transpose().DataP, dfr.SolutionY.Transpose().DataP
	field := SetTestField(X, Y, tf)
	for n := 0; n < 4; n++ {
		for k := 0; k < Kmax; k++ {
			for i := 0; i < Np; i++ {
				Q[n].Set(i, k, field[i+k*Np+n*Np*Kmax])
			}
		}
	}
	return
}

func setInteger(Np int, order float64) (field []float64) {
	field = make([]float64, Np*4)
	for i := 0; i < Np; i++ {
		for n := 0; n < 4; n++ {
			field[i+n*Np] = float64(n)
		}
	}
	return
}

func setRadial(X, Y []float64, order float64) (field []float64) {
	var (
		Np = len(X)
	)
	field = make([]float64, Np*4)
	var x, y float64
	for i := 0; i < Np; i++ {
		x, y = X[i], Y[i]
		for n := 0; n < 4; n++ {
			field[i+n*Np] = math.Pow(x, order) + math.Pow(y, order)
		}
	}
	return
}

func setShockConditions(X []float64, Mach, Alpha float64) (field []float64) {
	var (
		U1, U2 = ShockConditions(Mach, Alpha)
		Np     = len(X)
	)
	fmt.Printf("Mach: %.2f\n", Mach)
	fmt.Printf("U1: ")
	for n := 0; n < 4; n++ {
		fmt.Printf(" %.2f ", U1[n])
	}
	fmt.Printf("\n")
	fmt.Printf("U2: ")
	for n := 0; n < 4; n++ {
		fmt.Printf(" %.2f ", U2[n])
	}
	fmt.Printf("\n")
	field = make([]float64, Np*4)
	var x float64
	for i := 0; i < Np; i++ {
		x = X[i]
		for n := 0; n < 4; n++ {
			if x < 0 {
				field[i+n*Np] = U1[n]
			} else {
				field[i+n*Np] = U2[n]
			}
		}
	}
	return
}

func setIsoVortexConditions(X, Y []float64,
	iv *isentropic_vortex.IVortex) (field []float64) {
	var (
		Np = len(X)
	)
	field = make([]float64, Np*4)
	var x, y float64
	for i := 0; i < Np; i++ {
		x, y = X[i], Y[i]
		rho, rhoU, rhoV, rhoE := iv.GetStateC(0, x, y)
		field[i+0*Np] = rho
		field[i+1*Np] = rhoU
		field[i+2*Np] = rhoV
		field[i+3*Np] = rhoE
	}
	return
}

// ShockConditions computes post-shock properties given pre-shock Mach number (M1) and shock angle (alpha)
// where alpha = 0 corresponds to a normal shock.
func ShockConditions(M1, alpha float64) (U1, U2 [4]float64) {
	// Alpha = angle of incoming flow relative to the shock front (0 = normal shock).
	var (
		gamma = 1.4
	)

	// Convert angles to radians
	alphaRad := alpha * math.Pi / 180.0

	// Pre-shock conditions (non-dimensional)
	rho1 := 1.0
	p1 := 1.0
	u1 := M1
	v1 := 0.0 // Flow is initially aligned along x (free stream)

	// Resolve incoming velocity into normal/tangential components relative to the shock
	u1n := u1*math.Cos(alphaRad) + v1*math.Sin(alphaRad)  // Normal velocity component
	u1t := -u1*math.Sin(alphaRad) + v1*math.Cos(alphaRad) // Tangential component (unchanged)

	// Normal Mach number
	M1n := u1n / math.Sqrt(gamma*p1/rho1)

	// If subsonic normal component, no shock forms
	if M1n < 1.0 {
		En1 := p1/(gamma-1) + 0.5*rho1*(u1*u1+v1*v1)
		U1 = [4]float64{rho1, rho1 * u1, rho1 * v1, En1}
		U2 = U1
		return
	}

	// Rankine-Hugoniot relations (normal direction)
	rhoRatio := ((gamma + 1) * M1n * M1n) / ((gamma-1)*M1n*M1n + 2)
	pRatio := 1 + (2*gamma/(gamma+1))*(M1n*M1n-1)
	u2n_to_u1n := ((gamma-1)*M1n*M1n + 2) / ((gamma + 1) * M1n * M1n)

	// Post-shock normal values
	rho2 := rho1 * rhoRatio
	p2 := p1 * pRatio
	u2n := u1n * u2n_to_u1n

	// Tangential velocity is unchanged across the shock (slip condition for inviscid flow)
	u2t := u1t

	// Now transform back to x,y coordinates (undo rotation)
	u2 := u2n*math.Cos(alphaRad) - u2t*math.Sin(alphaRad)
	v2 := u2n*math.Sin(alphaRad) + u2t*math.Cos(alphaRad)

	// Total energy (consistent with post-shock pressure and velocity)
	En1 := p1/(gamma-1) + 0.5*rho1*(u1*u1+v1*v1)
	En2 := p2/(gamma-1) + 0.5*rho2*(u2*u2+v2*v2)

	// Pack conservative variables into U1 and U2
	U1 = [4]float64{rho1, rho1 * u1, rho1 * v1, En1}
	U2 = [4]float64{rho2, rho2 * u2, rho2 * v2, En2}

	return
}
