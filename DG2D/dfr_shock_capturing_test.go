package DG2D

import (
	"fmt"
	"math"
	"testing"
)

func TestPlotEquiTri(t *testing.T) {
	dfr := CreateEquiTriMesh(1, 50)
	_ = dfr
	if testing.Verbose() {
		PlotDFRElements(dfr)
	}
}

func TestMachConditions(t *testing.T) {
	U1, U2 := ShockConditions(2, 0)
	fmt.Printf("U1: %5.2f,%5.2f,%5.2f,%5.2f\n", U1[0], U1[1], U1[2], U1[3])
	fmt.Printf("U2: %5.2f,%5.2f,%5.2f,%5.2f\n", U2[0], U2[1], U2[2], U2[3])
	U1, U2 = ShockConditions(2, 15)
	fmt.Printf("U1: %5.2f,%5.2f,%5.2f,%5.2f\n", U1[0], U1[1], U1[2], U1[3])
	fmt.Printf("U2: %5.2f,%5.2f,%5.2f,%5.2f\n", U2[0], U2[1], U2[2], U2[3])
}

// ShockConditions computes post-shock properties given pre-shock Mach number (M1) and shock angle (alpha)
// where alpha = 0 corresponds to a normal shock.
func ShockConditions(M1, alpha float64) (U1, U2 [4]float64) {
	// Alpha = 0 for a normal shock.
	// Alpha can be used to determine oblique shock results in inviscid flow
	// accurately
	// For normal shock, set beta = alpha + 90 degrees.
	// Here, beta is the shock angle relative to the incoming flow.
	var (
		gamma   = 1.4
		beta    = alpha + 90.
		betaRad = beta * math.Pi / 180.0 // Convert beta to radians
	)

	// Compute normal component of Mach number.
	M1n := M1 * math.Sin(betaRad)

	// Compute density ratio from normal shock relations.
	rhoRatio := ((gamma + 1) * M1n * M1n) / ((gamma-1)*M1n*M1n + 2)

	// Compute pressure ratio from normal shock relations.
	pRatio := 1 + (2*gamma/(gamma+1))*(M1n*M1n-1)

	// Pre-shock conditions (assumed non-dimensional for testing).
	rho1 := 1.0 // Reference density.
	p1 := 1.0   // Reference pressure.
	u1 := M1    // Freestream velocity (assumed normalized so that a1 = 1).
	v1 := 0.0   // No initial vertical velocity.

	// Define an inline Newton solver to compute u2n.
	NewtonSolver := func(rho1, u1n, p1, rho2, p2, tol float64, maxIter int) float64 {
		u2n := u1n / 2.0 // Initial guess.
		for i := 0; i < maxIter; i++ {
			// Correct residual: F = (rho2*u2n^2 + p2) - (rho1*u1n^2 + p1)
			F := (rho2*u2n*u2n + p2) - (rho1*u1n*u1n + p1)
			// Derivative: F' = 2 * rho2 * u2n.
			Fprime := 2.0 * rho2 * u2n

			// Newton update.
			u2nNew := u2n - F/Fprime

			// Check convergence.
			if math.Abs(u2nNew-u2n) < tol {
				return u2nNew
			}
			u2n = u2nNew
		}
		fmt.Println("Warning: Newton solver did not converge!")
		return u2n
	}

	// Compute post-shock properties.
	rho2 := rho1 * rhoRatio
	p2 := p1 * pRatio
	// Tangential component remains unchanged.
	u1t := u1 * math.Cos(betaRad)

	// Solve for post-shock normal velocity using the Newton solver.
	u2n := NewtonSolver(rho1, M1n, p1, rho2, p2, 1e-6, 50)

	// Compute full post-shock velocity components.
	u2 := u2n*math.Sin(betaRad) + u1t*math.Cos(betaRad)
	v2 := u2n*math.Cos(betaRad) - u1t*math.Sin(betaRad)

	// Compute total energy (using post-shock pressure and velocity).
	En1 := p1/(gamma-1) + 0.5*rho1*(u1*u1+v1*v1)
	En2 := p2/(gamma-1) + 0.5*rho2*(u2*u2+v2*v2)

	// Pack pre- and post-shock states into U1 and U2.
	U1 = [4]float64{rho1, u1, v1, En1}
	U2 = [4]float64{rho2, u2, v2, En2}
	return
}
