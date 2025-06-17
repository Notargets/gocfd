package gonudg

import (
	"math"
)

// JacobiPNormalized evaluates normalized Jacobi polynomial of type (alpha,beta) > -1
// at points x for order N. This is a direct port of JacobiP.cpp.
// Note: They are normalized to be orthonormal.
func JacobiPNormalized(x []float64, alpha, beta float64, N int) []float64 {
	Nc := len(x)
	P := make([]float64, Nc)

	ab := alpha + beta
	ab1 := alpha + beta + 1.0
	a1 := alpha + 1.0
	b1 := beta + 1.0

	// Initial values P_0(x) and P_1(x)
	gamma0 := math.Pow(2.0, ab1) / ab1 * gamma(a1) * gamma(b1) / gamma(ab1)

	if N == 0 {
		sqrtGamma0 := 1.0 / math.Sqrt(gamma0)
		for i := 0; i < Nc; i++ {
			P[i] = sqrtGamma0
		}
		return P
	}

	// Store polynomials P_0 through P_N
	// PL[i][j] = P_i(x_j)
	PL := make([][]float64, N+1)
	for i := 0; i <= N; i++ {
		PL[i] = make([]float64, Nc)
	}

	// P_0
	sqrtGamma0 := 1.0 / math.Sqrt(gamma0)
	for j := 0; j < Nc; j++ {
		PL[0][j] = sqrtGamma0
	}

	if N >= 1 {
		// P_1
		gamma1 := a1 * b1 / (ab + 3.0) * gamma0
		sqrtGamma1 := 1.0 / math.Sqrt(gamma1)
		for j := 0; j < Nc; j++ {
			PL[1][j] = ((ab+2.0)*x[j]/2.0 + (alpha-beta)/2.0) * sqrtGamma1
		}

		if N == 1 {
			copy(P, PL[1])
			return P
		}
	}

	// Repeat value in recurrence
	aold := 2.0 / (2.0 + ab) * math.Sqrt(a1*b1/(ab+3.0))

	// Forward recurrence using the symmetry of the recurrence
	for i := 1; i <= N-1; i++ {
		h1 := 2.0*float64(i) + ab
		fi := float64(i)
		anew := 2.0 / (h1 + 2.0) * math.Sqrt((fi+1.0)*(fi+ab1)*(fi+a1)*(fi+b1)/(h1+1.0)/(h1+3.0))
		bnew := -(alpha*alpha - beta*beta) / h1 / (h1 + 2.0)

		// PL[i+1] = 1/anew * (-aold*PL[i-1] + (x-bnew)*PL[i])
		for j := 0; j < Nc; j++ {
			PL[i+1][j] = 1.0 / anew * (-aold*PL[i-1][j] + (x[j]-bnew)*PL[i][j])
		}

		aold = anew
	}

	copy(P, PL[N])
	return P
}

// JacobiPNormalizedSingle evaluates normalized Jacobi polynomial at a single point
func JacobiPNormalizedSingle(x, alpha, beta float64, N int) float64 {
	xArr := []float64{x}
	result := JacobiPNormalized(xArr, alpha, beta, N)
	return result[0]
}

// GradJacobiPNormalized evaluates the derivative of normalized Jacobi polynomial
func GradJacobiPNormalized(x []float64, alpha, beta float64, N int) []float64 {
	if N == 0 {
		return make([]float64, len(x))
	}

	// dP_n/dx = sqrt(n(n+alpha+beta+1)) * P_{n-1}^{(alpha+1,beta+1)}(x)
	p := JacobiPNormalized(x, alpha+1, beta+1, N-1)
	fN := float64(N)
	fac := math.Sqrt(fN * (fN + alpha + beta + 1))

	for i := range p {
		p[i] *= fac
	}
	return p
}

// GradJacobiPNormalizedSingle evaluates the derivative at a single point
func GradJacobiPNormalizedSingle(x, alpha, beta float64, N int) float64 {
	xArr := []float64{x}
	result := GradJacobiPNormalized(xArr, alpha, beta, N)
	return result[0]
}

// gamma wraps math.Gamma for clarity
func gamma(x float64) float64 {
	return math.Gamma(x)
}

// The following are compatibility wrappers that redirect to the normalized versions
// These replace the previous DG1D wrappers

// JacobiP evaluates Jacobi polynomial - redirects to normalized version
func JacobiP(x []float64, alpha, beta float64, N int) []float64 {
	return JacobiPNormalized(x, alpha, beta, N)
}

// JacobiPSingle evaluates Jacobi polynomial at a single point
func JacobiPSingle(x, alpha, beta float64, N int) float64 {
	return JacobiPNormalizedSingle(x, alpha, beta, N)
}

// GradJacobiP evaluates the derivative of Jacobi polynomial
func GradJacobiP(x []float64, alpha, beta float64, N int) []float64 {
	return GradJacobiPNormalized(x, alpha, beta, N)
}

// GradJacobiPSingle evaluates the derivative at a single point
func GradJacobiPSingle(x, alpha, beta float64, N int) float64 {
	return GradJacobiPNormalizedSingle(x, alpha, beta, N)
}
