package DG1D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/stretchr/testify/assert"
)

func TestJacobiGQ_PartitionAndFirstMoment(t *testing.T) {
	const (
		α   = 0.3
		β   = 0.7
		N   = 5
		tol = 1e-12
	)

	// 1) Generate nodes & weights
	X, Wvec := JacobiGQ(α, β, N)
	x, w := X.DataP, Wvec.DataP

	// 2) Compute the analytic total weight = ∫_{-1}^1 (1-x)^α (1+x)^β dx
	// This equals 2^{α+β+1} * B(α+1, β+1) where B is the Beta function
	// B(α+1, β+1) = Γ(α+1) * Γ(β+1) / Γ(α+β+2)
	exactZero := math.Pow(2, α+β+1) * beta(α+1, β+1)

	// 3) Compute the analytic first moment = ∫_{-1}^1 x * (1-x)^α (1+x)^β dx
	// This equals (β-α)/(α+β+2) * ∫_{-1}^1 (1-x)^α (1+x)^β dx
	exactOne := (β - α) / (α + β + 2) * exactZero

	// 4) Sum using discrete quadrature rule
	var sum0, sum1 float64
	for i := range x {
		sum0 += w[i]
		sum1 += x[i] * w[i]
	}

	assert.InDeltaf(t, exactZero, sum0, tol,
		"sum(w) = %v, want %v", sum0, exactZero)
	assert.InDeltaf(t, exactOne, sum1, tol,
		"sum(x*w) = %v, want %v", sum1, exactOne)
}

func TestJacobiGQ_RootsAndMoments(t *testing.T) {
	const (
		α   = 0.3
		β   = 0.7
		N   = 5
		tol = 1e-10
	)

	// 1) Build your JacobiGQ
	Xvec, Wvec := JacobiGQ(α, β, N)
	X, W := Xvec.DataP, Wvec.DataP
	assert.Equal(t, N+1, len(X))
	assert.Equal(t, N+1, len(W))

	// 2) Root test: P_{N+1}^{(α,β)}(X[i]) ≈ 0
	// The nodes should be roots of the (N+1)-th Jacobi polynomial
	for i, xi := range X {
		pi := JacobiP(
			utils.NewVector(1, []float64{xi}), α, β, N+1,
		)[0]
		assert.InDeltaf(t, 0, pi, 1e-10,
			"JacobiP_{%d}^{(%.1f,%.1f)}(%g) = %g ≠ 0 (node %d)", N+1, α, β, xi, pi, i)
	}

	// 3) Moment test: ∑ W[i]*X[i]^k == exactMoment(k) for k = 0, 1, ..., 2N+1
	// Jacobi-Gauss quadrature should integrate polynomials up to degree 2N+1 exactly
	for k := 0; k <= 2*N+1; k++ {
		var s float64
		for i, xi := range X {
			s += W[i] * math.Pow(xi, float64(k))
		}
		m := exactMoment(k, α, β)
		assert.InDeltaf(t, m, s, tol,
			"moment %d: got %g, want %g", k, s, m)
	}

	// 4) Additional validation: check that nodes are in (-1, 1)
	for i, xi := range X {
		assert.True(t, xi > -1 && xi < 1,
			"Node %d: x_%d = %g should be in (-1, 1)", i, i, xi)
	}

	// 5) Additional validation: weights should be positive (for α, β > -1)
	if α > -1 && β > -1 {
		for i, wi := range W {
			assert.True(t, wi > 0,
				"Weight %d: w_%d = %g should be positive", i, i, wi)
		}
	}
}

// exactMoment computes ∫_{-1}^1 x^k (1-x)^α (1+x)^β dx analytically
func exactMoment(k int, α, β float64) float64 {
	if k < 0 {
		return 0
	}

	// Use binomial expansion: (1+x)^β expanded and integrated term by term
	// ∫_{-1}^1 x^k (1-x)^α (1+x)^β dx
	//
	// Method: Use substitution u = (1+x)/2, so x = 2u-1, dx = 2du
	// Integral becomes: 2^{α+β+1} ∫_0^1 (2u-1)^k u^β (1-u)^α du
	//
	// Expand (2u-1)^k using binomial theorem:
	// (2u-1)^k = ∑_{j=0}^k C(k,j) (2u)^j (-1)^{k-j}
	//          = ∑_{j=0}^k C(k,j) 2^j (-1)^{k-j} u^j

	var result float64
	for j := 0; j <= k; j++ {
		// Coefficient from binomial expansion
		coeff := float64(choose(k, j)) * math.Pow(2, float64(j)) * math.Pow(-1, float64(k-j))

		// ∫_0^1 u^{j+β} (1-u)^α du = B(j+β+1, α+1)
		betaIntegral := beta(float64(j)+β+1, α+1)

		result += coeff * betaIntegral
	}

	// Multiply by the factor 2^{α+β+1} from the substitution
	result *= math.Pow(2, α+β+1)

	return result
}

func TestJacobiGQ_SanityChecks(t *testing.T) {
	const (
		α   = 0.3
		β   = 0.7
		N   = 5 // choose a small N
		tol = 1e-10
	)

	// 1) Build your rule
	Xv, Wv := JacobiGQ(α, β, N)
	x, w := Xv.DataP, Wv.DataP
	assert.Equal(t, N+1, len(x))
	assert.Equal(t, N+1, len(w))

	// 2) Partition of unity test
	// For Jacobi-Gauss quadrature, ∑w_i = ∫_{-1}^1 (1-x)^α (1+x)^β dx
	// This integral equals 2^{α+β+1} * B(α+1, β+1)
	exact0 := math.Pow(2, α+β+1) * beta(α+1, β+1)

	var sum0 float64
	for i := range x {
		sum0 += w[i]
	}
	assert.InDeltaf(t, exact0, sum0, tol,
		"∑w = %g, want %g", sum0, exact0)

	// 3) First moment test
	// ∫_{-1}^1 x * (1-x)^α (1+x)^β dx
	// Using the identity: ∫_{-1}^1 x * (1-x)^α (1+x)^β dx = (β-α)/(α+β+2) * ∫_{-1}^1 (1-x)^α (1+x)^β dx
	exact1 := (β - α) / (α + β + 2) * exact0

	var sum1 float64
	for i := range x {
		sum1 += x[i] * w[i]
	}
	assert.InDeltaf(t, exact1, sum1, tol,
		"∑xw = %g, want %g", sum1, exact1)

	// 4) Root test: P_{N+1}^{(α,β)}(x_i) ≈ 0
	// The nodes should be roots of the (N+1)-th Jacobi polynomial
	for i, xi := range x {
		Pi1 := JacobiP(
			utils.NewVector(1, []float64{xi}),
			α, β, N+1,
		)[0]
		assert.InDeltaf(t, 0, Pi1, 1e-10,
			"P_{%d}^{(%.1f,%.1f)}(%.6f) = %g ≠ 0 (node %d)", N+1, α, β, xi, Pi1, i)
	}

	// 5) Exactness test for polynomials up to degree 2N+1
	// Jacobi-Gauss quadrature should integrate polynomials of degree ≤ 2N+1 exactly
	// We test moments: ∫_{-1}^1 x^k * (1-x)^α * (1+x)^β dx
	for k := 0; k <= 2*N+1; k++ {
		// Analytical moment using binomial expansion:
		// (1+x)^β = ∑_{j=0}^β C(β,j) x^j, but since β is not integer, we use:
		// ∫_{-1}^1 x^k (1-x)^α (1+x)^β dx = ∑_{j=0}^k C(k,j) (-1)^{k-j} ∫_{-1}^1 x^j (1+x)^{α+β+k-j} dx
		//
		// Actually, the correct approach is to use the substitution u = (1+x)/2:
		// ∫_{-1}^1 x^k (1-x)^α (1+x)^β dx = 2^{α+β+1} ∫_0^1 (2u-1)^k u^β (1-u)^α du

		var exactK float64

		// Use binomial expansion of (2u-1)^k = ∑_{j=0}^k C(k,j) (2u)^j (-1)^{k-j}
		for j := 0; j <= k; j++ {
			coeff := float64(choose(k, j)) * math.Pow(-1, float64(k-j)) * math.Pow(2, float64(j))
			// ∫_0^1 u^{j+β} (1-u)^α du = B(j+β+1, α+1)
			betaIntegral := beta(float64(j)+β+1, α+1)
			exactK += coeff * betaIntegral
		}
		exactK *= math.Pow(2, α+β+1)

		// Discrete moment using quadrature
		var sumK float64
		for i := range x {
			sumK += w[i] * math.Pow(x[i], float64(k))
		}

		assert.InDeltaf(t, exactK, sumK, tol,
			"moment %d: ∑w_i x_i^k = %g, want %g", k, sumK, exactK)
	}

	// 6) Additional test: verify nodes are in (-1, 1)
	for i, xi := range x {
		assert.True(t, xi > -1 && xi < 1,
			"Node %d: x_%d = %g should be in (-1, 1)", i, i, xi)
	}

	// 7) Additional test: verify weights are positive (for α, β > -1)
	if α > -1 && β > -1 {
		for i, wi := range w {
			assert.True(t, wi > 0,
				"Weight %d: w_%d = %g should be positive", i, i, wi)
		}
	}

	// 8) Symmetry test (when α = β)
	if math.Abs(α-β) < 1e-14 {
		// When α = β, the quadrature should be symmetric about x = 0
		for i := 0; i < len(x)/2; i++ {
			j := len(x) - 1 - i
			assert.InDeltaf(t, -x[i], x[j], tol,
				"Symmetry: x_%d = %g should equal -x_%d = %g", i, x[i], j, -x[j])
			assert.InDeltaf(t, w[i], w[j], tol,
				"Symmetry: w_%d = %g should equal w_%d = %g", i, w[i], j, w[j])
		}
	}
}

// Helper function for binomial coefficients
func choose(n, k int) int {
	if k > n || k < 0 {
		return 0
	}
	if k == 0 || k == n {
		return 1
	}
	if k > n-k {
		k = n - k
	}

	result := 1
	for i := 0; i < k; i++ {
		result = result * (n - i) / (i + 1)
	}
	return result
}

// Helper function for Beta function
func beta(a, b float64) float64 {
	return math.Gamma(a) * math.Gamma(b) / math.Gamma(a+b)
}

func TestJacobiPOrthogonality_WithJacobiGQ(t *testing.T) {
	const (
		α    = 0.3
		β    = 0.7
		Nmax = 6
	)

	// 1) Build Gauss-Jacobi nodes & weights with sufficient accuracy
	// To test orthogonality for polynomials up to degree Nmax, we need to integrate
	// products P_m * P_n where m,n ≤ Nmax. The highest degree product is 2*Nmax.
	// For safety, we use more quadrature points than strictly necessary.
	quadOrder := Nmax + 3 // This gives us exactness up to degree 2*(Nmax+3)+1 = 2*Nmax+7
	X, Wvec := JacobiGQ(α, β, quadOrder)
	nodes := X.DataP // []float64 of length quadOrder+1
	weights := Wvec.DataP

	// 2) Evaluate P_n at each quadrature node for n = 0, 1, ..., Nmax
	Pvals := make([][]float64, Nmax+1)
	for n := 0; n <= Nmax; n++ {
		// JacobiP takes a utils.Vector, so wrap nodes
		xv := utils.NewVector(len(nodes), nodes)
		Pvals[n] = JacobiP(xv, α, β, n)
	}

	// 3) Form the Gram matrix G[m][n] = ∫_{-1}^1 P_m(x) P_n(x) (1-x)^α (1+x)^β dx
	// This is approximated by the quadrature sum: ∑_i w_i * P_m(x_i) * P_n(x_i)
	G := make([][]float64, Nmax+1)
	for m := 0; m <= Nmax; m++ {
		G[m] = make([]float64, Nmax+1)
		for n := 0; n <= Nmax; n++ {
			var sum float64
			for i := range nodes {
				sum += weights[i] * Pvals[m][i] * Pvals[n][i]
			}
			G[m][n] = sum
		}
	}

	// 4) Check orthogonality: ∫ P_m P_n w(x) dx = 0 for m ≠ n
	// and positivity of diagonal elements: ∫ P_n^2 w(x) dx > 0
	tol := 1e-10
	for m := 0; m <= Nmax; m++ {
		for n := 0; n <= Nmax; n++ {
			if m != n {
				assert.InDeltaf(t, 0, G[m][n], tol,
					"Orthogonality failed: ∫P_%d^{(%.1f,%.1f)} P_%d^{(%.1f,%.1f)} w(x) dx = %g ≠ 0",
					m, α, β, n, α, β, G[m][n])
			} else {
				assert.Truef(t, G[m][n] > 0,
					"Normalization failed: ∫[P_%d^{(%.1f,%.1f)}]^2 w(x) dx = %g should be positive",
					m, α, β, G[m][n])

				// Check that the polynomials are orthonormal (norm = 1)
				// Your JacobiP function appears to return orthonormalized polynomials
				expectedNorm := 1.0
				assert.InDeltaf(t, expectedNorm, G[m][n], tol,
					"Orthonormality failed: ∫[P_%d^{(%.1f,%.1f)}]^2 w(x) dx = %g, expected 1.0",
					m, α, β, G[m][n])
			}
		}
	}

	// 5) Additional validation: Check that the quadrature has sufficient accuracy
	// Test with a known polynomial of degree 2*Nmax
	testDegree := 2 * Nmax
	if testDegree <= 2*quadOrder+1 {
		// Integrate x^{testDegree} * (1-x)^α * (1+x)^β exactly vs numerically
		exactIntegral := exactMoment(testDegree, α, β)
		var numericalIntegral float64
		for i := range nodes {
			numericalIntegral += weights[i] * math.Pow(nodes[i], float64(testDegree))
		}
		assert.InDeltaf(t, exactIntegral, numericalIntegral, tol*math.Max(1, math.Abs(exactIntegral)),
			"Quadrature accuracy check failed for degree %d polynomial", testDegree)
	}

	t.Logf("Successfully verified orthonormality for Jacobi polynomials P_n^{(%.1f,%.1f)} with n = 0,...,%d",
		α, β, Nmax)
	t.Logf("Used %d quadrature points (exact up to degree %d)", len(nodes), 2*quadOrder+1)
}
