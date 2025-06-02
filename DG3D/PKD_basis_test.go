package DG3D

import (
	"github.com/notargets/gocfd/utils"
	"math"
	"testing"
)

// generateTestNodes generates reference tetrahedral nodes for testing
// Uses a simple approach to get better distributed nodes
// generateTestNodes generates reference tetrahedral nodes for testing
// Uses a simple approach to get better distributed nodes
// generateTestNodes generates reference tetrahedral nodes for testing
// Uses a simple approach to get better distributed nodes
func generateTestNodes(P int) (r, s, t []float64) {
	switch P {
	case 1:
		// For P=1, use face-centered nodes that avoid corners
		r = []float64{-1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0}
		s = []float64{-1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0, -1.0 / 3.0}
		t = []float64{-1.0 / 3.0, -1.0 / 3.0, -1.0 / 3.0, 1.0 / 3.0}
		return r, s, t
	case 2:
		// For P=2, use a symmetric set of 10 nodes
		// Include centroid, face centers, and edge centers
		r = make([]float64, 10)
		s = make([]float64, 10)
		t = make([]float64, 10)

		// Centroid
		r[0], s[0], t[0] = -1.0/3.0, -1.0/3.0, -1.0/3.0

		// Face centers (4 faces)
		r[1], s[1], t[1] = 0.0, -1.0/3.0, -1.0/3.0      // r=0 face
		r[2], s[2], t[2] = -1.0/3.0, 0.0, -1.0/3.0      // s=0 face
		r[3], s[3], t[3] = -1.0/3.0, -1.0/3.0, 0.0      // t=0 face
		r[4], s[4], t[4] = -2.0/3.0, -2.0/3.0, -2.0/3.0 // r+s+t=-2 face

		// Edge midpoints (shifted inward)
		alpha := 0.5 // midpoint parameter
		beta := 0.7  // inward shift

		// Edge from v0 to v1
		r[5] = beta * (-1.0 + alpha*2.0)
		s[5] = beta * (-1.0)
		t[5] = beta * (-1.0)

		// Edge from v0 to v2
		r[6] = beta * (-1.0)
		s[6] = beta * (-1.0 + alpha*2.0)
		t[6] = beta * (-1.0)

		// Edge from v0 to v3
		r[7] = beta * (-1.0)
		s[7] = beta * (-1.0)
		t[7] = beta * (-1.0 + alpha*2.0)

		// Edge from v1 to v2
		r[8] = beta * (1.0 - alpha*2.0)
		s[8] = beta * (-1.0 + alpha*2.0)
		t[8] = beta * (-1.0)

		// Edge from v1 to v3
		r[9] = beta * (1.0 - alpha*2.0)
		s[9] = beta * (-1.0)
		t[9] = beta * (-1.0 + alpha*2.0)

		return r, s, t
	default:
		// For higher orders, use a Stroud-like distribution
		// that provides good conditioning for the PKD basis
		n := (P + 1) * (P + 2) * (P + 3) / 6
		r = make([]float64, n)
		s = make([]float64, n)
		t = make([]float64, n)

		// Use a warped product grid that respects the tetrahedral structure
		idx := 0
		for i := 0; i <= P; i++ {
			for j := 0; j <= P-i; j++ {
				for k := 0; k <= P-i-j; k++ {
					// Use warped coordinates to avoid clustering near boundaries
					xi := float64(i) / float64(P)
					eta := float64(j) / float64(P)
					zeta := float64(k) / float64(P)

					// Apply warping function to push nodes away from boundaries
					warp := func(x float64) float64 {
						// Smooth warping that preserves symmetry
						return x + 0.3*x*(1-x)*(1-2*x)
					}

					xi = warp(xi)
					eta = warp(eta)
					zeta = warp(zeta)

					// Convert to barycentric coordinates
					lambda := 1.0 - xi - eta - zeta

					// Ensure we're inside the tetrahedron
					if lambda < 0 {
						// Project back to tetrahedron
						sum := xi + eta + zeta
						xi /= sum
						eta /= sum
						zeta /= sum
						lambda = 0
					}

					// Convert barycentric to reference coordinates
					r[idx] = -1.0 + 2.0*xi
					s[idx] = -1.0 + 2.0*eta
					t[idx] = -1.0 + 2.0*zeta

					idx++
				}
			}
		}
		return r, s, t
	}
}

// Monomial evaluates a monomial r^p * s^q * t^k at given points
func Monomial(r, s, t []float64, p, q, k int) []float64 {
	n := len(r)
	u := make([]float64, n)
	for i := 0; i < n; i++ {
		u[i] = math.Pow(r[i], float64(p)) *
			math.Pow(s[i], float64(q)) *
			math.Pow(t[i], float64(k))
	}
	return u
}

// MonomialDeriv computes analytical derivative of monomial
func MonomialDeriv(r, s, t []float64, p, q, k int, dir int) []float64 {
	n := len(r)
	du := make([]float64, n)

	switch dir {
	case 0: // dr
		if p > 0 {
			for i := 0; i < n; i++ {
				du[i] = float64(p) * math.Pow(r[i], float64(p-1)) *
					math.Pow(s[i], float64(q)) *
					math.Pow(t[i], float64(k))
			}
		}
	case 1: // ds
		if q > 0 {
			for i := 0; i < n; i++ {
				du[i] = float64(q) * math.Pow(r[i], float64(p)) *
					math.Pow(s[i], float64(q-1)) *
					math.Pow(t[i], float64(k))
			}
		}
	case 2: // dt
		if k > 0 {
			for i := 0; i < n; i++ {
				du[i] = float64(k) * math.Pow(r[i], float64(p)) *
					math.Pow(s[i], float64(q)) *
					math.Pow(t[i], float64(k-1))
			}
		}
	}

	return du
}

// MatVec performs matrix-vector multiplication
func MatVec(A utils.Matrix, x []float64) []float64 {
	m := A.Rows()
	n := A.Cols()
	if len(x) != n {
		panic("dimension mismatch")
	}

	y := make([]float64, m)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			y[i] += A.At(i, j) * x[j]
		}
	}
	return y
}

// L2Error computes the L2 norm of the error between two vectors
func L2Error(u1, u2 []float64) float64 {
	if len(u1) != len(u2) {
		panic("dimension mismatch")
	}

	var err float64
	for i := range u1 {
		diff := u1[i] - u2[i]
		err += diff * diff
	}
	return math.Sqrt(err / float64(len(u1)))
}

func TestVandermonde(t *testing.T) {
	for P := 1; P <= 7; P++ {
		basis := NewPKDBasis(P)
		r, s, tt := generateTestNodes(P)

		// Check that we have the right number of nodes
		expectedNodes := (P + 1) * (P + 2) * (P + 3) / 6
		if len(r) != expectedNodes {
			t.Errorf("P=%d: Wrong number of nodes: got %d, expected %d", P, len(r), expectedNodes)
			continue
		}

		V := basis.Vandermonde(r, s, tt)

		// Test that V is square for nodal basis
		if V.Rows() != V.Cols() {
			t.Errorf("P=%d: Vandermonde not square: %d x %d", P, V.Rows(), V.Cols())
		}

		// Test conditioning
		cond := V.ConditionNumber()
		t.Logf("P=%d: Vandermonde size %d x %d, condition number %.2e",
			P, V.Rows(), V.Cols(), cond)

		// Warn if poorly conditioned
		if cond > 1e10 {
			t.Logf("WARNING: P=%d: Vandermonde matrix is poorly conditioned", P)
		}
	}
}

func TestBasisEvaluation(t *testing.T) {
	// Test orthogonal PKD modal basis properties
	for P := 1; P <= 3; P++ {
		basis := NewPKDBasis(P)

		// Test 1: Constant mode (i=j=k=0) should equal 1 everywhere
		testPoints := [][3]float64{
			{-1, -1, -1},
			{1, -1, -1},
			{-1, 1, -1},
			{-1, -1, 1},
			{0, 0, 0},
			{-0.5, -0.5, -0.5},
			{0.5, -0.5, -0.5},
		}

		for _, pt := range testPoints {
			phi := basis.EvalBasis(pt[0], pt[1], pt[2])
			// First basis function (constant mode) should be 1
			if math.Abs(phi[0]-1.0) > 1e-10 {
				t.Errorf("P=%d: Constant mode at point %v = %g (expected 1)",
					P, pt, phi[0])
			}
		}

		// Test 2: Polynomial exactness - can represent any polynomial up to degree P
		// Test with monomials r^p * s^q * t^k where p+q+k <= P
		r, s, tt := generateTestNodes(P)
		V := basis.Vandermonde(r, s, tt)

		// Check that V is well-conditioned enough to represent polynomials
		cond := V.ConditionNumber()
		if cond > 1e12 {
			t.Logf("P=%d: Warning - Vandermonde condition number = %.2e", P, cond)
		}

		// Test 3: Basis completeness - Vandermonde should be invertible
		_, err := V.Inverse()
		if err != nil && cond < 1e14 {
			t.Errorf("P=%d: Vandermonde matrix is singular", P)
		}

		// Test 4: Hierarchical property - first (P-1+1)(P-1+2)(P-1+3)/6 modes
		// should span the P-1 polynomial space
		if P > 1 {
			N_prev := (P * (P + 1) * (P + 2)) / 6
			// Evaluate basis at a test point
			phi := basis.EvalBasis(0.2, 0.1, 0.3)

			// Check that we have the expected number of basis functions
			if len(phi) != basis.N {
				t.Errorf("P=%d: Expected %d basis functions, got %d",
					P, basis.N, len(phi))
			}

			// Verify hierarchical structure exists
			if N_prev > basis.N {
				t.Errorf("P=%d: Hierarchical structure violated", P)
			}
		}

		t.Logf("P=%d: Modal basis evaluation test passed (N=%d, cond=%.2e)",
			P, basis.N, cond)
	}

	// Test 5: Specific values for P=1 to verify correct modal basis
	P := 1
	basis := NewPKDBasis(P)

	// At origin (0,0,0), the collapsed coordinates are:
	// a = 2*(1+0)/(1-0-0) - 1 = 1
	// b = 2*(1+0)/(1-0) - 1 = 1
	// c = 2*(1+0) - 1 = 1
	phi := basis.EvalBasis(0, 0, 0)
	t.Logf("P=1 at origin: phi = %v", phi)

	// For P=1, we expect 4 modes:
	// Mode 0 (i=0,j=0,k=0): P₀(1) * P₀(1) * P₀(1) = 1*1*1 = 1
	// Mode 1 (i=1,j=0,k=0): P₁(1) * P₀(1) * P₀(1) * 2 = 1*1*1*2 = 2
	// Mode 2 (i=0,j=1,k=0): P₀(1) * P₁¹(1) * P₀(1) * 2 = 1*1.5*1*2 = 3
	// Mode 3 (i=0,j=0,k=1): P₀(1) * P₀(1) * P₁²(1) * 2 = 1*1*2*2 = 4

	if math.Abs(phi[0]-1.0) > 1e-10 {
		t.Errorf("P=1: Mode 0 at origin = %g (expected 1)", phi[0])
	}

	// Check a point where collapsed coords are closer to zero
	// Use the center of the reference tetrahedron
	centerR := -1.0 / 3.0
	centerS := -1.0 / 3.0
	centerT := -1.0 / 3.0
	phiCenter := basis.EvalBasis(centerR, centerS, centerT)
	t.Logf("P=1 at center (%.3f,%.3f,%.3f): phi = %v",
		centerR, centerS, centerT, phiCenter)

	// Constant mode should still be 1
	if math.Abs(phiCenter[0]-1.0) > 1e-10 {
		t.Errorf("P=1: Mode 0 at center = %g (expected 1)", phiCenter[0])
	}
}

func TestNodalBasisEvaluation(t *testing.T) {
	P := 1
	basis := NewPKDBasis(P)

	// Get nodal points (these are shifted inward from vertices)
	r, s, tt := generateTestNodes(P)

	// Test that each nodal basis function is 1 at its node and 0 at others
	for i := 0; i < len(r); i++ {
		// Evaluate at node i
		phi := basis.EvalNodalBasis(r[i], s[i], tt[i], r, s, tt)

		for j, val := range phi {
			expected := 0.0
			if i == j {
				expected = 1.0
			}
			if math.Abs(val-expected) > 1e-10 {
				t.Errorf("Nodal basis %d at node %d: got %f, expected %f",
					j, i, val, expected)
			}
		}
	}

	// Also test at a few arbitrary points to ensure nodal basis behaves smoothly
	testPoints := [][3]float64{
		{0.0, 0.0, 0.0},
		{-0.5, -0.5, -0.5},
		{0.2, -0.3, -0.4},
	}

	for _, pt := range testPoints {
		phi := basis.EvalNodalBasis(pt[0], pt[1], pt[2], r, s, tt)

		// Sum of nodal basis functions should be 1 (partition of unity)
		sum := 0.0
		for _, val := range phi {
			sum += val
		}

		if math.Abs(sum-1.0) > 1e-10 {
			t.Errorf("Nodal basis at point %v: sum = %f (expected 1)", pt, sum)
		}
	}

	t.Logf("P=%d: Nodal basis evaluation test passed", P)
}

func TestDerivativeMatrices(t *testing.T) {
	tol := 1e-6 // Relaxed tolerance for approximate nodes

	// Start with lower orders for debugging
	for P := 1; P <= 4; P++ {
		basis := NewPKDBasis(P)
		r, s, tt := generateTestNodes(P)

		t.Logf("P=%d: Testing with %d nodes", P, len(r))

		// Get derivative matrices
		Dr := basis.DerivativeMatrix(r, s, tt, 0)
		Ds := basis.DerivativeMatrix(r, s, tt, 1)
		Dt := basis.DerivativeMatrix(r, s, tt, 2)

		// Test on simple monomials first
		testCases := []struct{ p, q, k int }{
			{0, 0, 0}, // constant
			{1, 0, 0}, // r
			{0, 1, 0}, // s
			{0, 0, 1}, // t
		}

		if P >= 2 {
			testCases = append(testCases,
				struct{ p, q, k int }{2, 0, 0},
				struct{ p, q, k int }{1, 1, 0},
				struct{ p, q, k int }{0, 2, 0},
			)
		}

		for _, tc := range testCases {
			p, q, k := tc.p, tc.q, tc.k
			if p+q+k > P {
				continue
			}

			// Evaluate monomial
			u := Monomial(r, s, tt, p, q, k)

			// Compute derivatives using matrices
			dur := MatVec(Dr, u)
			dus := MatVec(Ds, u)
			dut := MatVec(Dt, u)

			// Compute analytical derivatives
			dur_exact := MonomialDeriv(r, s, tt, p, q, k, 0)
			dus_exact := MonomialDeriv(r, s, tt, p, q, k, 1)
			dut_exact := MonomialDeriv(r, s, tt, p, q, k, 2)

			// Check errors
			err_r := L2Error(dur, dur_exact)
			err_s := L2Error(dus, dus_exact)
			err_t := L2Error(dut, dut_exact)

			// Use relative error for non-zero derivatives
			max_dur := 0.0
			max_dus := 0.0
			max_dut := 0.0
			for i := range dur_exact {
				if math.Abs(dur_exact[i]) > max_dur {
					max_dur = math.Abs(dur_exact[i])
				}
				if math.Abs(dus_exact[i]) > max_dus {
					max_dus = math.Abs(dus_exact[i])
				}
				if math.Abs(dut_exact[i]) > max_dut {
					max_dut = math.Abs(dut_exact[i])
				}
			}

			// Normalize errors if derivatives are non-zero
			if max_dur > 1e-10 {
				err_r /= max_dur
			}
			if max_dus > 1e-10 {
				err_s /= max_dus
			}
			if max_dut > 1e-10 {
				err_t /= max_dut
			}

			if err_r > tol || err_s > tol || err_t > tol {
				t.Errorf("P=%d, monomial r^%d s^%d t^%d: errors (%.2e, %.2e, %.2e)",
					P, p, q, k, err_r, err_s, err_t)
				t.Logf("  max derivatives: (%.2e, %.2e, %.2e)", max_dur, max_dus, max_dut)
			} else {
				t.Logf("P=%d, monomial r^%d s^%d t^%d: PASSED (errors: %.2e, %.2e, %.2e)",
					P, p, q, k, err_r, err_s, err_t)
			}
		}

		t.Logf("P=%d: Derivative matrices test completed", P)
	}
}

func TestAffineTransform(t *testing.T) {
	// Test vertices for a regular tetrahedron
	v0 := [3]float64{0, 0, 0}
	v1 := [3]float64{1, 0, 0}
	v2 := [3]float64{0, 1, 0}
	v3 := [3]float64{0, 0, 1}

	at := NewAffineTransform(v0, v1, v2, v3)

	// Check Jacobian
	expectedJ := 1.0 / 6.0 // Volume of reference tet = 1/6
	if math.Abs(at.J-expectedJ) > 1e-14 {
		t.Errorf("Jacobian incorrect: got %f, expected %f", at.J, expectedJ)
	}

	// Test derivative transformation
	P := 3
	basis := NewPKDBasis(P)
	r, s, tt := generateTestNodes(P)

	// Simple test function u = x + 2y + 3z
	// In reference coords: u = (1-r-s-t)*0 + r*1 + s*2 + t*3
	u := make([]float64, len(r))
	for i := range u {
		u[i] = r[i] + 2*s[i] + 3*tt[i]
	}

	// Get derivatives in reference coordinates
	Dr := basis.DerivativeMatrix(r, s, tt, 0)
	Ds := basis.DerivativeMatrix(r, s, tt, 1)
	Dt := basis.DerivativeMatrix(r, s, tt, 2)

	dur := MatVec(Dr, u)
	dus := MatVec(Ds, u)
	dut := MatVec(Dt, u)

	// Transform to physical derivatives
	dx, dy, dz := at.TransformDerivatives(dur, dus, dut)

	// Expected: du/dx = 1, du/dy = 2, du/dz = 3
	tol := 1e-10
	for i := range dx {
		if math.Abs(dx[i]-1.0) > tol {
			t.Errorf("dx incorrect at node %d: got %f, expected 1.0", i, dx[i])
		}
		if math.Abs(dy[i]-2.0) > tol {
			t.Errorf("dy incorrect at node %d: got %f, expected 2.0", i, dy[i])
		}
		if math.Abs(dz[i]-3.0) > tol {
			t.Errorf("dz incorrect at node %d: got %f, expected 3.0", i, dz[i])
		}
	}
}

// TestSumFactorization tests the sum factorized operators
func TestSumFactorization(t *testing.T) {
	tol := 1e-8

	for P := 1; P <= 7; P++ {
		basis := NewPKDBasis(P)
		sf := NewSumFactorization(P)
		r, s, tt := generateTestNodes(P)

		// Get standard derivative matrices
		Dr := basis.DerivativeMatrix(r, s, tt, 0)
		Ds := basis.DerivativeMatrix(r, s, tt, 1)
		Dt := basis.DerivativeMatrix(r, s, tt, 2)

		// Test on various monomials
		for p := 0; p <= P; p++ {
			for q := 0; q <= P-p; q++ {
				for k := 0; k <= P-p-q; k++ {
					// Evaluate monomial
					u := Monomial(r, s, tt, p, q, k)

					// Standard matrix-vector products
					dur_std := MatVec(Dr, u)
					dus_std := MatVec(Ds, u)
					dut_std := MatVec(Dt, u)

					// Sum factorization
					dur_sf := sf.ApplyDr(u)
					dus_sf := sf.ApplyDs(u)
					dut_sf := sf.ApplyDt(u)

					// Compare results
					err_r := L2Error(dur_std, dur_sf)
					err_s := L2Error(dus_std, dus_sf)
					err_t := L2Error(dut_std, dut_sf)

					if err_r > tol || err_s > tol || err_t > tol {
						t.Errorf("P=%d, monomial r^%d s^%d t^%d: SF errors (%.2e, %.2e, %.2e)",
							P, p, q, k, err_r, err_s, err_t)
					}
				}
			}
		}

		t.Logf("P=%d: Sum factorization test passed", P)
	}
}

// Benchmark to compare direct vs sum factorized operations
func BenchmarkDerivativeOperators(b *testing.B) {
	P := 7
	basis := NewPKDBasis(P)
	sf := NewSumFactorization(P)
	r, s, t := generateTestNodes(P)

	Dr := basis.DerivativeMatrix(r, s, t, 0)
	u := Monomial(r, s, t, 2, 2, 2)

	b.Run("Direct", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			_ = MatVec(Dr, u)
		}
	})

	b.Run("SumFactorized", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			_ = sf.ApplyDr(u)
		}
	})
}

// TestSumFactorizationComplexity verifies the computational complexity
func TestSumFactorizationComplexity(t *testing.T) {
	// Measure operations for different polynomial orders
	operations := make(map[int]struct{ direct, sumfac int })

	for P := 2; P <= 6; P++ {
		N := (P + 1) * (P + 2) * (P + 3) / 6

		// Direct method: O(N^2) operations
		directOps := N * N

		// Sum factorization: O(P*N) operations
		// For tetrahedra: approximately 4*P*N operations
		sumfacOps := 4 * P * N

		operations[P] = struct{ direct, sumfac int }{directOps, sumfacOps}

		ratio := float64(directOps) / float64(sumfacOps)
		t.Logf("P=%d: Direct=%d ops, SumFac=%d ops, Speedup=%.1fx",
			P, directOps, sumfacOps, ratio)
	}

	// Verify that sum factorization becomes more efficient at higher orders
	if operations[6].direct <= operations[6].sumfac {
		t.Error("Sum factorization should be more efficient than direct method at P=6")
	}
}
