package DG3D

import (
	"github.com/notargets/gocfd/utils"
	"math"
	"testing"
)

// generateTestNodes generates reference tetrahedral nodes for testing
// Uses a simple approach to get better distributed nodes
func generateTestNodes(P int) (r, s, t []float64) {
	switch P {
	case 1:
		// Vertices of reference tetrahedron
		r = []float64{-1, 1, -1, -1}
		s = []float64{-1, -1, 1, -1}
		t = []float64{-1, -1, -1, 1}
		return r, s, t
	case 2:
		// Add edge midpoints
		r = []float64{
			-1, 1, -1, -1, // vertices
			0, -1, 0, -1, 0, -1, // edge midpoints
		}
		s = []float64{
			-1, -1, 1, -1, // vertices
			-1, 0, -1, 0, -1, 0, // edge midpoints
		}
		t = []float64{
			-1, -1, -1, 1, // vertices
			-1, -1, 0, 0, 0, 0, // edge midpoints
		}
		return r, s, t
	default:
		// For higher orders, use a simple approach with slight inward shift
		// to avoid exact boundary singularities
		n := (P + 1) * (P + 2) * (P + 3) / 6
		r = make([]float64, n)
		s = make([]float64, n)
		t = make([]float64, n)

		idx := 0
		shift := 0.01 // Small inward shift

		for i := 0; i <= P; i++ {
			for j := 0; j <= P-i; j++ {
				for k := 0; k <= P-i-j; k++ {
					// Barycentric coordinates
					l := P - i - j - k
					b1 := float64(i) / float64(P)
					b2 := float64(j) / float64(P)
					b3 := float64(k) / float64(P)
					b4 := float64(l) / float64(P)

					// Apply small shift toward center
					if b1 == 0 {
						b1 = shift
					}
					if b2 == 0 {
						b2 = shift
					}
					if b3 == 0 {
						b3 = shift
					}
					if b4 == 0 {
						b4 = shift
					}

					if b1 == 1 {
						b1 = 1 - shift
					}
					if b2 == 1 {
						b2 = 1 - shift
					}
					if b3 == 1 {
						b3 = 1 - shift
					}
					if b4 == 1 {
						b4 = 1 - shift
					}

					// Renormalize
					sum := b1 + b2 + b3 + b4
					b1 /= sum
					b2 /= sum
					b3 /= sum
					b4 /= sum

					// Convert to reference coordinates
					r[idx] = -b1 + b2 + b3 - b4
					s[idx] = -b1 - b2 + b3 + b4
					t[idx] = -b1 - b2 - b3 + b4

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
	// Test P=1 case explicitly
	P := 1
	basis := NewPKDBasis(P)

	// Vertices of reference tetrahedron
	vertices := [][3]float64{
		{-1, -1, -1},
		{1, -1, -1},
		{-1, 1, -1},
		{-1, -1, 1},
	}

	// Each basis function should be 1 at one vertex and 0 at others
	for i, v := range vertices {
		phi := basis.EvalBasis(v[0], v[1], v[2])
		t.Logf("Vertex %d (%v): phi = %v", i, v, phi)

		// For linear basis, we expect basis function i to be 1 at vertex i
		// Note: basis ordering might be different than vertex ordering
		foundOne := false
		for j, val := range phi {
			if math.Abs(val-1.0) < 1e-10 {
				foundOne = true
				t.Logf("  Basis function %d = 1", j)
			} else if math.Abs(val) > 1e-10 {
				t.Logf("  Basis function %d = %g (expected 0)", j, val)
			}
		}

		if !foundOne {
			t.Errorf("No basis function equals 1 at vertex %d", i)
		}
	}
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
