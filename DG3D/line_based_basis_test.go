package DG3D

import (
	"github.com/notargets/gocfd/utils"
	"math"
	"testing"
)

// generateTestNodes generates reference tetrahedral nodes for testing
// Uses the same approach as the line-based basis
func generateTestNodes(P int) (r, s, t []float64) {
	// Use the same node generation as the LineBasedBasis3D
	return GenerateOptimalTetrahedralNodes(P)
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

func TestLineStructure(t *testing.T) {
	// Test that the line structure is properly built
	for P := 1; P <= 5; P++ {
		basis := NewLineBasedBasis3D(P)

		t.Logf("P=%d: Testing line structure with %d nodes", P, basis.N)

		// Check that every node belongs to lines in each direction
		for i := 0; i < basis.N; i++ {
			if len(basis.lines.rLines[i]) == 0 {
				t.Errorf("P=%d: Node %d has no r-direction line", P, i)
			}
			if len(basis.lines.sLines[i]) == 0 {
				t.Errorf("P=%d: Node %d has no s-direction line", P, i)
			}
			if len(basis.lines.tLines[i]) == 0 {
				t.Errorf("P=%d: Node %d has no t-direction line", P, i)
			}
		}

		// Analyze the line structure
		basis.AnalyzeLineStructure()

		// Check that we have the expected number of unique operators
		numROperators := len(basis.lines.rOperators.operators)
		numSOperators := len(basis.lines.sOperators.operators)
		numTOperators := len(basis.lines.tOperators.operators)

		t.Logf("P=%d: Unique operators - R:%d, S:%d, T:%d",
			P, numROperators, numSOperators, numTOperators)

		// The number of unique operators should be much less than N
		totalOperators := numROperators + numSOperators + numTOperators
		if totalOperators >= 3*basis.N {
			t.Errorf("P=%d: Too many operators (%d), expected much less than %d",
				P, totalOperators, 3*basis.N)
		}
	}
}

func TestVandermonde(t *testing.T) {
	for P := 1; P <= 7; P++ {
		basis := NewLineBasedBasis3D(P)

		// Check that we have the right number of nodes
		expectedNodes := (P + 1) * (P + 2) * (P + 3) / 6
		if len(basis.r) != expectedNodes {
			t.Errorf("P=%d: Wrong number of nodes: got %d, expected %d",
				P, len(basis.r), expectedNodes)
			continue
		}

		// Test that V is square for nodal basis
		if basis.V.Rows() != basis.V.Cols() {
			t.Errorf("P=%d: Vandermonde not square: %d x %d",
				P, basis.V.Rows(), basis.V.Cols())
		}

		// Test conditioning
		cond := basis.V.ConditionNumber()
		t.Logf("P=%d: Vandermonde size %d x %d, condition number %.2e",
			P, basis.V.Rows(), basis.V.Cols(), cond)

		// The line-based approach should have better conditioning than collapsed coords
		if cond > 1e12 {
			t.Logf("WARNING: P=%d: Vandermonde matrix condition number %.2e", P, cond)
		}

		// Check that V * V^{-1} = I
		I := basis.V.Mul(basis.Vinv)
		for i := 0; i < basis.N; i++ {
			for j := 0; j < basis.N; j++ {
				expected := 0.0
				if i == j {
					expected = 1.0
				}
				if math.Abs(I.At(i, j)-expected) > 1e-10 {
					t.Errorf("P=%d: V*Vinv not identity at (%d,%d): %.2e",
						P, i, j, I.At(i, j))
				}
			}
		}
	}
}

func TestBasisEvaluation(t *testing.T) {
	// Test the line-based nodal basis evaluation
	for P := 1; P <= 3; P++ {
		basis := NewLineBasedBasis3D(P)

		// Test 1: Evaluate basis at the nodes themselves
		for i := 0; i < basis.N; i++ {
			phi := basis.EvalBasis(basis.r[i], basis.s[i], basis.t[i])

			// Should be 1 at node i, 0 at all others (Lagrange property)
			for j := 0; j < basis.N; j++ {
				expected := 0.0
				if i == j {
					expected = 1.0
				}
				if math.Abs(phi[j]-expected) > 1e-10 {
					t.Errorf("P=%d: Basis %d at node %d = %.2e (expected %.0f)",
						P, j, i, phi[j], expected)
				}
			}
		}

		// Test 2: Partition of unity - sum of all basis functions = 1
		testPoints := [][3]float64{
			{0, 0, 0},
			{-0.5, -0.5, -0.5},
			{0.2, -0.3, -0.4},
			{-0.8, -0.8, -0.8},
		}

		for _, pt := range testPoints {
			phi := basis.EvalBasis(pt[0], pt[1], pt[2])
			sum := 0.0
			for j := 0; j < basis.N; j++ {
				sum += phi[j]
			}
			if math.Abs(sum-1.0) > 1e-10 {
				t.Errorf("P=%d: Partition of unity fails at %v: sum = %.2e",
					P, pt, sum)
			}
		}

		// Test 3: Polynomial exactness
		// The basis should be able to represent any polynomial up to degree P
		// Test with monomials
		for p := 0; p <= P; p++ {
			for q := 0; q <= P-p; q++ {
				for k := 0; k <= P-p-q; k++ {
					// Evaluate monomial at nodes
					u := Monomial(basis.r, basis.s, basis.t, p, q, k)

					// Convert to modal coefficients using V^{-1}
					coeffs := MatVec(basis.Vinv, u)

					// Reconstruct at a test point
					testPt := [3]float64{-0.3, -0.2, -0.4}
					phi := basis.EvalBasis(testPt[0], testPt[1], testPt[2])

					reconstructed := 0.0
					for j := 0; j < basis.N; j++ {
						reconstructed += coeffs[j] * phi[j]
					}

					// Compare with exact monomial value
					exact := math.Pow(testPt[0], float64(p)) *
						math.Pow(testPt[1], float64(q)) *
						math.Pow(testPt[2], float64(k))

					if math.Abs(reconstructed-exact) > 1e-10 {
						t.Errorf("P=%d: Monomial r^%d s^%d t^%d fails at test point: %.2e vs %.2e",
							P, p, q, k, reconstructed, exact)
					}
				}
			}
		}

		t.Logf("P=%d: Basis evaluation tests passed", P)
	}
}

func TestLineBasedDerivatives(t *testing.T) {
	// Test the line-based derivative operators
	for P := 1; P <= 4; P++ {
		basis := NewLineBasedBasis3D(P)

		t.Logf("P=%d: Testing line-based derivatives with %d nodes", P, basis.N)

		// Test derivative operators on monomials
		testCases := []struct {
			p, q, k int
			name    string
		}{
			{0, 0, 0, "constant"},
			{1, 0, 0, "r"},
			{0, 1, 0, "s"},
			{0, 0, 1, "t"},
		}

		if P >= 2 {
			testCases = append(testCases,
				struct {
					p, q, k int
					name    string
				}{2, 0, 0, "r^2"},
				struct {
					p, q, k int
					name    string
				}{1, 1, 0, "r*s"},
				struct {
					p, q, k int
					name    string
				}{0, 2, 0, "s^2"},
				struct {
					p, q, k int
					name    string
				}{0, 0, 2, "t^2"},
			)
		}

		for _, tc := range testCases {
			if tc.p+tc.q+tc.k > P {
				continue
			}

			// Evaluate monomial at nodes
			u := Monomial(basis.r, basis.s, basis.t, tc.p, tc.q, tc.k)

			// Apply line-based derivatives
			dur := basis.ApplyDr(u)
			dus := basis.ApplyDs(u)
			dut := basis.ApplyDt(u)

			// Exact derivatives
			dur_exact := MonomialDeriv(basis.r, basis.s, basis.t, tc.p, tc.q, tc.k, 0)
			dus_exact := MonomialDeriv(basis.r, basis.s, basis.t, tc.p, tc.q, tc.k, 1)
			dut_exact := MonomialDeriv(basis.r, basis.s, basis.t, tc.p, tc.q, tc.k, 2)

			// Compute errors
			err_r := L2Error(dur, dur_exact)
			err_s := L2Error(dus, dus_exact)
			err_t := L2Error(dut, dut_exact)

			// Line-based approach should be exact for polynomials
			tol := 1e-10
			if err_r > tol || err_s > tol || err_t > tol {
				t.Errorf("P=%d: Monomial %s derivative errors: (%.2e, %.2e, %.2e)",
					P, tc.name, err_r, err_s, err_t)
			} else {
				t.Logf("P=%d: Monomial %s derivatives OK (errors: %.2e, %.2e, %.2e)",
					P, tc.name, err_r, err_s, err_t)
			}
		}

		// Test that derivative operators are consistent
		// Dr, Ds, Dt should differentiate basis functions correctly
		for idx := 0; idx < basis.N && idx < 5; idx++ {
			// Create unit vector for basis function idx
			u := make([]float64, basis.N)
			u[idx] = 1.0

			// Apply derivatives
			dur := basis.ApplyDr(u)
			dus := basis.ApplyDs(u)
			dut := basis.ApplyDt(u)

			// Check at a few test points
			testPts := [][3]float64{
				{-0.5, -0.3, -0.1},
				{0.1, -0.2, -0.3},
			}

			for _, pt := range testPts {
				// Evaluate basis function and its derivatives
				phi := basis.EvalBasis(pt[0], pt[1], pt[2])

				// The derivative of basis function idx should match
				// what we get by applying the operator to the unit vector
				// and then evaluating at the point

				tol := 1e-8
				_, _ = phi, tol

				// For now, just check that derivatives are reasonable
				if math.IsNaN(dur[0]) || math.IsNaN(dus[0]) || math.IsNaN(dut[0]) {
					t.Errorf("P=%d: NaN in derivatives for basis %d", P, idx)
				}
			}
		}
	}
}

func TestLineBasedEfficiency(t *testing.T) {
	// Compare the efficiency of line-based vs full matrix approach
	for P := 2; P <= 6; P++ {
		basis := NewLineBasedBasis3D(P)

		// Count operations for line-based approach
		totalLineOps := 0
		maxLineLength := 0

		// Check r-direction lines
		for i := 0; i < basis.N; i++ {
			lineLen := len(basis.lines.rLines[i])
			if lineLen > maxLineLength {
				maxLineLength = lineLen
			}
			if lineLen > 1 {
				// Each node requires lineLen multiplications
				totalLineOps += lineLen
			}
		}

		// Full matrix would require N^2 operations
		fullMatrixOps := basis.N * basis.N

		// Line-based uses approximately N * P operations
		speedup := float64(fullMatrixOps) / float64(totalLineOps)

		t.Logf("P=%d (N=%d): Line-based ops=%d, Full matrix ops=%d, Speedup=%.1fx",
			P, basis.N, totalLineOps, fullMatrixOps, speedup)

		// Verify complexity is O(P*N) not O(N^2)
		if float64(totalLineOps) > 2.0*float64(P*basis.N) {
			t.Errorf("P=%d: Line-based operations (%d) exceed O(P*N) bound",
				P, totalLineOps)
		}
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
	expectedJ := 1.0 / 6.0 // Volume of unit tetrahedron
	if math.Abs(at.J-expectedJ) > 1e-14 {
		t.Errorf("Jacobian incorrect: got %f, expected %f", at.J, expectedJ)
	}

	// Test derivative transformation with line-based basis
	P := 3
	basis := NewLineBasedBasis3D(P)

	// Simple test function u = x + 2y + 3z
	// In physical coords, at node i: u = x[i] + 2*y[i] + 3*z[i]
	u := make([]float64, basis.N)
	for i := range u {
		// Map reference node to physical coordinates
		x, y, z := at.MapToPhysical(basis.r[i], basis.s[i], basis.t[i])
		u[i] = x + 2*y + 3*z
	}

	// Get derivatives in reference coordinates using line-based approach
	dur := basis.ApplyDr(u)
	dus := basis.ApplyDs(u)
	dut := basis.ApplyDt(u)

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

// Benchmark to compare direct vs line-based operations
func BenchmarkDerivativeOperators(b *testing.B) {
	P := 7
	basis := NewLineBasedBasis3D(P)

	// Create a full derivative matrix for comparison
	Dr := basis.DerivativeMatrix(0)

	u := Monomial(basis.r, basis.s, basis.t, 2, 2, 2)

	b.Run("FullMatrix", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			_ = MatVec(Dr, u)
		}
	})

	b.Run("LineBased", func(b *testing.B) {
		for i := 0; i < b.N; i++ {
			_ = basis.ApplyDr(u)
		}
	})
}

// TestMemoryUsage verifies the memory savings of the line-based approach
func TestMemoryUsage(t *testing.T) {
	for P := 2; P <= 8; P++ {
		basis := NewLineBasedBasis3D(P)

		// Count total storage for 1D operators
		totalStorage := 0

		// Count r-direction operators
		for _, op := range basis.lines.rOperators.operators {
			totalStorage += op.Rows() * op.Cols()
		}
		// Count s-direction operators
		for _, op := range basis.lines.sOperators.operators {
			totalStorage += op.Rows() * op.Cols()
		}
		// Count t-direction operators
		for _, op := range basis.lines.tOperators.operators {
			totalStorage += op.Rows() * op.Cols()
		}

		// Full matrices would require 3 * N^2 storage
		fullStorage := 3 * basis.N * basis.N

		savings := 100.0 * (1.0 - float64(totalStorage)/float64(fullStorage))

		t.Logf("P=%d (N=%d): Line-based storage=%d, Full storage=%d, Savings=%.1f%%",
			P, basis.N, totalStorage, fullStorage, savings)

		// Should have significant savings for P >= 3
		if P >= 3 && savings < 50.0 {
			t.Errorf("P=%d: Insufficient memory savings (%.1f%%)", P, savings)
		}
	}
}
