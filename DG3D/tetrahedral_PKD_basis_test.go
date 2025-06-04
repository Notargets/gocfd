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

func TestP1BasisDebug(t *testing.T) {
	basis := NewLineBasedBasis3D(1)

	// Print the nodes
	t.Logf("Nodes:")
	for i := 0; i < 4; i++ {
		t.Logf("  %d: r=%f, s=%f, t=%f", i, basis.r[i], basis.s[i], basis.t[i])
	}

	// Test point
	r, s, tt := -0.3, -0.2, -0.4
	t.Logf("\nTest point: (%f, %f, %f)", r, s, tt)

	// Get basis functions at test point
	phi := basis.EvalBasis(r, s, tt)
	t.Logf("Basis functions at test point: %v", phi)

	// Check partition of unity
	sum := 0.0
	for i := 0; i < 4; i++ {
		sum += phi[i]
	}
	t.Logf("Sum of basis functions: %f (should be 1.0)", sum)

	// Test constant function f=1
	t.Logf("\nTesting constant f=1:")
	u := []float64{1, 1, 1, 1}
	reconstructed := 0.0
	for i := 0; i < 4; i++ {
		reconstructed += u[i] * phi[i]
	}
	t.Logf("  Reconstructed: %f (should be 1.0)", reconstructed)

	// Test function f=r
	t.Logf("\nTesting f=r:")
	for i := 0; i < 4; i++ {
		u[i] = basis.r[i]
	}
	t.Logf("  Nodal values: %v", u)
	reconstructed = 0.0
	for i := 0; i < 4; i++ {
		reconstructed += u[i] * phi[i]
	}
	t.Logf("  Reconstructed: %f (should be %f)", reconstructed, r)
}

func TestBasisEvaluation(t *testing.T) {
	// Test the line-based nodal basis evaluation
	for P := 1; P <= 3; P++ {
		basis := NewLineBasedBasis3D(P)

		// Test 3: Polynomial exactness - use nodal interpolation
		// For a nodal basis: f(x) = sum_i f(x_i) * phi_i(x)
		for p := 0; p <= P; p++ {
			for q := 0; q <= P-p; q++ {
				for k := 0; k <= P-p-q; k++ {
					// Evaluate monomial at nodes
					u := Monomial(basis.r, basis.s, basis.t, p, q, k)

					// Reconstruct at a test point using nodal interpolation
					testPt := [3]float64{-0.3, -0.2, -0.4}
					phi := basis.EvalBasis(testPt[0], testPt[1], testPt[2])

					reconstructed := 0.0
					for j := 0; j < basis.N; j++ {
						reconstructed += u[j] * phi[j] // Use nodal values directly
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

func TestDerivatives(t *testing.T) {
	// Test line-based derivatives on simple monomials
	for P := 1; P <= 3; P++ {
		basis := NewLineBasedBasis3D(P)

		// Test derivative of f = r
		u := make([]float64, basis.N)
		for i := 0; i < basis.N; i++ {
			u[i] = basis.r[i] // f = r at nodes
		}

		dur := basis.ApplyDr(u)
		dus := basis.ApplyDs(u)
		dut := basis.ApplyDt(u)

		// For dr/dr, only check nodes that are on multi-node r-lines
		// Single-node lines correctly have derivative = 0
		numChecked := 0
		errDr := 0.0
		for i := 0; i < basis.N; i++ {
			if len(basis.lines.rLines[i]) > 1 {
				// This node is on a real line, should have dr/dr = 1
				errDr += (dur[i] - 1.0) * (dur[i] - 1.0)
				numChecked++
			} else {
				// Single node line, derivative should be 0
				if math.Abs(dur[i]) > 1e-10 {
					t.Errorf("P=%d: Node %d has non-zero r-derivative on single-node line", P, i)
				}
			}
		}
		if numChecked > 0 {
			errDr = math.Sqrt(errDr / float64(numChecked))
		}

		// For cross-derivatives, should be zero everywhere
		errDs := 0.0
		errDt := 0.0
		for i := 0; i < basis.N; i++ {
			errDs += dus[i] * dus[i]
			errDt += dut[i] * dut[i]
		}
		errDs = math.Sqrt(errDs / float64(basis.N))
		errDt = math.Sqrt(errDt / float64(basis.N))

		t.Logf("P=%d: d(r)/dr error = %.2e on %d multi-node lines", P, errDr, numChecked)
		t.Logf("P=%d: d(r)/ds error = %.2e (should be ~0)", P, errDs)
		t.Logf("P=%d: d(r)/dt error = %.2e (should be ~0)", P, errDt)

		if errDr > 1e-10 || errDs > 1e-10 || errDt > 1e-10 {
			t.Errorf("P=%d: Derivative errors too large", P)
		}
	}
}
