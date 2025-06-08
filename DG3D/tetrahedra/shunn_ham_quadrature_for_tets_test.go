package tetrahedra

import (
	"fmt"
	"math"
	"testing"
)

// referenceTetrahedronMonomialIntegral computes the analytical integral of r^p * s^q * t^m
// over the reference tetrahedron with vertices at (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
func referenceTetrahedronMonomialIntegral(p, q, m int) float64 {
	// The transformation from barycentric to reference coordinates is:
	// r = -1 + 2(ξ₁ + ξ₂ + ξ₃)
	// s = -1 + 2(ξ₂ + ξ₃)
	// t = -1 + 2ξ₃
	// Jacobian = 8

	integral := 0.0

	// Expand r^p = (-1 + 2(ξ₁ + ξ₂ + ξ₃))^p using binomial theorem
	for i := 0; i <= p; i++ {
		coeffR := float64(binomial(p, i)) * math.Pow(-1, float64(p-i)) * math.Pow(2, float64(i))

		// Expand s^q = (-1 + 2(ξ₂ + ξ₃))^q
		for j := 0; j <= q; j++ {
			coeffS := float64(binomial(q, j)) * math.Pow(-1, float64(q-j)) * math.Pow(2, float64(j))

			// Expand t^m = (-1 + 2ξ₃)^m
			for k := 0; k <= m; k++ {
				coeffT := float64(binomial(m, k)) * math.Pow(-1, float64(m-k)) * math.Pow(2, float64(k))

				// Now expand (ξ₁ + ξ₂ + ξ₃)^i using multinomial theorem
				for a1 := 0; a1 <= i; a1++ {
					for a2 := 0; a2 <= i-a1; a2++ {
						a3 := i - a1 - a2
						coeffI := float64(multinomial(i, []int{a1, a2, a3}))

						// Expand (ξ₂ + ξ₃)^j using binomial theorem
						for b2 := 0; b2 <= j; b2++ {
							b3 := j - b2
							coeffJ := float64(binomial(j, b2))

							// Total powers of ξ₁, ξ₂, ξ₃
							pow1 := a1
							pow2 := a2 + b2
							pow3 := a3 + b3 + k

							// Integral of ξ₁^pow1 * ξ₂^pow2 * ξ₃^pow3 over unit simplex
							simplexIntegral := float64(factorial(pow1)*factorial(pow2)*factorial(pow3)) /
								float64(factorial(pow1+pow2+pow3+3))

							// Add contribution with Jacobian = 8
							integral += 8.0 * coeffR * coeffS * coeffT * coeffI * coeffJ * simplexIntegral
						}
					}
				}
			}
		}
	}

	// Clean up near-zero values
	if math.Abs(integral) < 1e-14 {
		return 0.0
	}

	return integral
}

// Helper functions

func factorial(n int) int {
	if n <= 1 {
		return 1
	}
	result := 1
	for i := 2; i <= n; i++ {
		result *= i
	}
	return result
}

// TestTetrahedralQuadratureMonomialIntegration tests monomial integration on both unit simplex and reference tet
func TestTetrahedralQuadratureMonomialIntegration(t *testing.T) {
	testCases := []struct {
		order     int
		maxDegree int
	}{
		{1, 1},
		{2, 2},
		{3, 3},
		{4, 5},
		{5, 6},
		{6, 8},
	}

	// Test both unit simplex and reference tetrahedron
	for _, tetType := range []string{"unit", "reference"} {
		for _, tc := range testCases {
			quad, err := NewShunnHamQuadrature(tc.order)
			if err != nil {
				t.Fatalf("Failed to create quadrature order %d: %v", tc.order, err)
			}

			// Get quadrature points and weights based on tet type
			var r, s, tt, w []float64
			if tetType == "unit" {
				r, s, tt, w = quad.GetRSTWUnit()
			} else {
				r, s, tt, w = quad.GetRSTWReference()
			}

			// Test all monomials up to maxDegree
			for totalDegree := 0; totalDegree <= tc.maxDegree; totalDegree++ {
				for p := 0; p <= totalDegree; p++ {
					for q := 0; q <= totalDegree-p; q++ {
						rr := totalDegree - p - q

						// Compute numerical integral using shared function
						numerical := computeIntegral(r, s, tt, w, p, q, rr)

						// Compute analytical integral based on tet type
						var analytical float64
						if tetType == "unit" {
							// For unit simplex: ∫∫∫ x^p * y^q * z^r dV = p! * q! * r! / (p+q+r+3)!
							analytical = float64(factorial(p)*factorial(q)*factorial(rr)) / float64(factorial(p+q+rr+3))
						} else {
							// For reference tet: use the transformation formula
							analytical = referenceTetrahedronMonomialIntegral(p, q, rr)
						}

						// Check accuracy with tolerance of 1e-12
						tolerance := 1e-12
						if math.Abs(numerical-analytical) > tolerance {
							t.Errorf("%s tet order %d: monomial x^%d*y^%d*z^%d integration failed. "+
								"Expected: %.15e, Got: %.15e, Diff: %.15e",
								tetType, tc.order, p, q, rr, analytical, numerical, math.Abs(numerical-analytical))
						}
					}
				}
			}
			t.Logf("%s tet order %d: successfully integrated all monomials up to degree %d",
				tetType, tc.order, tc.maxDegree)
		}
	}
}

// computeIntegral computes the numerical integral of x^p * y^q * z^r
func computeIntegral(r, s, t, w []float64, p, q, rr int) float64 {
	integral := 0.0
	for i := range r {
		integral += w[i] * math.Pow(r[i], float64(p)) *
			math.Pow(s[i], float64(q)) * math.Pow(t[i], float64(rr))
	}
	return integral
}

func binomial(n, k int) int {
	if k < 0 || k > n {
		return 0
	}
	if k == 0 || k == n {
		return 1
	}
	return factorial(n) / (factorial(k) * factorial(n-k))
}

func multinomial(n int, k []int) int {
	sum := 0
	for _, ki := range k {
		sum += ki
	}
	if sum != n {
		return 0
	}

	result := factorial(n)
	for _, ki := range k {
		result /= factorial(ki)
	}
	return result
}

// TestGetNodesShunnHam validates the nodes returned by GetNodesShunnHam
func TestGetNodesShunnHam(t *testing.T) {
	// Test each order from 1 to 6
	for P := 1; P <= 5; P++ {
		t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
			// Get nodes
			R, S, T := GetNodesShunnHam(P)

			// Check that all vectors have the same length
			if R.Len() != S.Len() || S.Len() != T.Len() {
				t.Errorf("N %d: Vector lengths don't match: R=%d, S=%d, T=%d",
					P, R.Len(), S.Len(), T.Len())
				return
			}

			// The number of points needed to interpolate a polynomial of degree P in a tetrahedron
			// is (P+1)(P+2)(P+3)/6
			expectedPoints := (P + 1) * (P + 2) * (P + 3) / 6

			if R.Len() != expectedPoints {
				t.Errorf("N %d: Expected %d points to interpolate polynomial degree %d, got %d",
					P, expectedPoints, P, R.Len())
			}

			t.Logf("N %d: %d points provided to interpolate polynomial degree %d",
				P, R.Len(), P)
		})
	}
}
