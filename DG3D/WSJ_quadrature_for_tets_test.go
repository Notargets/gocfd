package DG3D

import (
	"fmt"
	"math"
	"testing"
)

// TestMonomialIntegration tests monomial integration for orders 1-5
func TestMonomialIntegration(t *testing.T) {
	// Test orders 1 through 5
	for order := 1; order <= 5; order++ {
		t.Run(fmt.Sprintf("Order %d", order), func(t *testing.T) {
			// Create quadrature rule
			quad, err := NewUnitSimplexQuadrature(order)
			if err != nil {
				t.Fatalf("Failed to create order %d quadrature: %v", order, err)
			}

			// Test monomials up to degree = order
			for degree := 0; degree <= order; degree++ {
				// Test r^degree
				t.Run(fmt.Sprintf("r^%d", degree), func(t *testing.T) {
					numerical := 0.0
					for _, pt := range quad.Points {
						f := math.Pow(pt.R, float64(degree))
						numerical += pt.W * f
					}
					numerical /= 6.0 // Scale by 1/volume

					// Exact integral of r^degree over unit simplex
					exact := factorial(degree) / factorial(degree+3)

					relError := math.Abs(numerical-exact) / exact
					if relError > 1e-13 {
						t.Errorf("Failed: numerical=%.16f, exact=%.16f, relError=%.3e",
							numerical, exact, relError)
					}
				})

				// Test s^degree
				t.Run(fmt.Sprintf("s^%d", degree), func(t *testing.T) {
					numerical := 0.0
					for _, pt := range quad.Points {
						f := math.Pow(pt.S, float64(degree))
						numerical += pt.W * f
					}
					numerical /= 6.0 // Scale by 1/volume

					// Exact integral of s^degree over unit simplex
					exact := factorial(degree) / factorial(degree+3)

					relError := math.Abs(numerical-exact) / exact
					if relError > 1e-13 {
						t.Errorf("Failed: numerical=%.16f, exact=%.16f, relError=%.3e",
							numerical, exact, relError)
					}
				})

				// Test t^degree
				t.Run(fmt.Sprintf("t^%d", degree), func(t *testing.T) {
					numerical := 0.0
					for _, pt := range quad.Points {
						f := math.Pow(pt.T, float64(degree))
						numerical += pt.W * f
					}
					numerical /= 6.0 // Scale by 1/volume

					// Exact integral of t^degree over unit simplex
					exact := factorial(degree) / factorial(degree+3)

					relError := math.Abs(numerical-exact) / exact
					if relError > 1e-13 {
						t.Errorf("Failed: numerical=%.16f, exact=%.16f, relError=%.3e",
							numerical, exact, relError)
					}
				})
			}
		})
	}
}

// factorial computes n!
func factorial(n int) float64 {
	if n <= 0 {
		return 1.0
	}
	result := 1.0
	for i := 2; i <= n; i++ {
		result *= float64(i)
	}
	return result
}
