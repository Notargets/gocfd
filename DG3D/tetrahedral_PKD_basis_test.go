package DG3D

import (
	"fmt"
	"math"
	"math/rand/v2"
	"testing"

	"github.com/notargets/gocfd/utils"
)

// TestPKDBasisInterpolation verifies that the PKD basis can interpolate polynomials
// exactly up to order P
func TestPKDBasisInterpolation(t *testing.T) {
	var (
		tol = 1.e-12
	)
	for i, label := range []string{"Equispaced", "Shunn Ham"} {
		t.Logf("Using %s nodes\n", label)
		// Test orders 1 through 5
		for P := 1; P <= 5; P++ {
			t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
				// Create basis
				basis := NewTetBasis(P)

				// Get interpolation nodes
				var R, S, T utils.Vector
				switch i {
				case 0:
					R, S, T = EquispacedNodes3D(P)
				case 1:
					tol = 5.e-11
					R, S, T = GetNodesShunnHam(P)
				}

				// Compute Vandermonde matrix at the nodes
				V := basis.ComputeVandermonde(R, S, T)

				// Get the inverse for computing modal coefficients
				VInv := V.InverseWithCheck()

				// Define test polynomials with their total order
				testPolynomials := []struct {
					name     string
					f        func(r, s, t float64) float64
					maxOrder int // Maximum total order of polynomial
				}{
					// Order 0
					{"constant", func(r, s, t float64) float64 { return 1.0 }, 0},

					// Order 1
					{"linear_r", func(r, s, t float64) float64 { return r }, 1},
					{"linear_s", func(r, s, t float64) float64 { return s }, 1},
					{"linear_t", func(r, s, t float64) float64 { return t }, 1},
					{"linear_combo", func(r, s, t float64) float64 { return 2*r - 3*s + t }, 1},

					// Order 2
					{"quadratic_r2", func(r, s, t float64) float64 { return r * r }, 2},
					{"quadratic_s2", func(r, s, t float64) float64 { return s * s }, 2},
					{"quadratic_t2", func(r, s, t float64) float64 { return t * t }, 2},
					{"quadratic_rs", func(r, s, t float64) float64 { return r * s }, 2},
					{"quadratic_rt", func(r, s, t float64) float64 { return r * t }, 2},
					{"quadratic_st", func(r, s, t float64) float64 { return s * t }, 2},
					{"quadratic_combo", func(r, s, t float64) float64 { return r*r + 2*r*s - t*t }, 2},

					// Order 3
					{"cubic_r3", func(r, s, t float64) float64 { return r * r * r }, 3},
					{"cubic_s3", func(r, s, t float64) float64 { return s * s * s }, 3},
					{"cubic_rst", func(r, s, t float64) float64 { return r * s * t }, 3},
					{"cubic_r2s", func(r, s, t float64) float64 { return r * r * s }, 3},
					{"cubic_combo", func(r, s, t float64) float64 { return r*r*r - 2*r*s*t + s*s*t }, 3},

					// Order 4
					{"quartic_r4", func(r, s, t float64) float64 { return r * r * r * r }, 4},
					{"quartic_r2s2", func(r, s, t float64) float64 { return r * r * s * s }, 4},
					{"quartic_r3t", func(r, s, t float64) float64 { return r * r * r * t }, 4},
					{"quartic_combo", func(r, s, t float64) float64 { return r*r*r*r + r*r*s*s - 2*r*s*t*t }, 4},

					// Order 5
					{"quintic_r5", func(r, s, t float64) float64 { return r * r * r * r * r }, 5},
					{"quintic_r3st", func(r, s, t float64) float64 { return r * r * r * s * t }, 5},
					{"quintic_r2s2t", func(r, s, t float64) float64 { return r * r * s * s * t }, 5},
					{"quintic_combo", func(r, s, t float64) float64 { return r*r*r*r*r - 3*r*r*r*s*t + s*s*s*t*t }, 5},
				}

				for _, test := range testPolynomials {
					// Only test polynomials up to order P
					if test.maxOrder > P {
						continue
					}

					// Evaluate function at interpolation nodes
					fVals := make([]float64, R.Len())
					for i := 0; i < R.Len(); i++ {
						fVals[i] = test.f(R.At(i), S.At(i), T.At(i))
					}
					fVec := utils.NewVector(len(fVals), fVals)

					// Compute modal coefficients: coeffs = V^{-1} * f
					coeffs := VInv.Mul(fVec.ToMatrix())

					// Test interpolation at random points within the reference tetrahedron
					nTest := 50
					maxError := 0.0

					for iTest := 0; iTest < nTest; iTest++ {
						// Generate random point in reference tetrahedron using barycentric coordinates
						// L1, L2, L3, L4 >= 0 and L1 + L2 + L3 + L4 = 1
						L1 := rand.Float64()
						L2 := rand.Float64() * (1 - L1)
						L3 := rand.Float64() * (1 - L1 - L2)
						L4 := 1 - L1 - L2 - L3

						// Convert to reference coordinates
						r := -L1 - L2 - L3 + L4
						s := -L1 - L2 + L3 - L4
						t := -L1 + L2 - L3 - L4

						// Evaluate basis functions at test point
						phi := basis.EvaluateBasis(r, s, t)

						// Compute interpolated value
						interpVal := 0.0
						for i := 0; i < basis.Np; i++ {
							interpVal += coeffs.At(i, 0) * phi[i]
						}

						// Compare with exact value
						exactVal := test.f(r, s, t)
						err := math.Abs(interpVal - exactVal)
						if err > maxError {
							maxError = err
						}
					}

					// Check that interpolation is exact (within tolerance)
					// Use a tolerance that accounts for conditioning
					interpTol := 1e-10 * float64(P*P)
					if maxError > interpTol {
						t.Errorf("Order %d, polynomial %s: max interpolation error %g exceeds tolerance %g",
							P, test.name, maxError, interpTol)
					}
				}

				// Additional test: Verify that the basis can exactly represent all polynomials
				// up to order P by checking at the interpolation nodes themselves
				for _, test := range testPolynomials {
					if test.maxOrder > P {
						continue
					}

					// Evaluate function at nodes
					fVals := make([]float64, R.Len())
					for i := 0; i < R.Len(); i++ {
						fVals[i] = test.f(R.At(i), S.At(i), T.At(i))
					}
					fVec := utils.NewVector(len(fVals), fVals)

					// Compute modal coefficients and reconstruct at nodes
					coeffs := VInv.Mul(fVec.ToMatrix())
					fReconstructed := V.Mul(coeffs)

					// Check reconstruction error at nodes
					maxNodeError := tol
					for i := 0; i < R.Len(); i++ {
						err := math.Abs(fReconstructed.At(i, 0) - fVals[i])
						if err > maxNodeError {
							maxNodeError = err
						}
					}

					if maxNodeError > tol {
						t.Errorf("Order %d, polynomial %s: reconstruction error at nodes %g",
							P, test.name, maxNodeError)
					}
				}
			})
		}
	}
}
