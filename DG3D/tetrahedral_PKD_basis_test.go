package DG3D

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"
)

// TestPKDBasisInterpolation verifies that the PKD basis can interpolate polynomials
// exactly up to order P using quadrature integration
func _TestPKDBasisInterpolation(t *testing.T) {
	// Test orders 1 through 6
	for P := 1; P <= 5; P++ {
		t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
			// Create basis
			basis := NewTetBasis(P)

			// Get interpolation nodes using equispaced nodes for now
			// since GetNodesShunnHam returns quadrature points, not interpolation nodes
			R, S, T := EquispacedNodes3D(P)

			// Compute Vandermonde matrix at the nodes
			V := basis.ComputeVandermonde(R, S, T)

			// Get the inverse for later use
			VInv := V.InverseWithCheck()

			// Test 1: Check orthogonality - Mass matrix should be diagonal (not identity)
			// Get quadrature for exact integration
			quadOrder := P
			if 2*P > 6 {
				quadOrder = 6 // Max available
			}
			quad, err := NewShunnHamQuadrature(quadOrder)
			if err != nil {
				t.Fatalf("Failed to create quadrature: %v", err)
			}

			// Get quadrature points and weights for reference tetrahedron
			rq, sq, tq, wq := quad.GetRSTWReference()
			Rq := utils.NewVector(len(rq), rq)
			Sq := utils.NewVector(len(sq), sq)
			Tq := utils.NewVector(len(tq), tq)

			// Create cubature structure
			cubature := &TetCubature{
				R: Rq,
				S: Sq,
				T: Tq,
				W: utils.NewVector(len(wq), wq),
			}

			// Compute Vandermonde at quadrature points
			Vq := basis.ComputeVandermonde(Rq, Sq, Tq)

			// Compute mass matrix
			M := basis.ComputeMassMatrix(Vq, cubature)

			// Check that M is diagonal (for orthogonal basis)
			tol := 1e-10
			for i := 0; i < basis.Np; i++ {
				for j := 0; j < basis.Np; j++ {
					if i != j {
						// Off-diagonal elements should be zero
						if math.Abs(M.At(i, j)) > tol {
							t.Errorf("Order %d: Mass matrix not diagonal at (%d,%d): got %g, want 0",
								P, i, j, M.At(i, j))
						}
					} else {
						// Diagonal elements should be positive
						if M.At(i, i) <= 0 {
							t.Errorf("Order %d: Mass matrix diagonal at (%d,%d) not positive: got %g",
								P, i, i, M.At(i, i))
						}
					}
				}
			}

			// Test 2: Interpolate test polynomials and verify exactness
			testPolynomials := []struct {
				name     string
				f        func(r, s, t float64) float64
				maxOrder int // Maximum total order of polynomial
			}{
				{"constant", func(r, s, t float64) float64 { return 1.0 }, 0},
				{"linear_r", func(r, s, t float64) float64 { return r }, 1},
				{"linear_s", func(r, s, t float64) float64 { return s }, 1},
				{"linear_t", func(r, s, t float64) float64 { return t }, 1},
				{"quadratic_r2", func(r, s, t float64) float64 { return r * r }, 2},
				{"quadratic_rs", func(r, s, t float64) float64 { return r * s }, 2},
				{"quadratic_st", func(r, s, t float64) float64 { return s * t }, 2},
				{"cubic_r3", func(r, s, t float64) float64 { return r * r * r }, 3},
				{"cubic_rst", func(r, s, t float64) float64 { return r * s * t }, 3},
				{"quartic_r2s2", func(r, s, t float64) float64 { return r * r * s * s }, 4},
				{"quintic_r3st", func(r, s, t float64) float64 { return r * r * r * s * t }, 5},
				{"sextic_r2s2t2", func(r, s, t float64) float64 { return r * r * s * s * t * t }, 6},
			}

			for _, test := range testPolynomials {
				// Only test polynomials up to order P
				if test.maxOrder > P {
					continue
				}

				// Evaluate function at nodes
				fVals := make([]float64, R.Len())
				for i := 0; i < R.Len(); i++ {
					fVals[i] = test.f(R.At(i), S.At(i), T.At(i))
				}
				fVec := utils.NewVector(len(fVals), fVals)

				// Compute modal coefficients: coeffs = V^{-1} * f
				coeffs := VInv.Mul(fVec.ToMatrix())

				// Test at random points within the reference tetrahedron
				nTest := 20
				maxError := 0.0

				for iTest := 0; iTest < nTest; iTest++ {
					// Generate random point in reference tetrahedron
					// Use rejection sampling to ensure point is inside
					var r, s, t float64
					for {
						r = -1.0 + 2.0*rand.Float64()
						s = -1.0 + 2.0*rand.Float64()
						t = -1.0 + 2.0*rand.Float64()

						// Check if inside reference tetrahedron: r+s+t <= 1
						if r+s+t <= 1.0 {
							break
						}
					}

					// Evaluate basis at test point
					phi := basis.EvaluateBasis(r, s, t)

					// Compute interpolated value
					interpVal := 0.0
					for i := 0; i < basis.Np; i++ {
						interpVal += coeffs.At(i, 0) * phi[i]
					}

					// Compare with exact value
					exactVal := test.f(r, s, t)
					error := math.Abs(interpVal - exactVal)
					if error > maxError {
						maxError = error
					}
				}

				// Check that interpolation is exact (within tolerance)
				interpTol := 1e-10
				if maxError > interpTol {
					t.Errorf("Order %d, polynomial %s: max interpolation error %g exceeds tolerance %g",
						P, test.name, maxError, interpTol)
				}
			}

			// Test 3: Verify partition of unity
			// Sum of all basis functions should equal 1 everywhere
			nTest := 10
			for iTest := 0; iTest < nTest; iTest++ {
				// Generate random point in reference tetrahedron
				var rTest, sTest, tTest float64
				for {
					rTest = -1.0 + 2.0*rand.Float64()
					sTest = -1.0 + 2.0*rand.Float64()
					tTest = -1.0 + 2.0*rand.Float64()

					// Check if inside: r+s+t <= 1
					if rTest+sTest+tTest <= 1.0 {
						break
					}
				}

				// Evaluate all basis functions
				phi := basis.EvaluateBasis(rTest, sTest, tTest)

				// Compute coefficients for constant function f=1
				oneVals := make([]float64, R.Len())
				for i := 0; i < R.Len(); i++ {
					oneVals[i] = 1.0
				}
				oneVec := utils.NewVector(len(oneVals), oneVals)
				oneCoeffs := VInv.Mul(oneVec.ToMatrix())

				// Sum weighted basis functions
				sum := 0.0
				for i := 0; i < basis.Np; i++ {
					sum += oneCoeffs.At(i, 0) * phi[i]
				}

				if math.Abs(sum-1.0) > 1e-10 {
					t.Errorf("Order %d: partition of unity failed at (%g,%g,%g): sum=%g",
						P, rTest, sTest, tTest, sum)
				}
			}
		})
	}
}

// Helper function for random number generation
var rand = struct {
	Float64 func() float64
}{
	Float64: func() float64 {
		// Simple linear congruential generator for reproducibility
		var seed uint64 = 12345
		seed = (seed*1103515245 + 12345) & 0x7fffffff
		return float64(seed) / float64(0x7fffffff)
	},
}
