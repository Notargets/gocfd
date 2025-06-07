package DG3D

import (
	"fmt"
	"github.com/notargets/gocfd/DG2D"
	"github.com/stretchr/testify/assert"
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
	// for i, label := range []string{"Warp/Blend", "Equispaced", "Shunn Ham"} {
	// 	t.Logf("Using %s nodes\n", label)
	// Test orders 1 through 5
	for P := 1; P <= 5; P++ {
		t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
			// Create basis
			basis := NewTetBasis(P)

			// Get interpolation nodes
			var R, S, T utils.Vector
			// switch i {
			// case 0:
			R, S, T = Nodes3D(P)
			// case 1:
			// 	R, S, T = EquispacedNodes3D(P)
			// case 2:
			// 	tol = 6.e-11
			// 	R, S, T = GetNodesShunnHam(P)
			// }
			assert.False(t, utils.IsNan(R))
			assert.False(t, utils.IsNan(S))
			assert.False(t, utils.IsNan(T))

			// Compute Vandermonde matrix at the nodes
			V := basis.ComputeVandermonde(R, S, T)

			// Get the inverse for computing modal coefficients
			VInv := V.InverseWithCheck()
			// if testing.Verbose() {
			// 	V.Print("V")
			// 	VInv.Print("VInv")
			// }
			assert.False(t, utils.IsNan(V))
			assert.False(t, utils.IsNan(VInv))

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
	// }
}

// TestPKDBasisDerivatives verifies that the PKD basis can compute derivatives
// of polynomials exactly up to order P-1
func TestPKDBasisDerivatives(t *testing.T) {
	tol := 2.e-10

	// Test orders 1 through 5
	for P := 1; P <= 5; P++ {
		t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
			// Create basis
			basis := NewTetBasis(P)

			R, S, T := Nodes3D(P)
			assert.False(t, utils.IsNan(R))
			assert.False(t, utils.IsNan(S))
			assert.False(t, utils.IsNan(T))

			// Compute Vandermonde matrix and its inverse
			V := basis.ComputeVandermonde(R, S, T)
			VInv := V.InverseWithCheck()
			assert.False(t, utils.IsNan(V))
			assert.False(t, utils.IsNan(VInv))

			// Compute gradient Vandermonde matrices
			Vr, Vs, Vt := basis.ComputeGradVandermonde(R, S, T)

			// Compute differentiation matrices: Dr = Vr * V^{-1}, etc.
			Dr := Vr.Mul(VInv)
			Ds := Vs.Mul(VInv)
			Dt := Vt.Mul(VInv)

			// After computing differentiation matrices
			assert.False(t, utils.IsNan(Dr), "Dr contains NaN")
			assert.False(t, utils.IsNan(Ds), "Ds contains NaN")
			assert.False(t, utils.IsNan(Dt), "Dt contains NaN")

			// Define test polynomials with their exact derivatives
			testPolynomials := []struct {
				name     string
				f        func(r, s, t float64) float64
				dfdr     func(r, s, t float64) float64
				dfds     func(r, s, t float64) float64
				dfdt     func(r, s, t float64) float64
				maxOrder int
			}{
				// Order 0 - derivatives should be zero
				{
					"constant",
					func(r, s, t float64) float64 { return 1.0 },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					0,
				},

				// Order 1 - derivatives are constants
				{
					"linear_r",
					func(r, s, t float64) float64 { return r },
					func(r, s, t float64) float64 { return 1.0 },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					1,
				},
				{
					"linear_s",
					func(r, s, t float64) float64 { return s },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 1.0 },
					func(r, s, t float64) float64 { return 0.0 },
					1,
				},
				{
					"linear_t",
					func(r, s, t float64) float64 { return t },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 1.0 },
					1,
				},
				{
					"linear_combo",
					func(r, s, t float64) float64 { return 2*r - 3*s + 4*t },
					func(r, s, t float64) float64 { return 2.0 },
					func(r, s, t float64) float64 { return -3.0 },
					func(r, s, t float64) float64 { return 4.0 },
					1,
				},

				// Order 2 - derivatives are linear
				{
					"quadratic_r2",
					func(r, s, t float64) float64 { return r * r },
					func(r, s, t float64) float64 { return 2 * r },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					2,
				},
				{
					"quadratic_rs",
					func(r, s, t float64) float64 { return r * s },
					func(r, s, t float64) float64 { return s },
					func(r, s, t float64) float64 { return r },
					func(r, s, t float64) float64 { return 0.0 },
					2,
				},
				{
					"quadratic_rst_sum",
					func(r, s, t float64) float64 { return r*r + s*s + t*t + r*s + r*t + s*t },
					func(r, s, t float64) float64 { return 2*r + s + t },
					func(r, s, t float64) float64 { return 2*s + r + t },
					func(r, s, t float64) float64 { return 2*t + r + s },
					2,
				},

				// Order 3 - derivatives are quadratic
				{
					"cubic_r3",
					func(r, s, t float64) float64 { return r * r * r },
					func(r, s, t float64) float64 { return 3 * r * r },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					3,
				},
				{
					"cubic_rst",
					func(r, s, t float64) float64 { return r * s * t },
					func(r, s, t float64) float64 { return s * t },
					func(r, s, t float64) float64 { return r * t },
					func(r, s, t float64) float64 { return r * s },
					3,
				},
				{
					"cubic_combo",
					func(r, s, t float64) float64 { return r*r*r + 2*r*r*s - 3*r*s*t + s*s*t },
					func(r, s, t float64) float64 { return 3*r*r + 4*r*s - 3*s*t },
					func(r, s, t float64) float64 { return 2*r*r - 3*r*t + 2*s*t },
					func(r, s, t float64) float64 { return -3*r*s + s*s },
					3,
				},

				// Order 4 - derivatives are cubic
				{
					"quartic_r4",
					func(r, s, t float64) float64 { return r * r * r * r },
					func(r, s, t float64) float64 { return 4 * r * r * r },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					4,
				},
				{
					"quartic_r2s2",
					func(r, s, t float64) float64 { return r * r * s * s },
					func(r, s, t float64) float64 { return 2 * r * s * s },
					func(r, s, t float64) float64 { return 2 * r * r * s },
					func(r, s, t float64) float64 { return 0.0 },
					4,
				},

				// Order 5 - derivatives are quartic
				{
					"quintic_r5",
					func(r, s, t float64) float64 { return r * r * r * r * r },
					func(r, s, t float64) float64 { return 5 * r * r * r * r },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					5,
				},
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

				// Compute derivatives using differentiation matrices
				fMat := fVec.ToMatrix()
				dfdr_computed := Dr.Mul(fMat)
				dfds_computed := Ds.Mul(fMat)
				dfdt_computed := Dt.Mul(fMat)
				// After computing derivatives
				assert.False(t, utils.IsNan(dfdr_computed), "Computed dr contains NaN")
				assert.False(t, utils.IsNan(dfds_computed), "Computed ds contains NaN")
				assert.False(t, utils.IsNan(dfdt_computed), "Computed dt contains NaN")

				// Check derivative accuracy at the nodes
				maxErrorDr := 0.0
				maxErrorDs := 0.0
				maxErrorDt := 0.0

				for i := 0; i < R.Len(); i++ {
					// Exact derivatives at node i
					exactDr := test.dfdr(R.At(i), S.At(i), T.At(i))
					exactDs := test.dfds(R.At(i), S.At(i), T.At(i))
					exactDt := test.dfdt(R.At(i), S.At(i), T.At(i))

					// Computed derivatives at node i
					computedDr := dfdr_computed.At(i, 0)
					computedDs := dfds_computed.At(i, 0)
					computedDt := dfdt_computed.At(i, 0)

					// Errors
					errDr := math.Abs(computedDr - exactDr)
					errDs := math.Abs(computedDs - exactDs)
					errDt := math.Abs(computedDt - exactDt)

					if errDr > maxErrorDr {
						maxErrorDr = errDr
					}
					if errDs > maxErrorDs {
						maxErrorDs = errDs
					}
					if errDt > maxErrorDt {
						maxErrorDt = errDt
					}
				}

				// Check errors are within tolerance
				if maxErrorDr > tol {
					t.Errorf("Order %d, polynomial %s: dr/dr error %g exceeds tolerance %g",
						P, test.name, maxErrorDr, tol)
				}
				if maxErrorDs > tol {
					t.Errorf("Order %d, polynomial %s: dr/ds error %g exceeds tolerance %g",
						P, test.name, maxErrorDs, tol)
				}
				if maxErrorDt > tol {
					t.Errorf("Order %d, polynomial %s: dr/dt error %g exceeds tolerance %g",
						P, test.name, maxErrorDt, tol)
				}
			}

			// Additional test: verify derivatives at random points using basis gradient evaluation
			for _, test := range testPolynomials {
				if test.maxOrder > P {
					continue
				}

				// Get modal coefficients
				fVals := make([]float64, R.Len())
				for i := 0; i < R.Len(); i++ {
					fVals[i] = test.f(R.At(i), S.At(i), T.At(i))
				}
				fVec := utils.NewVector(len(fVals), fVals)
				coeffs := VInv.Mul(fVec.ToMatrix())

				// Test at random points
				nTest := 20
				maxErrorDr := 0.0
				maxErrorDs := 0.0
				maxErrorDt := 0.0

				for iTest := 0; iTest < nTest; iTest++ {
					// Generate random point in reference tetrahedron
					L1 := rand.Float64()
					L2 := rand.Float64() * (1 - L1)
					L3 := rand.Float64() * (1 - L1 - L2)
					L4 := 1 - L1 - L2 - L3

					r := -L1 - L2 - L3 + L4
					s := -L1 - L2 + L3 - L4
					tt := -L1 + L2 - L3 - L4

					// Evaluate basis gradients at test point
					dphidr, dphids, dphidt := basis.EvaluateBasisGrad(r, s, tt)
					// After computing gradients at random points
					assert.False(t, utils.IsNan(dphidr), "Basis gradient dphidr contains NaN")
					assert.False(t, utils.IsNan(dphids), "Basis gradient dphids contains NaN")
					assert.False(t, utils.IsNan(dphidt), "Basis gradient dphidt contains NaN")

					// Compute derivatives using modal coefficients
					computedDr := 0.0
					computedDs := 0.0
					computedDt := 0.0
					for i := 0; i < basis.Np; i++ {
						computedDr += coeffs.At(i, 0) * dphidr[i]
						computedDs += coeffs.At(i, 0) * dphids[i]
						computedDt += coeffs.At(i, 0) * dphidt[i]
					}

					// Exact derivatives
					exactDr := test.dfdr(r, s, tt)
					exactDs := test.dfds(r, s, tt)
					exactDt := test.dfdt(r, s, tt)

					// Update max errors
					errDr := math.Abs(computedDr - exactDr)
					errDs := math.Abs(computedDs - exactDs)
					errDt := math.Abs(computedDt - exactDt)

					if errDr > maxErrorDr {
						maxErrorDr = errDr
					}
					if errDs > maxErrorDs {
						maxErrorDs = errDs
					}
					if errDt > maxErrorDt {
						maxErrorDt = errDt
					}
				}

				// Use slightly relaxed tolerance for random point evaluation
				randomTol := tol * 10
				if maxErrorDr > randomTol {
					t.Errorf("Order %d, polynomial %s: random point dr/dr error %g exceeds tolerance %g",
						P, test.name, maxErrorDr, randomTol)
				}
				if maxErrorDs > randomTol {
					t.Errorf("Order %d, polynomial %s: random point dr/ds error %g exceeds tolerance %g",
						P, test.name, maxErrorDs, randomTol)
				}
				if maxErrorDt > randomTol {
					t.Errorf("Order %d, polynomial %s: random point dr/dt error %g exceeds tolerance %g",
						P, test.name, maxErrorDt, randomTol)
				}
			}
		})
	}
}

func TestPKDBasisFaceLift(t *testing.T) {
	// TODO: Figure out how we're going to handle point distributions,
	//  include face, edge and vertices like Hesthaven or only interior
	for P := 1; P <= 5; P++ {
		t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
			// Create basis
			basis := NewTetBasis(P)

			// Get interpolation nodes
			R, S, T := Nodes3D(P)

			// Compute Vandermonde matrix
			V := basis.ComputeVandermonde(R, S, T)

			// Get face masks
			fmask := BuildFmask3D(R, S, T, P)

			// Compute Lift matrix
			LIFT := Lift3D(P, fmask, V, R, S, T)
			assert.False(t, LIFT.IsEmpty(), "LIFT is empty")
		})
	}
}

func parked(t *testing.T) {
	P := 1
	t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
		// Create basis
		basis := NewTetBasis(P)

		// Get interpolation nodes
		R, S, T := Nodes3D(P)

		// Compute Vandermonde matrix
		V := basis.ComputeVandermonde(R, S, T)

		// Get face masks
		fmask := BuildFmask3D(R, S, T, P)

		// Compute Lift matrix
		LIFT := Lift3D(P, fmask, V, R, S, T)

		assert.False(t, LIFT.IsEmpty(), "LIFT is empty") // Get WSJ quadrature for faces (2D triangle quadrature)

		faceQuadNodes, faceWeights := DG2D.WilliamsShunnJamesonWithWeights(P)
		Nfq := len(faceWeights) // Number of face quadrature points
		RFace, SFace := DG2D.MakeRSFromPoints(faceQuadNodes)

		// Create 2D basis for face interpolation
		jb2d := DG2D.NewJacobiBasis2D(P, RFace, SFace)

		// Test the lift operator properties
		Nfp := (P + 1) * (P + 2) / 2
		Nfaces := 4

		// Test 1: Verify LIFT matrix dimensions
		liftRows, liftCols := LIFT.Dims()
		expectedCols := Nfaces * Nfp
		if liftRows != basis.Np || liftCols != expectedCols {
			t.Errorf("Order %d: LIFT matrix has wrong dimensions (%d,%d), expected (%d,%d)",
				P, liftRows, liftCols, basis.Np, expectedCols)
		}

		// Test 2: Apply LIFT to constant function on all faces
		// This tests that surface integrals are correctly lifted to volume
		ones := utils.NewVector(Nfaces * Nfp)
		for i := 0; i < Nfaces*Nfp; i++ {
			ones.Set(i, 1.0)
		}

		liftOnes := LIFT.Mul(ones.ToMatrix())

		// Check that result is finite
		for i := 0; i < basis.Np; i++ {
			val := liftOnes.At(i, 0)
			if math.IsNaN(val) || math.IsInf(val, 0) {
				t.Errorf("Order %d: LIFT of constant produced invalid value at node %d: %g",
					P, i, val)
			}
		}

		// Test 3: Test with polynomial that's zero on some faces
		// This verifies face locality
		testVec := utils.NewVector(Nfaces * Nfp)
		// Set non-zero values only on face 0
		for i := 0; i < Nfp; i++ {
			testVec.Set(i, 1.0)
		}

		liftTest := LIFT.Mul(testVec.ToMatrix())

		// Nodes not connected to face 0 should have small contributions
		for i := 0; i < basis.Np; i++ {
			isOnFace0 := false
			for _, idx := range fmask[0] {
				if idx == i {
					isOnFace0 = true
					break
				}
			}

			val := liftTest.At(i, 0)
			if !isOnFace0 && math.Abs(val) > 1e-10 {
				// Interior nodes should have non-zero values due to lifting
				// but they should be reasonable
				if math.Abs(val) > 1e3 {
					t.Errorf("Order %d: Interior node %d has unreasonably large lift value: %g",
						P, i, val)
				}
			}
		}

		// Test 4: Verify using face quadrature
		// The WSJ quadrature can exactly integrate polynomials of order 2P
		// We'll verify that lifted face data integrates correctly

		// Create a test polynomial on each face
		for face := 0; face < Nfaces; face++ {
			// Get face coordinates
			var faceR, faceS utils.Vector
			switch face {
			case 0: // t = -1
				faceR = R.SubsetIndex(fmask[face])
				faceS = S.SubsetIndex(fmask[face])
			case 1: // s = -1
				faceR = R.SubsetIndex(fmask[face])
				faceS = T.SubsetIndex(fmask[face])
			case 2: // r+s+t = -1
				faceR = S.SubsetIndex(fmask[face])
				faceS = T.SubsetIndex(fmask[face])
			case 3: // r = -1
				faceR = S.SubsetIndex(fmask[face])
				faceS = T.SubsetIndex(fmask[face])
			}

			// Create face Vandermonde for interpolation to quadrature points
			faceV := jb2d.Vandermonde2D(P, faceR, faceS)
			faceVq := jb2d.Vandermonde2D(P, RFace, SFace)

			// Interpolation matrix from face nodes to quadrature points
			faceInterp := faceVq.Mul(faceV.InverseWithCheck())

			// Test polynomial f(r,s) = r + s
			faceData := utils.NewVector(Nfp)
			for i := 0; i < Nfp; i++ {
				faceData.Set(i, faceR.At(i)+faceS.At(i))
			}

			// Interpolate to quadrature points
			faceDataQ := faceInterp.Mul(faceData.ToMatrix())

			// Compute integral using quadrature
			integral := 0.0
			for q := 0; q < Nfq; q++ {
				integral += faceDataQ.At(q, 0) * faceWeights[q]
			}

			// For reference triangle, integral of (r+s) should be 0
			// since the triangle is symmetric about r=0, s=0
			if math.Abs(integral) > 1e-10 {
				t.Logf("Order %d, face %d: Integral of (r+s) = %g (should be near 0)",
					P, face, integral)
			}
		}

		// Test 5: Verify LIFT preserves polynomial exactness
		// For a polynomial g on the faces, LIFT(g) should be exact
		// up to polynomial order P

		// Create linear function on all faces
		gFace := utils.NewVector(Nfaces * Nfp)
		for face := 0; face < Nfaces; face++ {
			for i := 0; i < Nfp; i++ {
				nodeIdx := fmask[face][i]
				// Use a linear function: f = r + 2s + 3t
				gFace.Set(face*Nfp+i, R.At(nodeIdx)+2*S.At(nodeIdx)+3*T.At(nodeIdx))
			}
		}

		// Apply LIFT
		liftG := LIFT.Mul(gFace.ToMatrix())

		// Verify the lifted function is reasonable
		maxVal := 0.0
		for i := 0; i < basis.Np; i++ {
			val := math.Abs(liftG.At(i, 0))
			if val > maxVal {
				maxVal = val
			}
		}

		if maxVal > 1e3 {
			t.Errorf("Order %d: LIFT of linear function has unreasonably large values: max = %g",
				P, maxVal)
		}

		t.Logf("Order %d: LIFT operator test completed. Max lifted value: %g", P, maxVal)
	})
}
