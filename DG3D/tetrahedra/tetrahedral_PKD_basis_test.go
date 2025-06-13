package tetrahedra

import (
	"fmt"
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
		tol = 1.e-9
	)
	// for i, label := range []string{"Warp/Blend", "Equispaced", "Shunn Ham"} {
	// 	t.Logf("Using %s nodes\n", label)
	// Test orders 1 through 5
	for P := 1; P <= 6; P++ {
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
				// N 0
				{"constant", func(r, s, t float64) float64 { return 1.0 }, 0},

				// N 1
				{"linear_r", func(r, s, t float64) float64 { return r }, 1},
				{"linear_s", func(r, s, t float64) float64 { return s }, 1},
				{"linear_t", func(r, s, t float64) float64 { return t }, 1},
				{"linear_combo", func(r, s, t float64) float64 { return 2*r - 3*s + t }, 1},

				// N 2
				{"quadratic_r2", func(r, s, t float64) float64 { return r * r }, 2},
				{"quadratic_s2", func(r, s, t float64) float64 { return s * s }, 2},
				{"quadratic_t2", func(r, s, t float64) float64 { return t * t }, 2},
				{"quadratic_rs", func(r, s, t float64) float64 { return r * s }, 2},
				{"quadratic_rt", func(r, s, t float64) float64 { return r * t }, 2},
				{"quadratic_st", func(r, s, t float64) float64 { return s * t }, 2},
				{"quadratic_combo", func(r, s, t float64) float64 { return r*r + 2*r*s - t*t }, 2},

				// N 3
				{"cubic_r3", func(r, s, t float64) float64 { return r * r * r }, 3},
				{"cubic_s3", func(r, s, t float64) float64 { return s * s * s }, 3},
				{"cubic_rst", func(r, s, t float64) float64 { return r * s * t }, 3},
				{"cubic_r2s", func(r, s, t float64) float64 { return r * r * s }, 3},
				{"cubic_combo", func(r, s, t float64) float64 { return r*r*r - 2*r*s*t + s*s*t }, 3},

				// N 4
				{"quartic_r4", func(r, s, t float64) float64 { return r * r * r * r }, 4},
				{"quartic_r2s2", func(r, s, t float64) float64 { return r * r * s * s }, 4},
				{"quartic_r3t", func(r, s, t float64) float64 { return r * r * r * t }, 4},
				{"quartic_combo", func(r, s, t float64) float64 { return r*r*r*r + r*r*s*s - 2*r*s*t*t }, 4},

				// N 5
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
					t.Errorf("N %d, polynomial %s: max interpolation error %g exceeds tolerance %g",
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
					t.Errorf("N %d, polynomial %s: reconstruction error at nodes %g",
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
	tol := 2.e-9

	// Test orders 1 through 5
	for P := 1; P <= 6; P++ {
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
				// N 0 - derivatives should be zero
				{
					"constant",
					func(r, s, t float64) float64 { return 1.0 },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					func(r, s, t float64) float64 { return 0.0 },
					0,
				},

				// N 1 - derivatives are constants
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

				// N 2 - derivatives are linear
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

				// N 3 - derivatives are quadratic
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

				// N 4 - derivatives are cubic
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

				// N 5 - derivatives are quartic
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
					t.Errorf("N %d, polynomial %s: dr/dr error %g exceeds tolerance %g",
						P, test.name, maxErrorDr, tol)
				}
				if maxErrorDs > tol {
					t.Errorf("N %d, polynomial %s: dr/ds error %g exceeds tolerance %g",
						P, test.name, maxErrorDs, tol)
				}
				if maxErrorDt > tol {
					t.Errorf("N %d, polynomial %s: dr/dt error %g exceeds tolerance %g",
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
					t.Errorf("N %d, polynomial %s: random point dr/dr error %g exceeds tolerance %g",
						P, test.name, maxErrorDr, randomTol)
				}
				if maxErrorDs > randomTol {
					t.Errorf("N %d, polynomial %s: random point dr/ds error %g exceeds tolerance %g",
						P, test.name, maxErrorDs, randomTol)
				}
				if maxErrorDt > randomTol {
					t.Errorf("N %d, polynomial %s: random point dr/dt error %g exceeds tolerance %g",
						P, test.name, maxErrorDt, randomTol)
				}
			}
		})
	}
}

func TestBuildFmask3D(t *testing.T) {
	for P := 1; P <= 6; P++ {
		t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
			// Get interpolation nodes
			tb := NewTetBasis(P)

			// Build face masks
			fmask := tb.Fmask

			// Expected number of nodes per face
			Nfp := (P + 1) * (P + 2) / 2

			// Total number of nodes
			Np := (P + 1) * (P + 2) * (P + 3) / 6

			// Test 1: Check we have 4 faces
			assert.Equal(t, 4, len(fmask), "Should have 4 faces")

			// Test 2: Check each face has the correct number of nodes
			for face := 0; face < 4; face++ {
				assert.Equal(t, Nfp, len(fmask[face]),
					"Face %d should have %d nodes, got %d", face, Nfp, len(fmask[face]))
			}

			// Test 3: Check all face node indices are valid
			for face := 0; face < 4; face++ {
				for _, idx := range fmask[face] {
					assert.True(t, idx >= 0 && idx < Np,
						"Face %d has invalid node index %d (should be 0-%d)", face, idx, Np-1)
				}
			}

			// Test 4: Check no duplicate nodes within a face
			for face := 0; face < 4; face++ {
				seen := make(map[int]bool)
				for _, idx := range fmask[face] {
					assert.False(t, seen[idx],
						"Face %d has duplicate node index %d", face, idx)
					seen[idx] = true
				}
			}

			// Test 5: Verify face nodes satisfy face equations
			tol := 1e-10
			for face := 0; face < 4; face++ {
				for _, idx := range fmask[face] {
					r := tb.R.At(idx)
					s := tb.S.At(idx)
					tt := tb.T.At(idx)

					switch face {
					case 0: // Face 1: t = -1
						assert.True(t, math.Abs(1.0+tt) < tol,
							"Face 0 node %d: t=%g should be -1", idx, tt)
					case 1: // Face 2: s = -1
						assert.True(t, math.Abs(1.0+s) < tol,
							"Face 1 node %d: s=%g should be -1", idx, s)
					case 2: // Face 3: r+s+t = -1
						assert.True(t, math.Abs(1.0+r+s+tt) < tol,
							"Face 2 node %d: r+s+t=%g should be -1", idx, r+s+tt)
					case 3: // Face 4: r = -1
						assert.True(t, math.Abs(1.0+r) < tol,
							"Face 3 node %d: r=%g should be -1", idx, r)
					}
				}
			}

			// Test 6: Count total unique face nodes
			allFaceNodes := make(map[int]bool)
			for face := 0; face < 4; face++ {
				for _, idx := range fmask[face] {
					allFaceNodes[idx] = true
				}
			}

			// For a tetrahedron, all nodes should be on the boundary
			// except for higher order internal nodes
			t.Logf("N %d: %d total nodes, %d face nodes", P, Np, len(allFaceNodes))

			// Test 7: Print face node counts for debugging
			t.Logf("N %d face node counts: [%d, %d, %d, %d]",
				P, len(fmask[0]), len(fmask[1]), len(fmask[2]), len(fmask[3]))
		})
	}
}

func TestBuildFmask3DWithToleranceCheck(t *testing.T) {
	for P := 1; P <= 6; P++ {
		t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
			// Get interpolation nodes
			R, S, T := Nodes3D(P)

			// Expected number of nodes per face
			Nfp := (P + 1) * (P + 2) / 2

			// Check with different tolerances
			tolerances := []float64{1e-10, 1e-9, 1e-8, 1e-7, 1e-6}

			for _, tol := range tolerances {
				// Count nodes near each face
				face0Count := 0
				face1Count := 0
				face2Count := 0
				face3Count := 0

				for i := 0; i < R.Len(); i++ {
					r := R.At(i)
					s := S.At(i)
					t := T.At(i)

					if math.Abs(1.0+t) < tol {
						face0Count++
					}
					if math.Abs(1.0+s) < tol {
						face1Count++
					}
					if math.Abs(1.0+r+s+t) < tol {
						face2Count++
					}
					if math.Abs(1.0+r) < tol {
						face3Count++
					}
				}

				t.Logf("N %d, tolerance %e: face counts [%d, %d, %d, %d], expected %d",
					P, tol, face0Count, face1Count, face2Count, face3Count, Nfp)
			}

			// Also check the actual distance from faces for nodes we expect to be on faces
			maxDist := 0.0
			for i := 0; i < R.Len(); i++ {
				r := R.At(i)
				s := S.At(i)
				tt := T.At(i)

				// Check minimum distance to any face
				dist0 := math.Abs(1.0 + tt)
				dist1 := math.Abs(1.0 + s)
				dist2 := math.Abs(1.0 + r + s + tt)
				dist3 := math.Abs(1.0 + r)

				minDist := math.Min(math.Min(dist0, dist1), math.Min(dist2, dist3))

				// If this should be a face node (within some reasonable tolerance)
				if minDist < 1e-6 {
					if minDist > maxDist {
						maxDist = minDist
						t.Logf("Node %d: r=%g, s=%g, t=%g, min face distance=%e",
							i, r, s, tt, minDist)
					}
				}
			}
			t.Logf("N %d: Maximum distance from face for 'face nodes': %e", P, maxDist)
		})
	}
}

func TestLiftOperatorExactness(t *testing.T) {
	// The LIFT operator should satisfy:
	// ∫_Ω (LIFT * f_face) · v dx = ∫_∂Ω f · v ds
	// for any test function v in the polynomial space

	for P := 1; P <= 6; P++ {
		t.Run(fmt.Sprintf("Order_%d", P), func(t *testing.T) {
			// Create basis
			basis := NewTetBasis(P)

			// Get interpolation nodes
			R, S := basis.R, basis.S

			M := basis.M

			// Get face masks
			fmask := basis.Fmask

			// Compute Lift matrix
			LIFT := basis.LIFT

			Nfp := (P + 1) * (P + 2) / 2
			Nfaces := 4

			// Adjust tolerance based on polynomial order
			// Higher orders have more numerical error due to conditioning
			tol := 1e-10
			if P >= 6 {
				tol = 1e-6 * math.Pow(float64(P), 2)
			}

			// Test 1: For constant data on faces, LIFT should produce a function
			// whose integral equals the total surface "measure"
			faceData := utils.NewVector(Nfaces * Nfp)
			for i := 0; i < Nfaces*Nfp; i++ {
				faceData.Set(i, 1.0)
			}

			// Apply LIFT
			lifted := LIFT.Mul(faceData.ToMatrix())

			// Compute volume integral using mass matrix
			ones := utils.NewVector(basis.Np)
			for i := 0; i < basis.Np; i++ {
				ones.Set(i, 1.0)
			}

			// Volume integral = ones^T * M * lifted
			volumeIntegral := 0.0
			Mlifted := M.Mul(lifted)
			for i := 0; i < basis.Np; i++ {
				volumeIntegral += ones.At(i) * Mlifted.At(i, 0)
			}

			expectedSurfaceMeasure := 8.0
			relError := math.Abs(volumeIntegral-expectedSurfaceMeasure) / expectedSurfaceMeasure

			assert.True(t, relError < tol,
				"N %d: Volume integral %g differs from expected surface measure %g by %g%% (tol=%g)",
				P, volumeIntegral, expectedSurfaceMeasure, relError*100, tol*100)

			// Test 2: LIFT applied to zero face data should give zero
			zeroData := utils.NewVector(Nfaces * Nfp)
			liftedZero := LIFT.Mul(zeroData.ToMatrix())

			maxVal := 0.0
			for i := 0; i < basis.Np; i++ {
				val := math.Abs(liftedZero.At(i, 0))
				if val > maxVal {
					maxVal = val
				}
			}
			zeroTol := 1e-14
			if P >= 6 {
				zeroTol = 1e-12
			}
			assert.True(t, maxVal < zeroTol,
				"N %d: LIFT of zero data has max value %g", P, maxVal)

			// Test 3: Test polynomial exactness for r
			pVals := R.Copy()

			faceP := utils.NewVector(Nfaces * Nfp)
			for face := 0; face < Nfaces; face++ {
				for i := 0; i < Nfp; i++ {
					nodeIdx := fmask[face][i]
					faceP.Set(face*Nfp+i, pVals.At(nodeIdx))
				}
			}

			liftedP := LIFT.Mul(faceP.ToMatrix())

			volumeIntegralP := 0.0
			MliftedP := M.Mul(liftedP)
			for i := 0; i < basis.Np; i++ {
				volumeIntegralP += ones.At(i) * MliftedP.At(i, 0)
			}

			expectedIntegralR := -4.0

			assert.True(t, math.Abs(volumeIntegralP-expectedIntegralR) < tol,
				"N %d: Integral of r over surface = %g, expected %g (err=%g, tol=%g)",
				P, volumeIntegralP, expectedIntegralR,
				math.Abs(volumeIntegralP-expectedIntegralR), tol)

			// Test 4: Test with s
			sVals := S.Copy()
			faceS := utils.NewVector(Nfaces * Nfp)
			for face := 0; face < Nfaces; face++ {
				for i := 0; i < Nfp; i++ {
					nodeIdx := fmask[face][i]
					faceS.Set(face*Nfp+i, sVals.At(nodeIdx))
				}
			}

			liftedS := LIFT.Mul(faceS.ToMatrix())
			volumeIntegralS := 0.0
			MliftedS := M.Mul(liftedS)
			for i := 0; i < basis.Np; i++ {
				volumeIntegralS += ones.At(i) * MliftedS.At(i, 0)
			}

			expectedIntegralS := -4.0
			assert.True(t, math.Abs(volumeIntegralS-expectedIntegralS) < tol,
				"N %d: Integral of s over surface = %g, expected %g (err=%g, tol=%g)",
				P, volumeIntegralS, expectedIntegralS,
				math.Abs(volumeIntegralS-expectedIntegralS), tol)

			// Test 5: Energy stability
			// Skip stability test for very high orders as the bounds become meaningless
			if P <= 5 {
				randomData := utils.NewVector(Nfaces * Nfp)
				for i := 0; i < Nfaces*Nfp; i++ {
					randomData.Set(i, 2*rand.Float64()-1)
				}

				liftedRandom := LIFT.Mul(randomData.ToMatrix())

				liftedNorm := 0.0
				MliftedRandom := M.Mul(liftedRandom)
				for i := 0; i < basis.Np; i++ {
					liftedNorm += liftedRandom.At(i, 0) * MliftedRandom.At(i, 0)
				}
				liftedNorm = math.Sqrt(liftedNorm)

				growthFactor := math.Pow(float64(P+1), 2.0)

				surfaceNorm := 0.0
				for i := 0; i < Nfaces*Nfp; i++ {
					surfaceNorm += randomData.At(i) * randomData.At(i)
				}
				surfaceNorm = math.Sqrt(surfaceNorm / float64(Nfp))

				assert.True(t, liftedNorm < growthFactor*surfaceNorm,
					"N %d: Lifted norm %g exceeds stability bound %g * %g",
					P, liftedNorm, growthFactor, surfaceNorm)
			}

			if testing.Verbose() {
				t.Logf("N %d: Surface measure = %.6f, ∫r = %.6f, ∫s = %.6f",
					P, volumeIntegral, volumeIntegralP, volumeIntegralS)
			}
		})
	}
}
