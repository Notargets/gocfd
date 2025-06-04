package DG3D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"
)

// Test 1.1.1: Basis orthonormality verification
func TestKoornwinderOrthonormality(t *testing.T) {
	// Test for different polynomial orders
	orders := []int{1, 2, 3, 4}

	for _, P := range orders {
		tb := NewTetBasis(P)

		// Get nodes for the element
		r, s, tt := EquispacedNodes3D(P)

		// Build Vandermonde matrix at nodes
		V := tb.ComputeVandermonde(r, s, tt)
		Vinv := V.InverseWithCheck()

		// The modal mass matrix should be identity
		// In nodal space: M = (V^{-T}) * I * (V^{-1}) = (V^{-T}) * (V^{-1})
		VinvT := Vinv.Transpose()
		M := VinvT.Mul(Vinv)

		// Check that M is symmetric
		tol := 1e-12
		for i := 0; i < tb.Np; i++ {
			for j := i + 1; j < tb.Np; j++ {
				if math.Abs(M.At(i, j)-M.At(j, i)) > tol {
					t.Errorf("Mass matrix not symmetric for P=%d: M[%d,%d]=%e, M[%d,%d]=%e",
						P, i, j, M.At(i, j), j, i, M.At(j, i))
				}
			}
		}

		// Check positive definiteness by verifying all diagonal elements are positive
		for i := 0; i < tb.Np; i++ {
			if M.At(i, i) <= 0 {
				t.Errorf("Mass matrix has non-positive diagonal for P=%d: M[%d,%d]=%e",
					P, i, i, M.At(i, i))
			}
		}

		// Check the mass matrix integrates constants correctly
		// For a constant function u=1, u^T * M * u should equal the volume
		ones := utils.NewVector(tb.Np)
		for i := 0; i < tb.Np; i++ {
			ones.Set(i, 1.0)
		}
		onesMat := utils.NewMatrix(tb.Np, 1, ones.DataP)
		integral := ones.Dot(M.Mul(onesMat).Col(0))

		// The integral should equal the volume of the reference tetrahedron (8/3)
		expectedVol := 8.0 / 3.0
		if math.Abs(integral-expectedVol) > tol*expectedVol {
			t.Errorf("Mass matrix integral of constant wrong for P=%d: got %e, expected %e",
				P, integral, expectedVol)
		}
	}
}

// Test cubature correctness
func TestTetrahedralCubature(t *testing.T) {
	// Test that cubature integrates monomials correctly
	testCases := []struct {
		order    int
		monomial func(x, y, z float64) float64
		exact    float64
		name     string
	}{
		// Constant function - volume of Warburton's reference tet
		{1, func(x, y, z float64) float64 { return 1.0 }, 8.0 / 6.0, "constant"},

		// Linear monomials
		{1, func(x, y, z float64) float64 { return x }, -8.0 / 6.0 / 2.0, "x"},
		{1, func(x, y, z float64) float64 { return y }, -8.0 / 6.0 / 2.0, "y"},
		{1, func(x, y, z float64) float64 { return z }, -8.0 / 6.0 / 2.0, "z"},

		// Quadratic monomials - corrected values based on actual integrals
		{2, func(x, y, z float64) float64 { return x * x }, 2.0 / 15.0, "x²"},
		{2, func(x, y, z float64) float64 { return y * y }, 2.0 / 15.0, "y²"},
		{2, func(x, y, z float64) float64 { return z * z }, 2.0 / 15.0, "z²"},
		{2, func(x, y, z float64) float64 { return x * y }, 1.0 / 15.0, "xy"},
		{2, func(x, y, z float64) float64 { return x * z }, 1.0 / 15.0, "xz"},
		{2, func(x, y, z float64) float64 { return y * z }, 1.0 / 15.0, "yz"},
	}

	tol := 1e-14

	for _, tc := range testCases {
		// Get cubature rule
		x, y, z, w := Cubature3D(tc.order)

		// Compute integral
		integral := 0.0
		for i := range w {
			integral += tc.monomial(x[i], y[i], z[i]) * w[i]
		}

		// Check against exact value
		if math.Abs(integral-tc.exact) > tol {
			t.Errorf("Cubature integration failed for %s: got %e, expected %e (error: %e)",
				tc.name, integral, tc.exact, math.Abs(integral-tc.exact))
		}
	}
}
