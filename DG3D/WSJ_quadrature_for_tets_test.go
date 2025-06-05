package DG3D

import (
	"math"
	"testing"
)

// TestQuadratureOrder4PointCount verifies we have exactly 20 points
func TestQuadratureOrder4PointCount(t *testing.T) {
	quad, err := NewTetrahedralQuadrature(4)
	if err != nil {
		t.Fatalf("Failed to create quadrature: %v", err)
	}

	if len(quad.Points) != 20 {
		t.Errorf("Expected 20 points, got %d", len(quad.Points))
	}
}

// TestQuadratureWeightSum verifies that weights sum to the tetrahedron volume
func TestQuadratureWeightSum(t *testing.T) {
	quad, err := NewTetrahedralQuadrature(4)
	if err != nil {
		t.Fatalf("Failed to create quadrature: %v", err)
	}

	weightSum := 0.0
	for _, pt := range quad.Points {
		weightSum += pt.W
	}

	expectedVolume := ReferenceTetrahedronVolume()
	tolerance := 1e-10

	if math.Abs(weightSum-expectedVolume) > tolerance {
		t.Errorf("Weight sum %f does not match tetrahedron volume %f", weightSum, expectedVolume)
	}
}

// TestQuadraturePointsInTetrahedron verifies all points are inside the reference tetrahedron
func TestQuadraturePointsInTetrahedron(t *testing.T) {
	quad, err := NewTetrahedralQuadrature(4)
	if err != nil {
		t.Fatalf("Failed to create quadrature: %v", err)
	}

	// For the [-1,1]³ reference tetrahedron with vertices at
	// (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
	// A point (r,s,t) is inside if:
	// r + s + t <= 1 and r,s,t >= -1
	for i, pt := range quad.Points {
		if pt.R < -1 || pt.S < -1 || pt.T < -1 {
			t.Errorf("Point %d at (%f,%f,%f) has coordinate < -1", i, pt.R, pt.S, pt.T)
		}

		if pt.R+pt.S+pt.T > 1.0+1e-10 {
			t.Errorf("Point %d at (%f,%f,%f) violates r+s+t <= 1 constraint", i, pt.R, pt.S, pt.T)
		}
	}
}

// TestQuadraturePolynomialExactness verifies exact integration of polynomials up to degree 4
func TestQuadraturePolynomialExactness(t *testing.T) {
	quad, err := NewTetrahedralQuadrature(4)
	if err != nil {
		t.Fatalf("Failed to create quadrature: %v", err)
	}

	tolerance := 1e-10

	// Test constant function (degree 0)
	t.Run("Constant", func(t *testing.T) {
		f := func(r, s, t float64) float64 { return 1.0 }
		result := quad.Integrate(f)
		expected := ReferenceTetrahedronVolume()
		if math.Abs(result-expected) > tolerance {
			t.Errorf("Constant integration failed: got %f, expected %f", result, expected)
		}
	})

	// Test linear functions (degree 1)
	t.Run("Linear_R", func(t *testing.T) {
		f := func(r, s, t float64) float64 { return r }
		result := quad.Integrate(f)
		// For our reference tetrahedron with vertices at
		// (1,-1,-1), (-1,1,-1), (-1,-1,1), (-1,-1,-1)
		// The centroid is at (-0.5, -0.5, -0.5)
		// So ∫r dV = Volume * r_centroid = 4/3 * (-0.5) = -2/3
		expected := -2.0 / 3.0
		if math.Abs(result-expected) > tolerance {
			t.Errorf("Linear r integration failed: got %f, expected %f", result, expected)
		}
	})

	// Test quadratic functions (degree 2)
	t.Run("Quadratic_R2", func(t *testing.T) {
		f := func(r, s, t float64) float64 { return r * r }
		result := quad.Integrate(f)
		// For the reference tetrahedron, ∫r² dV can be computed analytically
		// Using the formula for second moments
		expected := ReferenceTetrahedronVolume() / 5.0
		if math.Abs(result-expected) > tolerance {
			t.Errorf("Quadratic r² integration failed: got %f, expected %f", result, expected)
		}
	})

	// Test cubic function (degree 3)
	t.Run("Cubic_RST", func(t *testing.T) {
		f := func(r, s, t float64) float64 { return r * s * t }
		result := quad.Integrate(f)
		// For symmetry reasons and the specific reference tetrahedron
		// ∫rst dV = -Volume/60 for our reference element
		expected := -ReferenceTetrahedronVolume() / 60.0
		if math.Abs(result-expected) > tolerance {
			t.Errorf("Cubic rst integration failed: got %f, expected %f", result, expected)
		}
	})

	// Test quartic function (degree 4)
	t.Run("Quartic_R2S2", func(t *testing.T) {
		f := func(r, s, t float64) float64 { return r * r * s * s }
		result := quad.Integrate(f)
		// Analytical result for ∫r²s² over reference tetrahedron
		expected := ReferenceTetrahedronVolume() / 105.0
		if math.Abs(result-expected) > tolerance {
			t.Errorf("Quartic r²s² integration failed: got %f, expected %f", result, expected)
		}
	})
}

// TestMonomialIntegration provides a systematic test of all monomials up to degree 4
func TestMonomialIntegration(t *testing.T) {
	quad, err := NewTetrahedralQuadrature(4)
	if err != nil {
		t.Fatalf("Failed to create quadrature: %v", err)
	}

	// Precomputed exact integrals for monomials r^i * s^j * t^k
	// over the reference tetrahedron [-1,1]³
	// These are computed as V * (-1)^(i+j+k) / ((i+j+k+3) * C(i+j+k+2, 2))
	// where V = 4/3 is the volume and C is binomial coefficient

	type monomial struct {
		i, j, k int
		exact   float64
		name    string
	}

	V := ReferenceTetrahedronVolume()

	monomials := []monomial{
		// Degree 0
		{0, 0, 0, V, "1"},

		// Degree 1
		{1, 0, 0, -V / 3.0, "r"},
		{0, 1, 0, -V / 3.0, "s"},
		{0, 0, 1, -V / 3.0, "t"},

		// Degree 2
		{2, 0, 0, V / 5.0, "r²"},
		{1, 1, 0, V / 12.0, "rs"},
		{1, 0, 1, V / 12.0, "rt"},
		{0, 2, 0, V / 5.0, "s²"},
		{0, 1, 1, V / 12.0, "st"},
		{0, 0, 2, V / 5.0, "t²"},

		// Degree 3
		{3, 0, 0, -V / 7.0, "r³"},
		{2, 1, 0, -V / 20.0, "r²s"},
		{2, 0, 1, -V / 20.0, "r²t"},
		{1, 2, 0, -V / 20.0, "rs²"},
		{1, 1, 1, -V / 60.0, "rst"},
		{1, 0, 2, -V / 20.0, "rt²"},
		{0, 3, 0, -V / 7.0, "s³"},
		{0, 2, 1, -V / 20.0, "s²t"},
		{0, 1, 2, -V / 20.0, "st²"},
		{0, 0, 3, -V / 7.0, "t³"},

		// Degree 4
		{4, 0, 0, V / 9.0, "r⁴"},
		{3, 1, 0, V / 30.0, "r³s"},
		{3, 0, 1, V / 30.0, "r³t"},
		{2, 2, 0, V / 105.0, "r²s²"},
		{2, 1, 1, V / 120.0, "r²st"},
		{2, 0, 2, V / 105.0, "r²t²"},
		{1, 3, 0, V / 30.0, "rs³"},
		{1, 2, 1, V / 120.0, "rs²t"},
		{1, 1, 2, V / 120.0, "rst²"},
		{1, 0, 3, V / 30.0, "rt³"},
		{0, 4, 0, V / 9.0, "s⁴"},
		{0, 3, 1, V / 30.0, "s³t"},
		{0, 2, 2, V / 105.0, "s²t²"},
		{0, 1, 3, V / 30.0, "st³"},
		{0, 0, 4, V / 9.0, "t⁴"},
	}

	tolerance := 1e-9

	for _, m := range monomials {
		t.Run(m.name, func(t *testing.T) {
			f := func(r, s, t float64) float64 {
				result := 1.0
				for p := 0; p < m.i; p++ {
					result *= r
				}
				for p := 0; p < m.j; p++ {
					result *= s
				}
				for p := 0; p < m.k; p++ {
					result *= t
				}
				return result
			}

			computed := quad.Integrate(f)
			if math.Abs(computed-m.exact) > tolerance {
				t.Errorf("Monomial %s: got %f, expected %f (error: %e)",
					m.name, computed, m.exact, math.Abs(computed-m.exact))
			}
		})
	}
}

// TestBarycentricConversion verifies the coordinate transformation
func TestBarycentricConversion(t *testing.T) {
	testCases := []struct {
		name    string
		xi      [4]float64
		r, s, t float64
	}{
		{
			name: "Vertex 1",
			xi:   [4]float64{1, 0, 0, 0},
			r:    1, s: -1, t: -1,
		},
		{
			name: "Vertex 2",
			xi:   [4]float64{0, 1, 0, 0},
			r:    -1, s: 1, t: -1,
		},
		{
			name: "Vertex 3",
			xi:   [4]float64{0, 0, 1, 0},
			r:    -1, s: -1, t: 1,
		},
		{
			name: "Vertex 4",
			xi:   [4]float64{0, 0, 0, 1},
			r:    -1, s: -1, t: -1,
		},
		{
			name: "Center",
			xi:   [4]float64{0.25, 0.25, 0.25, 0.25},
			r:    -0.5, s: -0.5, t: -0.5,
		},
	}

	tolerance := 1e-10

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			r, s, tt := barycentricToRST(tc.xi[0], tc.xi[1], tc.xi[2], tc.xi[3])
			if math.Abs(r-tc.r) > tolerance || math.Abs(s-tc.s) > tolerance || math.Abs(tt-tc.t) > tolerance {
				t.Errorf("Conversion failed: got (%f,%f,%f), expected (%f,%f,%f)",
					r, s, tt, tc.r, tc.s, tc.t)
			}
		})
	}
}

// TestSymmetry verifies that the quadrature rule respects tetrahedral symmetry
func TestSymmetry(t *testing.T) {
	quad, err := NewTetrahedralQuadrature(4)
	if err != nil {
		t.Fatalf("Failed to create quadrature: %v", err)
	}

	// Count points by weight
	weightCounts := make(map[float64]int)
	for _, pt := range quad.Points {
		// Round weight to avoid floating point comparison issues
		w := math.Round(pt.W*1e10) / 1e10
		weightCounts[w]++
	}

	// We expect 3 groups: 4 points, 12 points, 4 points
	if len(weightCounts) != 2 { // Note: groups 2 and 3 have same weight
		t.Errorf("Expected 2 distinct weights, got %d", len(weightCounts))
	}

	// Verify we have the right number of points in each group
	expectedCounts := map[int]bool{4: true, 16: true} // 12+4=16 for same weight
	for _, count := range weightCounts {
		if !expectedCounts[count] {
			t.Errorf("Unexpected point count %d for a weight group", count)
		}
	}
}
