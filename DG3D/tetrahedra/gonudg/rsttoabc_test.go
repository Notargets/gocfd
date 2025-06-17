package gonudg

import (
	"fmt"
	"math"
	"testing"
)

func TestRSTtoABC(t *testing.T) {
	tests := []struct {
		name                string
		r, s, t             float64
		wantA, wantB, wantC float64
		tol                 float64
	}{
		// REMOVED "origin" test case - (0,0,0) is outside the reference tetrahedron
		{
			name: "vertex1",
			r:    -1, s: -1, t: -1,
			wantA: -1, wantB: -1, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "vertex2",
			r:    1, s: -1, t: -1,
			wantA: 1, wantB: -1, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "vertex3",
			r:    -1, s: 1, t: -1,
			wantA: -1, wantB: 1, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "vertex4",
			r:    -1, s: -1, t: 1,
			wantA: -1, wantB: -1, wantC: 1,
			tol: 1e-14,
		},
		// REMOVED invalid test cases that were outside the reference tetrahedron
		// The following are valid interior points:
		{
			name: "face_center_bottom",
			r:    -0.5, s: -0.5, t: -1,
			wantA: -1.0 / 3.0, wantB: -0.5, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "interior_point_1",
			r:    0, s: -0.5, t: -0.5,
			wantA: 1.0, wantB: -1.0 / 3.0, wantC: -0.5,
			tol: 1e-14,
		},
		{
			name: "interior_point_2",
			r:    -0.6, s: -0.2, t: -0.2,
			wantA: 1.0, wantB: 1.0 / 3.0, wantC: -0.2,
			tol: 1e-14,
		},
		// Test actual degenerate cases with VALID points
		{
			name: "degenerate_s_plus_t_zero_valid",
			r:    0.5, s: -0.75, t: -0.75, // r+s+t = -1, s+t = -1.5
			wantA: 1.0, wantB: -5.0 / 7.0, wantC: -0.75, // b = -5/7 ≈ -0.714
			tol: 1e-14,
		},
		{
			name: "edge_midpoint_v1_v2",
			r:    0, s: -1, t: -1, // Midpoint of edge between vertices 1 and 2
			wantA: 0, wantB: -1, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "edge_midpoint_v1_v3",
			r:    -1, s: 0, t: -1, // Midpoint of edge between vertices 1 and 3
			wantA: -1, wantB: 0, wantC: -1,
			tol: 1e-14,
		},
		{
			name: "edge_midpoint_v1_v4",
			r:    -1, s: -1, t: 0, // Midpoint of edge between vertices 1 and 4
			wantA: -1, wantB: -1, wantC: 0,
			tol: 1e-14,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			// First verify the point is valid
			if !isValidRST(tc.r, tc.s, tc.t) {
				t.Errorf("Test point (r,s,t)=(%v,%v,%v) is outside reference tetrahedron", tc.r, tc.s, tc.t)
				return
			}

			// Test single point version
			a, b, c := RSTtoABCSingle(tc.r, tc.s, tc.t)
			if math.Abs(a-tc.wantA) > tc.tol {
				t.Errorf("RSTtoABCSingle: a = %v, want %v", a, tc.wantA)
			}
			if math.Abs(b-tc.wantB) > tc.tol {
				t.Errorf("RSTtoABCSingle: b = %v, want %v", b, tc.wantB)
			}
			if math.Abs(c-tc.wantC) > tc.tol {
				t.Errorf("RSTtoABCSingle: c = %v, want %v", c, tc.wantC)
			}

			// Test array version
			rArr := []float64{tc.r}
			sArr := []float64{tc.s}
			tArr := []float64{tc.t}
			aArr, bArr, cArr := RSTtoABC(rArr, sArr, tArr)
			if math.Abs(aArr[0]-tc.wantA) > tc.tol {
				t.Errorf("RSTtoABC: a[0] = %v, want %v", aArr[0], tc.wantA)
			}
			if math.Abs(bArr[0]-tc.wantB) > tc.tol {
				t.Errorf("RSTtoABC: b[0] = %v, want %v", bArr[0], tc.wantB)
			}
			if math.Abs(cArr[0]-tc.wantC) > tc.tol {
				t.Errorf("RSTtoABC: c[0] = %v, want %v", cArr[0], tc.wantC)
			}
		})
	}
}

// Helper function to check if (r,s,t) is inside the reference tetrahedron
// The reference tetrahedron has vertices at:
// v1: (-1, -1, -1)
// v2: ( 1, -1, -1)
// v3: (-1,  1, -1)
// v4: (-1, -1,  1)
// Using barycentric coordinates, a point is inside if:
// λ₁ = -(1+r+s+t)/2, λ₂ = (1+r)/2, λ₃ = (1+s)/2, λ₄ = (1+t)/2
// with all λᵢ ≥ 0 and Σλᵢ = 1
func isValidRST(r, s, t float64) bool {
	const tol = 1e-10

	// Check bounds
	if r < -1-tol || s < -1-tol || t < -1-tol {
		return false
	}

	// Check the constraint from λ₁ ≥ 0
	// λ₁ = -(1+r+s+t)/2 ≥ 0 means 1+r+s+t ≤ 0, or r+s+t ≤ -1
	if r+s+t > -1+tol {
		return false
	}

	return true
}

func TestRSTtoABCFix(t *testing.T) {
	// Test the specific case that was failing
	testCases := []struct {
		name                            string
		r, s, t                         float64
		expectedA, expectedB, expectedC float64
	}{
		{
			name: "Node 6",
			r:    -1, s: -1, t: 0,
			expectedA: -1, expectedB: -1, expectedC: 0,
		},
		{
			name: "Node 8 (was failing)",
			r:    -1, s: 0, t: 0,
			expectedA: -1, expectedB: 1, expectedC: 0,
		},
		{
			name: "Vertex 1",
			r:    -1, s: -1, t: -1,
			expectedA: -1, expectedB: -1, expectedC: -1,
		},
		{
			name: "Vertex 2",
			r:    1, s: -1, t: -1,
			expectedA: 1, expectedB: -1, expectedC: -1,
		},
		{
			name: "Vertex 3",
			r:    -1, s: 1, t: -1,
			expectedA: -1, expectedB: 1, expectedC: -1,
		},
		{
			name: "Vertex 4",
			r:    -1, s: -1, t: 1,
			expectedA: -1, expectedB: -1, expectedC: 1,
		},
	}

	for _, tc := range testCases {
		a, b, c := RSTtoABCSingle(tc.r, tc.s, tc.t)

		passed := math.Abs(a-tc.expectedA) < 1e-10 &&
			math.Abs(b-tc.expectedB) < 1e-10 &&
			math.Abs(c-tc.expectedC) < 1e-10

		status := "PASS"
		if !passed {
			status = "FAIL"
		}

		fmt.Printf("%s: (r,s,t)=(%.1f,%.1f,%.1f) -> (a,b,c)=(%.1f,%.1f,%.1f) [%s]\n",
			tc.name, tc.r, tc.s, tc.t, a, b, c, status)

		if !passed {
			fmt.Printf("  Expected: (a,b,c)=(%.1f,%.1f,%.1f)\n",
				tc.expectedA, tc.expectedB, tc.expectedC)
			t.Errorf("%s: transformation incorrect", tc.name)
		}
	}
}

// TestRSTtoABCComprehensive verifies all cases work correctly
func TestRSTtoABCComprehensive(t *testing.T) {
	testCases := []struct {
		name    string
		r, s, t float64
		a, b, c float64
	}{
		// The problematic nodes from N=2
		{
			name: "Node 6",
			r:    -1, s: -1, t: 0,
			a: -1, b: -1, c: 0,
		},
		{
			name: "Node 8 (was broken)",
			r:    -1, s: 0, t: 0,
			a: -1, b: 1, c: 0,
		},
		// Vertices
		{
			name: "Vertex 1",
			r:    -1, s: -1, t: -1,
			a: -1, b: -1, c: -1,
		},
		{
			name: "Vertex 2",
			r:    1, s: -1, t: -1,
			a: 1, b: -1, c: -1,
		},
		{
			name: "Vertex 3",
			r:    -1, s: 1, t: -1,
			a: -1, b: 1, c: -1,
		},
		{
			name: "Vertex 4",
			r:    -1, s: -1, t: 1,
			a: -1, b: -1, c: 1,
		},
		// Edge and face centers
		{
			name: "Edge center",
			r:    0, s: -1, t: -1,
			a: 0, b: -1, c: -1,
		},
		{
			name: "Face center",
			r:    -1.0 / 3.0, s: -1.0 / 3.0, t: -1,
			a: 0, b: -1.0 / 3.0, c: -1, // FIXED expected values
		},
	}

	fmt.Println("Comprehensive RSTtoABC transformation test:")
	fmt.Println("==========================================")

	allPassed := true
	for _, tc := range testCases {
		// Test single point version
		a, b, c := RSTtoABCSingle(tc.r, tc.s, tc.t)

		errA := math.Abs(a - tc.a)
		errB := math.Abs(b - tc.b)
		errC := math.Abs(c - tc.c)

		passed := errA < 1e-10 && errB < 1e-10 && errC < 1e-10

		status := "PASS"
		if !passed {
			status = "FAIL"
			allPassed = false
		}

		fmt.Printf("%-20s: (%.3f,%.3f,%.3f) -> (%.3f,%.3f,%.3f) [%s]\n",
			tc.name, tc.r, tc.s, tc.t, a, b, c, status)

		if !passed {
			fmt.Printf("  Expected: (%.3f,%.3f,%.3f)\n", tc.a, tc.b, tc.c)
			fmt.Printf("  Errors: a=%.2e, b=%.2e, c=%.2e\n", errA, errB, errC)
			t.Errorf("%s: transformation incorrect", tc.name)
		}

		// Also test array version
		rArr := []float64{tc.r}
		sArr := []float64{tc.s}
		tArr := []float64{tc.t}
		aArr, bArr, cArr := RSTtoABC(rArr, sArr, tArr)

		if math.Abs(aArr[0]-tc.a) > 1e-10 ||
			math.Abs(bArr[0]-tc.b) > 1e-10 ||
			math.Abs(cArr[0]-tc.c) > 1e-10 {
			t.Errorf("%s: array version differs from single version", tc.name)
		}
	}

	if allPassed {
		fmt.Println("\nAll transformation tests passed!")
	}
}

// TestVandermonde3DWithAllFixes tests the complete solution
func TestVandermonde3DWithAllFixes(t *testing.T) {
	fmt.Println("\nTesting Vandermonde3D with all fixes applied:")
	fmt.Println("=============================================")

	for N := 1; N <= 3; N++ {
		fmt.Printf("\nTesting N=%d:\n", N)

		// Generate nodes
		X, Y, Z := Nodes3D(N)
		r, s, tt := XYZtoRST(X, Y, Z)

		// Build Vandermonde matrix
		V := Vandermonde3D(N, r, s, tt)

		// Try to invert
		nr, nc := V.Dims()
		fmt.Printf("  Matrix size: %dx%d\n", nr, nc)

		// Check for duplicates in (a,b,c) space
		a, b, c := RSTtoABC(r, s, tt)
		duplicates := make(map[string][]int)

		for i := 0; i < len(a); i++ {
			key := fmt.Sprintf("%.10f,%.10f,%.10f", a[i], b[i], c[i])
			duplicates[key] = append(duplicates[key], i)
		}

		hasDuplicates := false
		for key, nodes := range duplicates {
			if len(nodes) > 1 {
				fmt.Printf("  DUPLICATE nodes in (a,b,c): %v -> %s\n", nodes, key)
				hasDuplicates = true
			}
		}

		if !hasDuplicates {
			fmt.Printf("  No duplicate nodes in (a,b,c) space\n")

			// Try to invert
			defer func() {
				if r := recover(); r != nil {
					fmt.Printf("  Matrix inversion failed: %v\n", r)
					t.Errorf("N=%d: Vandermonde matrix is still singular", N)
				}
			}()

			Vinv := V.InverseWithCheck()
			fmt.Printf("  Matrix inverted successfully!\n")

			// Verify it's a good inverse
			I := V.Mul(Vinv)
			sumCols := I.SumCols()
			var sum float64
			for i := 0; i < sumCols.Len(); i++ {
				sum += sumCols.AtVec(i)
			}
			fmt.Printf("  Sum check: %.6f (expected %d)\n", sum, nr)
		} else {
			t.Errorf("N=%d: Still have duplicate nodes", N)
		}
	}
}
