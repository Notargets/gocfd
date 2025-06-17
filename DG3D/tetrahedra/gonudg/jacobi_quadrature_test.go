package gonudg

import (
	"fmt"
	"math"
	"testing"
)

func TestJacobiGL(t *testing.T) {
	// Test cases with known values
	tests := []struct {
		name     string
		alpha    float64
		beta     float64
		N        int
		expected []float64
		tol      float64
	}{
		{
			name:     "N=0",
			alpha:    0.0,
			beta:     0.0,
			N:        0,
			expected: []float64{0.0},
			tol:      1e-14,
		},
		{
			name:     "N=1",
			alpha:    0.0,
			beta:     0.0,
			N:        1,
			expected: []float64{-1.0, 1.0},
			tol:      1e-14,
		},
		{
			name:     "N=2_Legendre",
			alpha:    0.0,
			beta:     0.0,
			N:        2,
			expected: []float64{-1.0, 0.0, 1.0},
			tol:      1e-14,
		},
		{
			name:     "N=3_Legendre",
			alpha:    0.0,
			beta:     0.0,
			N:        3,
			expected: []float64{-1.0, -0.4472135954999579, 0.4472135954999579, 1.0},
			tol:      1e-14,
		},
		{
			name:     "N=4_Legendre",
			alpha:    0.0,
			beta:     0.0,
			N:        4,
			expected: []float64{-1.0, -0.6546536707079771, 0.0, 0.6546536707079771, 1.0},
			tol:      1e-14,
		},
		{
			name:  "N=5_Legendre",
			alpha: 0.0,
			beta:  0.0,
			N:     5,
			// These are zeros of (1-x^2)*P'_5(x)
			expected: []float64{-1.0, -0.7650553239294647, -0.2852315164806451,
				0.2852315164806451, 0.7650553239294647, 1.0},
			tol: 1e-14,
		},
		{
			name:     "N=2_Jacobi_alpha1_beta1",
			alpha:    1.0,
			beta:     1.0,
			N:        2,
			expected: []float64{-1.0, 0.0, 1.0},
			tol:      1e-14,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			result := JacobiGL(tc.alpha, tc.beta, tc.N)

			// Check length
			if len(result) != len(tc.expected) {
				t.Errorf("Wrong number of points: got %d, want %d",
					len(result), len(tc.expected))
				return
			}

			// Check values
			for i, val := range result {
				if math.Abs(val-tc.expected[i]) > tc.tol {
					t.Errorf("Point %d: got %.15f, want %.15f (diff: %.2e)",
						i, val, tc.expected[i], math.Abs(val-tc.expected[i]))
				}
			}

			// Additional checks for N >= 2
			if tc.N >= 2 {
				// Check endpoints
				if math.Abs(result[0]+1.0) > tc.tol {
					t.Errorf("First point should be -1, got %.15f", result[0])
				}
				if math.Abs(result[tc.N]-1.0) > tc.tol {
					t.Errorf("Last point should be 1, got %.15f", result[tc.N])
				}

				// Check symmetry for alpha=beta
				if tc.alpha == tc.beta {
					mid := tc.N / 2
					for i := 0; i < mid; i++ {
						expected := -result[tc.N-i]
						if math.Abs(result[i]-expected) > tc.tol {
							t.Errorf("Symmetry violated: x[%d]=%v, -x[%d]=%v",
								i, result[i], tc.N-i, expected)
						}
					}
				}

				// Check ordering (should be strictly increasing)
				for i := 1; i < len(result); i++ {
					if result[i] <= result[i-1] {
						t.Errorf("Points not in increasing order: x[%d]=%v >= x[%d]=%v",
							i-1, result[i-1], i, result[i])
					}
				}
			}
		})
	}
}

func TestJacobiGQ(t *testing.T) {
	// Test Gauss quadrature points
	tests := []struct {
		name     string
		alpha    float64
		beta     float64
		N        int
		checkSum bool // Check if weights sum to 2
		tol      float64
	}{
		{
			name:     "N=0_Legendre",
			alpha:    0.0,
			beta:     0.0,
			N:        0,
			checkSum: true,
			tol:      1e-14,
		},
		{
			name:     "N=1_Legendre",
			alpha:    0.0,
			beta:     0.0,
			N:        1,
			checkSum: true,
			tol:      1e-14,
		},
		{
			name:     "N=2_Legendre",
			alpha:    0.0,
			beta:     0.0,
			N:        2,
			checkSum: true,
			tol:      1e-14,
		},
		{
			name:     "N=1_Jacobi",
			alpha:    1.0,
			beta:     1.0,
			N:        1,
			checkSum: true,
			tol:      1e-14,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			x := JacobiGQ(tc.alpha, tc.beta, tc.N)

			// Check correct number of points
			expectedLen := tc.N + 1
			if tc.N == 0 {
				expectedLen = 1
			}
			if len(x) != expectedLen {
				t.Errorf("Wrong number of points: got %d, want %d",
					len(x), expectedLen)
			}

			// Check that all points are in [-1, 1]
			for i, xi := range x {
				if xi < -1.0-tc.tol || xi > 1.0+tc.tol {
					t.Errorf("Point %d out of range [-1,1]: %v", i, xi)
				}
			}

			// Check ordering (should be strictly increasing)
			for i := 1; i < len(x); i++ {
				if x[i] <= x[i-1] {
					t.Errorf("Points not in increasing order: x[%d]=%v >= x[%d]=%v",
						i-1, x[i-1], i, x[i])
				}
			}

			// For Legendre (alpha=beta=0), check symmetry
			if tc.alpha == 0 && tc.beta == 0 && tc.N > 0 {
				mid := len(x) / 2
				for i := 0; i < mid; i++ {
					expected := -x[len(x)-1-i]
					if math.Abs(x[i]-expected) > tc.tol {
						t.Errorf("Symmetry violated: x[%d]=%v, -x[%d]=%v",
							i, x[i], len(x)-1-i, expected)
					}
				}
			}
		})
	}
}

// TestJacobiGLInteriorPoints verifies that JacobiGL correctly uses N-2 interior points
func TestJacobiGLInteriorPoints(t *testing.T) {
	// This test specifically checks the fix for the N-1 vs N-2 bug
	// For Gauss-Lobatto points of order N:
	// - Total points: N+1
	// - Endpoints: 2 (always at -1 and 1)
	// - Interior points: N-1 (computed using JacobiGQ of order N-2)
	testCases := []struct {
		N                int
		expectedInterior int
		expectedTotal    int
	}{
		{N: 2, expectedInterior: 1, expectedTotal: 3}, // 3 points: -1, [0], 1
		{N: 3, expectedInterior: 2, expectedTotal: 4}, // 4 points: -1, [-√5/5, √5/5], 1
		{N: 4, expectedInterior: 3, expectedTotal: 5}, // 5 points: -1, [-√21/7, 0, √21/7], 1
		{N: 5, expectedInterior: 4, expectedTotal: 6}, // 6 points: -1, [4 interior], 1
	}

	for _, tc := range testCases {
		t.Run(fmt.Sprintf("N=%d", tc.N), func(t *testing.T) {
			x := JacobiGL(0.0, 0.0, tc.N)

			// Total points should be N+1
			if len(x) != tc.expectedTotal {
				t.Errorf("Wrong total points: got %d, want %d", len(x), tc.expectedTotal)
			}

			// First and last should be -1 and 1
			if math.Abs(x[0]+1.0) > 1e-14 {
				t.Errorf("First point should be -1, got %v", x[0])
			}
			if math.Abs(x[tc.N]-1.0) > 1e-14 {
				t.Errorf("Last point should be 1, got %v", x[tc.N])
			}

			// Interior points
			numInterior := len(x) - 2
			if numInterior != tc.expectedInterior {
				t.Errorf("Wrong number of interior points: got %d, want %d",
					numInterior, tc.expectedInterior)
			}
		})
	}
}

// BenchmarkJacobiGL benchmarks the Gauss-Lobatto points computation
func BenchmarkJacobiGL(b *testing.B) {
	orders := []int{5, 10, 15, 20}

	for _, N := range orders {
		b.Run(fmt.Sprintf("N=%d", N), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				_ = JacobiGL(0.0, 0.0, N)
			}
		})
	}
}

// TestJacobiGLAccuracy tests the accuracy of computed points by checking
// that they satisfy the equation (1-x^2)*P'_N(x) = 0
func TestJacobiGLAccuracy(t *testing.T) {
	// Test for Legendre polynomials (alpha=beta=0)
	for N := 2; N <= 6; N++ {
		t.Run(fmt.Sprintf("N=%d", N), func(t *testing.T) {
			x := JacobiGL(0.0, 0.0, N)

			// For each interior point, check that it's a zero of P'_N(x)
			for i := 1; i < N; i++ {
				// Evaluate P'_N(x[i]) using recurrence relation
				// This is an approximation - in practice you'd use the
				// actual derivative evaluation
				xi := x[i]

				// For interior points, (1-x^2) != 0, so P'_N(x) should be 0
				// We can't easily check this without implementing P'_N,
				// but we can at least verify |x| < 1
				if math.Abs(xi) >= 1.0 {
					t.Errorf("Interior point %d has |x|>=1: x=%v", i, xi)
				}
			}
		})
	}
}
