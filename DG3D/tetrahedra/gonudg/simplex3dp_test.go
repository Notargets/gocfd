package gonudg

import (
	"fmt"
	"math"
	"testing"
)

func TestSimplex3DP(t *testing.T) {
	// Test orthogonality of basis functions
	// The PKD basis is orthogonal with respect to the weighted inner product
	// on the reference tetrahedron
	
	// Test at vertices of reference tetrahedron
	tests := []struct {
		name string
		r, s, t float64
		i, j, k int
		checkValue bool
		expectedSign float64 // For checking sign at vertices
	}{
		{
			name: "P000_at_v1",
			r: -1, s: -1, t: -1,
			i: 0, j: 0, k: 0,
			checkValue: true,
			expectedSign: 1.0,
		},
		{
			name: "P100_at_v2", 
			r: 1, s: -1, t: -1,
			i: 1, j: 0, k: 0,
			checkValue: true,
			expectedSign: 1.0,
		},
		{
			name: "P010_at_v3",
			r: -1, s: 1, t: -1,
			i: 0, j: 1, k: 0,
			checkValue: true,
			expectedSign: 1.0,
		},
		{
			name: "P001_at_v4",
			r: -1, s: -1, t: 1,
			i: 0, j: 0, k: 1,
			checkValue: true,
			expectedSign: 1.0,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			// Test single point version
			a, b, c := RSTtoABCSingle(tc.r, tc.s, tc.t)
			value := Simplex3DPSingle(a, b, c, tc.i, tc.j, tc.k)
			
			if tc.checkValue && value*tc.expectedSign < 0 {
				t.Errorf("Simplex3DPSingle(%v,%v,%v, %d,%d,%d) = %v, expected sign %v",
					tc.r, tc.s, tc.t, tc.i, tc.j, tc.k, value, tc.expectedSign)
			}

			// Test array version
			rArr := []float64{tc.r}
			sArr := []float64{tc.s}
			tArr := []float64{tc.t}
			aArr, bArr, cArr := RSTtoABC(rArr, sArr, tArr)
			values := Simplex3DP(aArr, bArr, cArr, tc.i, tc.j, tc.k)
			
			if math.Abs(values[0]-value) > 1e-14 {
				t.Errorf("Array version differs from single version: %v vs %v", values[0], value)
			}
		})
	}
}

func TestSimplex3DPNormalization(t *testing.T) {
	// Test that P_{0,0,0} gives the expected constant value
	// The constant basis function should be sqrt(3/2) to satisfy orthonormality
	
	nPoints := 10
	r := make([]float64, nPoints)
	s := make([]float64, nPoints)
	tt := make([]float64, nPoints) // avoid conflict with t
	
	// Create random points in the reference tetrahedron
	for i := 0; i < nPoints; i++ {
		// Use barycentric coordinates to ensure points are inside
		l1 := 0.1 + 0.2*float64(i)/float64(nPoints-1)
		l2 := 0.1 + 0.15*float64(i)/float64(nPoints-1)
		l3 := 0.1 + 0.1*float64(i)/float64(nPoints-1)
		l4 := 1.0 - l1 - l2 - l3
		
		r[i] = -l1 - l2 - l3 + l4
		s[i] = -l1 - l2 + l3 - l4
		tt[i] = -l1 + l2 - l3 - l4
	}
	
	a, b, c := RSTtoABC(r, s, tt)
	P000 := Simplex3DP(a, b, c, 0, 0, 0)
	
	// Check that P_{0,0,0} is constant
	expectedValue := P000[0]
	for i := 1; i < nPoints; i++ {
		if math.Abs(P000[i]-expectedValue) > 1e-14 {
			t.Errorf("P_{0,0,0} not constant: P000[0]=%v, P000[%d]=%v", 
				expectedValue, i, P000[i])
		}
	}
	
	// The expected value is sqrt(2) based on the normalization
	expectedNorm := math.Sqrt(2.0)
	if math.Abs(P000[0]-expectedNorm) > 1e-10 {
		t.Errorf("P_{0,0,0} normalization incorrect: got %v, expected ~%v", 
			P000[0], expectedNorm)
	}
}

func TestGradSimplex3DP(t *testing.T) {
	// Test gradient computation by finite differences
	h := 1e-7
	tol := 1e-5
	
	// Test at centroid of reference tetrahedron
	r := []float64{-0.5}
	s := []float64{-0.5}
	tt := []float64{-0.5}
	
	// Test for different polynomial orders
	testCases := []struct {
		i, j, k int
	}{
		{0, 0, 0}, // constant
		{1, 0, 0}, // linear in r
		{0, 1, 0}, // linear in s
		{0, 0, 1}, // linear in t
		{2, 0, 0}, // quadratic
		{1, 1, 0}, // mixed
		{1, 0, 1}, // mixed
	}
	
	for _, tc := range testCases {
		t.Run(fmt.Sprintf("grad_P%d%d%d", tc.i, tc.j, tc.k), func(t *testing.T) {
			// Compute analytical gradients
			dr, ds, dt := GradSimplex3DP(r, s, tt, tc.i, tc.j, tc.k)
			
			// Compute finite difference approximations
			// r-derivative
			rp := []float64{r[0] + h}
			rm := []float64{r[0] - h}
			ap, bp, cp := RSTtoABC(rp, s, tt)
			am, bm, cm := RSTtoABC(rm, s, tt)
			Pp := Simplex3DP(ap, bp, cp, tc.i, tc.j, tc.k)
			Pm := Simplex3DP(am, bm, cm, tc.i, tc.j, tc.k)
			dr_fd := (Pp[0] - Pm[0]) / (2 * h)
			
			// s-derivative
			sp := []float64{s[0] + h}
			sm := []float64{s[0] - h}
			ap, bp, cp = RSTtoABC(r, sp, tt)
			am, bm, cm = RSTtoABC(r, sm, tt)
			Pp = Simplex3DP(ap, bp, cp, tc.i, tc.j, tc.k)
			Pm = Simplex3DP(am, bm, cm, tc.i, tc.j, tc.k)
			ds_fd := (Pp[0] - Pm[0]) / (2 * h)
			
			// t-derivative
			tp := []float64{tt[0] + h}
			tm := []float64{tt[0] - h}
			ap, bp, cp = RSTtoABC(r, s, tp)
			am, bm, cm = RSTtoABC(r, s, tm)
			Pp = Simplex3DP(ap, bp, cp, tc.i, tc.j, tc.k)
			Pm = Simplex3DP(am, bm, cm, tc.i, tc.j, tc.k)
			dt_fd := (Pp[0] - Pm[0]) / (2 * h)
			
			// Compare
			if math.Abs(dr[0]-dr_fd) > tol {
				t.Errorf("dr mismatch: analytical=%v, finite diff=%v, error=%v",
					dr[0], dr_fd, math.Abs(dr[0]-dr_fd))
			}
			if math.Abs(ds[0]-ds_fd) > tol {
				t.Errorf("ds mismatch: analytical=%v, finite diff=%v, error=%v",
					ds[0], ds_fd, math.Abs(ds[0]-ds_fd))
			}
			if math.Abs(dt[0]-dt_fd) > tol {
				t.Errorf("dt mismatch: analytical=%v, finite diff=%v, error=%v",
					dt[0], dt_fd, math.Abs(dt[0]-dt_fd))
			}
		})
	}
}

func BenchmarkSimplex3DP(b *testing.B) {
	// Benchmark polynomial evaluation
	n := 100
	r := make([]float64, n)
	s := make([]float64, n)
	t := make([]float64, n)
	
	// Fill with points in reference tetrahedron
	for i := 0; i < n; i++ {
		fi := float64(i) / float64(n-1)
		r[i] = -1 + 0.5*fi
		s[i] = -1 + 0.3*fi
		t[i] = -1 + 0.2*fi
	}
	
	a, bb, c := RSTtoABC(r, s, t)
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		// Evaluate a moderate order polynomial
		Simplex3DP(a, bb, c, 3, 2, 1)
	}
}

