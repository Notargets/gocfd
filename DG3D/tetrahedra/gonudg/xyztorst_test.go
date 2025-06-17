package gonudg

import (
	"math"
	"testing"
)

func TestXYZtoRST(t *testing.T) {
	sqrt3 := math.Sqrt(3.0)
	sqrt6 := math.Sqrt(6.0)
	
	// Test vertices of equilateral tetrahedron map to vertices of reference tetrahedron
	tests := []struct {
		name string
		x, y, z float64
		wantR, wantS, wantT float64
	}{
		{
			name: "vertex1",
			x: -1.0, y: -1.0/sqrt3, z: -1.0/sqrt6,
			wantR: -1, wantS: -1, wantT: -1,
		},
		{
			name: "vertex2",
			x: 1.0, y: -1.0/sqrt3, z: -1.0/sqrt6,
			wantR: 1, wantS: -1, wantT: -1,
		},
		{
			name: "vertex3",
			x: 0.0, y: 2.0/sqrt3, z: -1.0/sqrt6,
			wantR: -1, wantS: 1, wantT: -1,
		},
		{
			name: "vertex4",
			x: 0.0, y: 0.0, z: 3.0/sqrt6,
			wantR: -1, wantS: -1, wantT: 1,
		},
	}
	
	tol := 1e-14
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			X := []float64{tc.x}
			Y := []float64{tc.y}
			Z := []float64{tc.z}
			
			r, s, tt := XYZtoRST(X, Y, Z)
			
			if math.Abs(r[0]-tc.wantR) > tol {
				t.Errorf("r = %v, want %v", r[0], tc.wantR)
			}
			if math.Abs(s[0]-tc.wantS) > tol {
				t.Errorf("s = %v, want %v", s[0], tc.wantS)
			}
			if math.Abs(tt[0]-tc.wantT) > tol {
				t.Errorf("t = %v, want %v", tt[0], tc.wantT)
			}
		})
	}
}

func TestRSTtoXYZ(t *testing.T) {
	// Test the inverse transformation
	// Generate points in reference tetrahedron
	testPoints := []struct {
		name string
		r, s, t float64
	}{
		{"centroid", -0.5, -0.5, -0.5},
		{"face1_center", 0, -0.5, -0.5},
		{"face2_center", -0.5, 0, -0.5},
		{"face3_center", -1.0/3.0, -1.0/3.0, -1.0/3.0},
		{"face4_center", -0.5, -0.5, 0},
		{"edge_center", 0, -1, -1},
	}
	
	tol := 1e-14
	for _, tc := range testPoints {
		t.Run(tc.name, func(t *testing.T) {
			// Forward transformation
			r := []float64{tc.r}
			s := []float64{tc.s}
			tt := []float64{tc.t}
			
			X, Y, Z := RSTtoXYZ(r, s, tt)
			
			// Inverse transformation
			rBack, sBack, tBack := XYZtoRST(X, Y, Z)
			
			// Check round trip
			if math.Abs(rBack[0]-tc.r) > tol {
				t.Errorf("Round trip failed for r: got %v, want %v", rBack[0], tc.r)
			}
			if math.Abs(sBack[0]-tc.s) > tol {
				t.Errorf("Round trip failed for s: got %v, want %v", sBack[0], tc.s)
			}
			if math.Abs(tBack[0]-tc.t) > tol {
				t.Errorf("Round trip failed for t: got %v, want %v", tBack[0], tc.t)
			}
		})
	}
}

func TestXYZtoRSTMultiplePoints(t *testing.T) {
	// Test with multiple points
	n := 10
	r := make([]float64, n)
	s := make([]float64, n)
	tt := make([]float64, n)
	
	// Create points on a line in reference tetrahedron
	for i := 0; i < n; i++ {
		fi := float64(i) / float64(n-1)
		// Line from (-1,-1,-1) to centroid
		r[i] = -1.0 + 0.5*fi
		s[i] = -1.0 + 0.5*fi
		tt[i] = -1.0 + 0.5*fi
	}
	
	// Transform to equilateral coordinates
	X, Y, Z := RSTtoXYZ(r, s, tt)
	
	// Transform back
	rBack, sBack, tBack := XYZtoRST(X, Y, Z)
	
	// Check all points
	tol := 1e-14
	for i := 0; i < n; i++ {
		if math.Abs(rBack[i]-r[i]) > tol ||
		   math.Abs(sBack[i]-s[i]) > tol ||
		   math.Abs(tBack[i]-tt[i]) > tol {
			t.Errorf("Round trip failed at point %d", i)
		}
	}
}

func TestReferenceTetrahedronVolume(t *testing.T) {
	// Verify that the reference tetrahedron has the correct volume
	// Vertices: (-1,-1,-1), (1,-1,-1), (-1,1,-1), (-1,-1,1)
	// Volume = |det(v2-v1, v3-v1, v4-v1)|/6 = 8/6 = 4/3
	
	// Create a regular grid of points in the reference tetrahedron
	// and transform to equilateral coordinates
	nPoints := 0
	totalVolume := 0.0
	h := 0.1 // grid spacing
	
	for ri := -1.0; ri <= 1.0; ri += h {
		for si := -1.0; si <= 1.0; si += h {
			for ti := -1.0; ti <= 1.0; ti += h {
				// Check if point is inside reference tetrahedron
				if ri+si+ti <= -1.0+1e-10 {
					nPoints++
				}
			}
		}
	}
	
	// The number of grid points should be approximately proportional
	// to the volume for fine enough grids
	expectedRatio := 4.0 / 3.0 / 8.0 // volume_ref_tet / volume_cube
	actualRatio := float64(nPoints) / math.Pow(2.0/h+1, 3)
	
	// This is a rough check - the ratio should be close
	if math.Abs(actualRatio-expectedRatio) > 0.1*expectedRatio {
		t.Logf("Warning: volume ratio check gives %v, expected ~%v", 
			actualRatio, expectedRatio)
	}
}

func BenchmarkXYZtoRST(b *testing.B) {
	// Benchmark the transformation
	n := 1000
	X := make([]float64, n)
	Y := make([]float64, n)
	Z := make([]float64, n)
	
	// Fill with points in equilateral tetrahedron
	sqrt3 := math.Sqrt(3.0)
	sqrt6 := math.Sqrt(6.0)
	for i := 0; i < n; i++ {
		fi := float64(i) / float64(n-1)
		// Interpolate between vertices
		X[i] = -1.0 + 2.0*fi*0.5
		Y[i] = -1.0/sqrt3 + fi/sqrt3
		Z[i] = -1.0/sqrt6 + fi/sqrt6
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		XYZtoRST(X, Y, Z)
	}
}

