package gonudg

import (
	"math"
	"testing"
)

func TestEquiNodes3D(t *testing.T) {
	tests := []struct {
		name string
		p    int
		wantNp int
	}{
		{"p=0", 0, 1},
		{"p=1", 1, 4},
		{"p=2", 2, 10},
		{"p=3", 3, 20},
		{"p=4", 4, 35},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			r, s, tt := EquiNodes3D(tc.p)
			
			// Check number of nodes
			if len(r) != tc.wantNp || len(s) != tc.wantNp || len(tt) != tc.wantNp {
				t.Errorf("Wrong number of nodes: got %d, want %d", len(r), tc.wantNp)
			}
			
			// Check that all nodes are in the reference tetrahedron
			for i := 0; i < len(r); i++ {
				if !IsValidRST(r[i], s[i], tt[i]) {
					t.Errorf("Node %d at (%f,%f,%f) is outside reference tetrahedron",
						i, r[i], s[i], tt[i])
				}
			}
			
			// Check that nodes are equidistributed
			// For p=1, should have nodes at the 4 vertices
			if tc.p == 1 {
				vertices := [][]float64{
					{-1, -1, -1},
					{1, -1, -1},
					{-1, 1, -1},
					{-1, -1, 1},
				}
				
				for _, v := range vertices {
					found := false
					for i := 0; i < len(r); i++ {
						if math.Abs(r[i]-v[0]) < 1e-14 &&
						   math.Abs(s[i]-v[1]) < 1e-14 &&
						   math.Abs(tt[i]-v[2]) < 1e-14 {
							found = true
							break
						}
					}
					if !found {
						t.Errorf("Vertex (%f,%f,%f) not found in nodes", v[0], v[1], v[2])
					}
				}
			}
		})
	}
}

func TestNodes3D(t *testing.T) {
	tests := []struct {
		name string
		p    int
	}{
		{"p=1", 1},
		{"p=2", 2},
		{"p=3", 3},
		{"p=4", 4},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			X, Y, Z := Nodes3D(tc.p)
			
			expectedNp := (tc.p + 1) * (tc.p + 2) * (tc.p + 3) / 6
			
			// Check number of nodes
			if len(X) != expectedNp || len(Y) != expectedNp || len(Z) != expectedNp {
				t.Errorf("Wrong number of nodes: got %d, want %d", len(X), expectedNp)
			}
			
			// Transform to reference coordinates and verify
			r, s, tt := XYZtoRST(X, Y, Z)
			
			// All nodes should be in the reference tetrahedron
			for i := 0; i < len(r); i++ {
				if !IsValidRST(r[i], s[i], tt[i]) {
					t.Errorf("Node %d at rst=(%f,%f,%f) is outside reference tetrahedron",
						i, r[i], s[i], tt[i])
				}
			}
			
			// For p=1, check that we get the vertices of the equilateral tetrahedron
			if tc.p == 1 {
				sqrt3 := math.Sqrt(3.0)
				sqrt6 := math.Sqrt(6.0)
				
				vertices := [][]float64{
					{-1.0, -1.0/sqrt3, -1.0/sqrt6},
					{1.0, -1.0/sqrt3, -1.0/sqrt6},
					{0.0, 2.0/sqrt3, -1.0/sqrt6},
					{0.0, 0.0, 3.0/sqrt6},
				}
				
				for _, v := range vertices {
					found := false
					tol := 1e-10
					for i := 0; i < len(X); i++ {
						if math.Abs(X[i]-v[0]) < tol &&
						   math.Abs(Y[i]-v[1]) < tol &&
						   math.Abs(Z[i]-v[2]) < tol {
							found = true
							break
						}
					}
					if !found {
						t.Errorf("Vertex (%f,%f,%f) not found in nodes", v[0], v[1], v[2])
					}
				}
			}
		})
	}
}

func TestWarpShiftFace3D(t *testing.T) {
	// Test basic properties of warp shift
	p := 3
	n := 10
	
	// Create test barycentric coordinates
	La := make([]float64, n)
	Lb := make([]float64, n)
	Lc := make([]float64, n)
	Ld := make([]float64, n)
	
	for i := 0; i < n; i++ {
		// Create points on a face (La = 0)
		La[i] = 0.0
		Lb[i] = 0.3 + 0.4*float64(i)/float64(n-1)
		Lc[i] = 0.3 + 0.2*float64(i)/float64(n-1)
		Ld[i] = 1.0 - La[i] - Lb[i] - Lc[i]
	}
	
	alpha := 1.0
	warpx, warpy := WarpShiftFace3D(p, alpha, alpha, La, Lb, Lc, Ld)
	
	// Check that warp is defined
	if len(warpx) != n || len(warpy) != n {
		t.Errorf("Wrong warp dimensions: got (%d,%d), want (%d,%d)",
			len(warpx), len(warpy), n, n)
	}
	
	// Warp should be zero at vertices
	// When two of Lb, Lc, Ld are zero, warp should be zero
	La2 := []float64{0.0}
	Lb2 := []float64{1.0}
	Lc2 := []float64{0.0}
	Ld2 := []float64{0.0}
	
	wx2, wy2 := WarpShiftFace3D(p, alpha, alpha, La2, Lb2, Lc2, Ld2)
	
	tol := 1e-10
	if math.Abs(wx2[0]) > tol || math.Abs(wy2[0]) > tol {
		t.Errorf("Warp at vertex should be zero: got (%f,%f)", wx2[0], wy2[0])
	}
}

func TestEvalwarp(t *testing.T) {
	// Test the warp evaluation function
	p := 3
	gaussX := JacobiGL(0, 0, p)
	
	// Negate (as in the code)
	for i := range gaussX {
		gaussX[i] = -gaussX[i]
	}
	
	// Test at equidistant nodes - warp should give the difference
	xeq := make([]float64, p+1)
	for i := 0; i <= p; i++ {
		xeq[i] = -1.0 + 2.0*float64(i)/float64(p)
	}
	
	warp := evalwarp(p, gaussX, xeq)
	
	// At endpoints, warp should be zero (sf = 0)
	tol := 1e-14
	if math.Abs(warp[0]) > tol || math.Abs(warp[p]) > tol {
		t.Errorf("Warp at endpoints should be zero: got %f and %f", warp[0], warp[p])
	}
}

func BenchmarkNodes3D(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Nodes3D(5)
	}
}

