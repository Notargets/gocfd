package gonudg

import (
	"math"
	"testing"
)

func TestNewDG3D(t *testing.T) {
	// Test single tetrahedron
	VX := []float64{0, 1, 0, 0}
	VY := []float64{0, 0, 1, 0}
	VZ := []float64{0, 0, 0, 1}
	EToV := [][]int{{0, 1, 2, 3}}
	
	tests := []struct {
		name string
		N    int
		wantNp int
		wantNfp int
	}{
		{"N=1", 1, 4, 3},
		{"N=2", 2, 10, 6},
		{"N=3", 3, 20, 10},
		{"N=4", 4, 35, 15},
	}
	
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			dg, err := NewDG3D(tc.N, VX, VY, VZ, EToV)
			if err != nil {
				t.Fatalf("Failed to create DG3D: %v", err)
			}
			
			if dg.Np != tc.wantNp {
				t.Errorf("Np = %d, want %d", dg.Np, tc.wantNp)
			}
			
			if dg.Nfp != tc.wantNfp {
				t.Errorf("Nfp = %d, want %d", dg.Nfp, tc.wantNfp)
			}
			
			// Check arrays are initialized
			if len(dg.r) != tc.wantNp {
				t.Errorf("len(r) = %d, want %d", len(dg.r), tc.wantNp)
			}
			
			if len(dg.V) != tc.wantNp || len(dg.V[0]) != tc.wantNp {
				t.Errorf("V dimensions wrong")
			}
			
			if len(dg.Fmask) != 4 {
				t.Errorf("Fmask should have 4 faces")
			}
		})
	}
}

func TestDG3DCoordinates(t *testing.T) {
	// Unit tetrahedron
	VX := []float64{0, 1, 0, 0}
	VY := []float64{0, 0, 1, 0}
	VZ := []float64{0, 0, 0, 1}
	EToV := [][]int{{0, 1, 2, 3}}
	
	dg, err := NewDG3D(1, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}
	
	// Check that physical coordinates are correct
	// For N=1, nodes should be at vertices
	tol := 1e-10
	
	// Check some nodes
	for i := 0; i < dg.Np; i++ {
		x := dg.x[i][0]
		y := dg.y[i][0]
		z := dg.z[i][0]
		
		// Each node should be at a vertex
		foundVertex := false
		for v := 0; v < 4; v++ {
			if math.Abs(x-VX[v]) < tol && 
			   math.Abs(y-VY[v]) < tol && 
			   math.Abs(z-VZ[v]) < tol {
				foundVertex = true
				break
			}
		}
		
		if !foundVertex {
			t.Errorf("Node %d at (%f,%f,%f) not at a vertex", i, x, y, z)
		}
	}
}

func TestBuildFmask(t *testing.T) {
	VX := []float64{0, 1, 0, 0}
	VY := []float64{0, 0, 1, 0}
	VZ := []float64{0, 0, 0, 1}
	EToV := [][]int{{0, 1, 2, 3}}
	
	dg, err := NewDG3D(2, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}
	
	// Check Fmask dimensions
	if len(dg.Fmask) != 4 {
		t.Errorf("Fmask should have 4 faces")
	}
	
	// Each face should have Nfp nodes
	for f := 0; f < 4; f++ {
		if len(dg.Fmask[f]) != dg.Nfp {
			t.Errorf("Face %d has %d nodes, expected %d", 
				f, len(dg.Fmask[f]), dg.Nfp)
		}
	}
	
	// Check that face nodes satisfy face conditions
	tol := 1e-10
	
	// Face 0: t = -1
	for _, idx := range dg.Fmask[0] {
		if math.Abs(dg.t[idx]+1.0) > tol {
			t.Errorf("Face 0 node %d: t=%f, expected -1", idx, dg.t[idx])
		}
	}
	
	// Face 1: s = -1
	for _, idx := range dg.Fmask[1] {
		if math.Abs(dg.s[idx]+1.0) > tol {
			t.Errorf("Face 1 node %d: s=%f, expected -1", idx, dg.s[idx])
		}
	}
	
	// Face 2: r+s+t = -1
	for _, idx := range dg.Fmask[2] {
		sum := dg.r[idx] + dg.s[idx] + dg.t[idx]
		if math.Abs(sum+1.0) > tol {
			t.Errorf("Face 2 node %d: r+s+t=%f, expected -1", idx, sum)
		}
	}
	
	// Face 3: r = -1
	for _, idx := range dg.Fmask[3] {
		if math.Abs(dg.r[idx]+1.0) > tol {
			t.Errorf("Face 3 node %d: r=%f, expected -1", idx, dg.r[idx])
		}
	}
}

func BenchmarkStartUp3D(b *testing.B) {
	// Create a simple mesh
	VX := []float64{0, 1, 0, 0, 2, 1}
	VY := []float64{0, 0, 1, 0, 0, 1}
	VZ := []float64{0, 0, 0, 1, 0, 0}
	EToV := [][]int{
		{0, 1, 2, 3},
		{1, 4, 5, 3},
	}
	
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		NewDG3D(3, VX, VY, VZ, EToV)
	}
}
