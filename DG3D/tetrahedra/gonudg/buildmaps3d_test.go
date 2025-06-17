package gonudg

import (
	"math"
	"testing"
)

func TestBuildMaps3D(t *testing.T) {
	// Test with two connected tetrahedra sharing a face
	// This is a simple test case where we can verify the connectivity manually

	// Setup for order N=1 (linear elements)
	N := 1
	Np := (N + 1) * (N + 2) * (N + 3) / 6  // 4 nodes per tet
	Nfp := (N + 1) * (N + 2) / 2           // 3 nodes per face
	Nfaces := 4
	K := 2  // 2 elements

	// Create two tetrahedra that share a face
	// Tet 1: vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
	// Tet 2: vertices at (1,0,0), (1,1,0), (0,1,0), (1,0,1)
	// They share face with vertices (1,0,0), (0,1,0), and one more

	// Coordinates for all nodes (2 tets * 4 nodes = 8 nodes)
	x := []float64{
		// Tet 1
		0.0, 1.0, 0.0, 0.0,
		// Tet 2
		1.0, 1.0, 0.0, 1.0,
	}
	y := []float64{
		// Tet 1
		0.0, 0.0, 1.0, 0.0,
		// Tet 2
		0.0, 1.0, 1.0, 0.0,
	}
	z := []float64{
		// Tet 1
		0.0, 0.0, 0.0, 1.0,
		// Tet 2
		0.0, 0.0, 0.0, 1.0,
	}

	// Face masks - which nodes are on each face (0-based)
	// For reference tetrahedron with nodes at standard positions
	Fmask := [][]int{
		{0, 1, 2}, // Face 0: nodes where t = -1
		{0, 1, 3}, // Face 1: nodes where s = -1
		{1, 2, 3}, // Face 2: nodes where r+s+t = -1
		{0, 2, 3}, // Face 3: nodes where r = -1
	}

	// Element to element connectivity
	// Tet 0 face 2 connects to Tet 1 face 3
	EToE := [][]int{
		{0, 0, 1, 0}, // Tet 0: self-connected except face 2 -> Tet 1
		{1, 1, 1, 0}, // Tet 1: self-connected except face 3 -> Tet 0
	}
	EToF := [][]int{
		{0, 1, 3, 3}, // Tet 0: face 2 -> face 3 of Tet 1
		{0, 1, 2, 2}, // Tet 1: face 3 -> face 2 of Tet 0
	}

	// Call BuildMaps3D
	vmapM, vmapP, mapB, vmapB := BuildMaps3D(K, Np, Nfp, Nfaces, x, y, z, EToE, EToF, Fmask)

	// Test 1: Check array sizes
	expectedSize := Nfp * Nfaces * K
	if len(vmapM) != expectedSize || len(vmapP) != expectedSize {
		t.Errorf("Wrong array sizes: vmapM=%d, vmapP=%d, expected=%d",
			len(vmapM), len(vmapP), expectedSize)
	}

	// Test 2: vmapM should map face nodes to volume nodes correctly
	for k := 0; k < K; k++ {
		for f := 0; f < Nfaces; f++ {
			for i := 0; i < Nfp; i++ {
				idx := k*Nfaces*Nfp + f*Nfp + i
				volumeNode := Fmask[f][i] + k*Np
				if vmapM[idx] != volumeNode {
					t.Errorf("vmapM[%d] = %d, expected %d (elem %d, face %d, node %d)",
						idx, vmapM[idx], volumeNode, k, f, i)
				}
			}
		}
	}

	// Test 3: For connected faces, vmapP should point to matching nodes
	// We know Tet 0 face 2 connects to Tet 1 face 3
	// The nodes should match geometrically
	
	// Test 4: For boundary faces, vmapP should equal vmapM
	for k := 0; k < K; k++ {
		for f := 0; f < Nfaces; f++ {
			if EToE[k][f] == k && EToF[k][f] == f {
				// This is a boundary face
				for i := 0; i < Nfp; i++ {
					idx := k*Nfaces*Nfp + f*Nfp + i
					if vmapP[idx] != vmapM[idx] {
						t.Errorf("Boundary face should have vmapP=vmapM at idx %d", idx)
					}
				}
			}
		}
	}

	// Test 5: Check boundary node list
	if len(mapB) != len(vmapB) {
		t.Errorf("mapB and vmapB should have same length: %d vs %d", len(mapB), len(vmapB))
	}

	// Verify all boundary nodes are correctly identified
	for i, idx := range mapB {
		if vmapP[idx] != vmapM[idx] {
			t.Errorf("mapB[%d]=%d is not a boundary node", i, idx)
		}
		if vmapB[i] != vmapM[idx] {
			t.Errorf("vmapB[%d]=%d doesn't match vmapM[%d]=%d", i, vmapB[i], idx, vmapM[idx])
		}
	}
}

func TestBuildMaps3DSingleTet(t *testing.T) {
	// Test with a single tetrahedron (all boundary)
	N := 2  // quadratic element
	Np := (N + 1) * (N + 2) * (N + 3) / 6  // 10 nodes
	Nfp := (N + 1) * (N + 2) / 2           // 6 nodes per face
	Nfaces := 4
	K := 1

	// Create node coordinates for a single tet
	r, s, tt := EquiNodes3D(N)
	X, Y, Z := RSTtoXYZ(r, s, tt)
	
	// Create Fmask using the face detection logic
	Fmask := make([][]int, Nfaces)
	tol := 1e-10
	
	// Face 0: t = -1
	Fmask[0] = []int{}
	for i := 0; i < Np; i++ {
		if math.Abs(1.0+tt[i]) < tol {
			Fmask[0] = append(Fmask[0], i)
		}
	}
	
	// Face 1: s = -1
	Fmask[1] = []int{}
	for i := 0; i < Np; i++ {
		if math.Abs(1.0+s[i]) < tol {
			Fmask[1] = append(Fmask[1], i)
		}
	}
	
	// Face 2: r+s+t = -1
	Fmask[2] = []int{}
	for i := 0; i < Np; i++ {
		if math.Abs(1.0+r[i]+s[i]+tt[i]) < tol {
			Fmask[2] = append(Fmask[2], i)
		}
	}
	
	// Face 3: r = -1
	Fmask[3] = []int{}
	for i := 0; i < Np; i++ {
		if math.Abs(1.0+r[i]) < tol {
			Fmask[3] = append(Fmask[3], i)
		}
	}

	// Single element - all self-connected
	EToE := [][]int{{0, 0, 0, 0}}
	EToF := [][]int{{0, 1, 2, 3}}

	// Call BuildMaps3D
	vmapM, vmapP, mapB, _ := BuildMaps3D(K, Np, Nfp, Nfaces, X, Y, Z, EToE, EToF, Fmask)

	// For a single element, all faces are boundary
	// So vmapP should equal vmapM everywhere
	for i := 0; i < len(vmapM); i++ {
		if vmapP[i] != vmapM[i] {
			t.Errorf("Single element should have vmapP[%d]=vmapM[%d], got %d vs %d",
				i, i, vmapP[i], vmapM[i])
		}
	}

	// All face nodes should be boundary nodes
	expectedBoundaryNodes := Nfp * Nfaces * K
	if len(mapB) != expectedBoundaryNodes {
		t.Errorf("Expected %d boundary nodes, got %d", expectedBoundaryNodes, len(mapB))
	}
}

func BenchmarkBuildMaps3D(b *testing.B) {
	// Setup for a moderate size problem
	N := 3
	K := 100  // 100 elements
	Np := (N + 1) * (N + 2) * (N + 3) / 6
	Nfp := (N + 1) * (N + 2) / 2
	Nfaces := 4

	// Create dummy data
	x := make([]float64, Np*K)
	y := make([]float64, Np*K)
	z := make([]float64, Np*K)
	
	// Simple connectivity - mostly self-connected
	EToE := make([][]int, K)
	EToF := make([][]int, K)
	for i := 0; i < K; i++ {
		EToE[i] = []int{i, i, i, i}
		EToF[i] = []int{0, 1, 2, 3}
	}
	
	// Create Fmask
	Fmask := make([][]int, Nfaces)
	for f := 0; f < Nfaces; f++ {
		Fmask[f] = make([]int, Nfp)
		for i := 0; i < Nfp; i++ {
			Fmask[f][i] = i  // Simplified
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		BuildMaps3D(K, Np, Nfp, Nfaces, x, y, z, EToE, EToF, Fmask)
	}
}

