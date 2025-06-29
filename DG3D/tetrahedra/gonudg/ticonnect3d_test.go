package gonudg

import (
	"testing"

	"github.com/notargets/gocfd/utils"
)

// Helper function to convert CompleteMesh to EToV format for tiConnect3D
func meshToEToV(mesh utils.CompleteMesh) [][]int {
	// Only process tet elements
	var EToV [][]int

	for _, elemSet := range mesh.Elements {
		if elemSet.Type == utils.Tet {
			for _, elem := range elemSet.Elements {
				// Convert logical node names to indices
				tetConn := make([]int, 4)
				for i, nodeName := range elem {
					tetConn[i] = mesh.Nodes.NodeMap[nodeName]
				}
				EToV = append(EToV, tetConn)
			}
		}
	}

	return EToV
}

// TestTiConnect3D_SingleTet tests the simplest case: a single tetrahedron
// Following Unit Testing Principle: Start with fundamentals
func TestTiConnect3D_SingleTet(t *testing.T) {
	tm := utils.GetStandardTestMeshes()

	// Create a single tet mesh
	singleTetMesh := utils.CompleteMesh{
		Nodes:    tm.TetraNodes,
		Elements: []utils.ElementSet{tm.SingleTet},
	}

	EToV := meshToEToV(singleTetMesh)

	EToE, EToF := tiConnect3D(EToV)

	// Verify dimensions
	if len(EToE) != 1 || len(EToF) != 1 {
		t.Fatalf("Expected 1 element, got %d", len(EToE))
	}

	if len(EToE[0]) != 4 || len(EToF[0]) != 4 {
		t.Fatalf("Expected 4 faces per element")
	}

	// Mathematical property: Single element should be self-connected on all faces
	for f := 0; f < 4; f++ {
		if EToE[0][f] != 0 {
			t.Errorf("Face %d: Expected self-connection to element 0, got %d", f, EToE[0][f])
		}
		if EToF[0][f] != f {
			t.Errorf("Face %d: Expected self-connection to face %d, got %d", f, f, EToF[0][f])
		}
	}
}

// TestTiConnect3D_TwoTets tests two tetrahedra sharing exactly one face
// Following Unit Testing Principle: Build systematically
func TestTiConnect3D_TwoTets(t *testing.T) {
	tm := utils.GetStandardTestMeshes()

	// Use the standard two-tet mesh
	EToV := meshToEToV(tm.TwoTetMesh)

	if len(EToV) != 2 {
		t.Fatalf("Expected 2 elements in TwoTetMesh, got %d", len(EToV))
	}

	EToE, EToF := tiConnect3D(EToV)

	// Count connections
	sharedFaces := 0
	boundaryFaces := 0

	for elem := 0; elem < 2; elem++ {
		for face := 0; face < 4; face++ {
			if EToE[elem][face] == elem {
				boundaryFaces++
			} else {
				sharedFaces++
			}
		}
	}

	// Mathematical property: Should have exactly 2 shared faces (one from each element)
	if sharedFaces != 2 {
		t.Errorf("Expected 2 shared face connections, got %d", sharedFaces)
	}

	// Mathematical property: Should have 6 boundary faces (4*2 - 2 shared)
	if boundaryFaces != 6 {
		t.Errorf("Expected 6 boundary faces, got %d", boundaryFaces)
	}

	// Test reciprocity: If elem A face f connects to elem B face g,
	// then elem B face g must connect to elem A face f
	testReciprocity(t, EToE, EToF)
}

// TestTiConnect3D_CubeMesh tests a more complex configuration
// Following Unit Testing Principle: Progressive complexity
func TestTiConnect3D_CubeMesh(t *testing.T) {
	tm := utils.GetStandardTestMeshes()

	// Use the standard cube mesh (6 tetrahedra)
	EToV := meshToEToV(tm.CubeMesh)

	if len(EToV) != 6 {
		t.Fatalf("Expected 6 elements in CubeMesh, got %d", len(EToV))
	}

	EToE, EToF := tiConnect3D(EToV)

	// Verify all elements have correct number of faces
	for elem := 0; elem < 6; elem++ {
		if len(EToE[elem]) != 4 || len(EToF[elem]) != 4 {
			t.Errorf("Element %d has wrong number of faces", elem)
		}
	}

	// Test mathematical properties
	testConnectivityProperties(t, EToE, EToF)

	// Count face types for cube mesh
	boundaryFaces := 0
	interiorFaces := 0

	for elem := 0; elem < 6; elem++ {
		for face := 0; face < 4; face++ {
			if EToE[elem][face] == elem {
				boundaryFaces++
			} else {
				interiorFaces++
			}
		}
	}

	// Mathematical property: A cube has 12 boundary faces (6 faces * 2 triangles per face)
	// and some number of interior faces
	if boundaryFaces != 12 {
		t.Errorf("Cube mesh should have 12 boundary faces, got %d", boundaryFaces)
	}

	t.Logf("Cube mesh statistics: %d boundary faces, %d interior face connections",
		boundaryFaces, interiorFaces)
}

// TestTiConnect3D_NodeNumbering tests robustness to different node numbering schemes
// Following Unit Testing Principle: Test degeneracies and special cases
func TestTiConnect3D_NodeNumbering(t *testing.T) {
	// Test with high node numbers
	EToV := [][]int{
		{100, 101, 102, 103},
		{101, 102, 103, 104}, // Shares face {101,102,103}
	}

	EToE, EToF := tiConnect3D(EToV)

	// Should find exactly one shared face between the two tets
	sharedCount := 0
	for elem := 0; elem < 2; elem++ {
		for face := 0; face < 4; face++ {
			if EToE[elem][face] != elem {
				sharedCount++
			}
		}
	}

	if sharedCount != 2 {
		t.Errorf("High node numbers: Expected 2 shared face connections, got %d", sharedCount)
	}

	// Test reciprocity
	testReciprocity(t, EToE, EToF)
}

// TestTiConnect3D_EdgeCases tests edge cases and error conditions
// Following Unit Testing Principle: Detailed Coverage
func TestTiConnect3D_EdgeCases(t *testing.T) {
	testCases := []struct {
		name     string
		EToV     [][]int
		validate func(t *testing.T, EToE, EToF [][]int)
	}{
		{
			name: "Empty mesh",
			EToV: [][]int{},
			validate: func(t *testing.T, EToE, EToF [][]int) {
				if len(EToE) != 0 || len(EToF) != 0 {
					t.Error("Empty mesh should return empty arrays")
				}
			},
		},
		{
			name: "Three tets sharing a vertex",
			EToV: [][]int{
				{0, 1, 2, 3},
				{0, 1, 2, 4},
				{0, 1, 3, 4},
			},
			validate: func(t *testing.T, EToE, EToF [][]int) {
				// Each tet shares one face with another
				for elem := 0; elem < 3; elem++ {
					sharedFaces := 0
					for face := 0; face < 4; face++ {
						if EToE[elem][face] != elem {
							sharedFaces++
						}
					}
					if sharedFaces == 0 {
						t.Errorf("Element %d should have at least one shared face", elem)
					}
				}
				testReciprocity(t, EToE, EToF)
			},
		},
		{
			name: "Four tets forming a pyramid",
			EToV: [][]int{
				{0, 1, 2, 4}, // Base tet 1
				{1, 2, 3, 4}, // Base tet 2
				{2, 3, 0, 4}, // Base tet 3
				{3, 0, 1, 4}, // Base tet 4
			},
			validate: func(t *testing.T, EToE, EToF [][]int) {
				testReciprocity(t, EToE, EToF)
				// Each tet should connect to exactly 3 others
				for elem := 0; elem < 4; elem++ {
					neighbors := make(map[int]bool)
					for face := 0; face < 4; face++ {
						if EToE[elem][face] != elem {
							neighbors[EToE[elem][face]] = true
						}
					}
					if len(neighbors) != 3 {
						t.Errorf("Element %d should connect to 3 neighbors, got %d",
							elem, len(neighbors))
					}
				}
			},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			EToE, EToF := tiConnect3D(tc.EToV)
			tc.validate(t, EToE, EToF)
		})
	}
}

// TestTiConnect3D_FaceOrdering verifies correct face numbering convention
// Following Unit Testing Principle: Complete basis testing
func TestTiConnect3D_FaceOrdering(t *testing.T) {
	// Create specific test case where we know which faces should match
	// Tet 0: (0,0,0), (1,0,0), (0,1,0), (0,0,1) - standard orientation
	// Tet 1: (1,0,0), (0,1,0), (0,0,1), (1,1,1) - shares face 2 of Tet 0
	EToV := [][]int{
		{0, 1, 2, 3}, // Tet 0
		{1, 2, 3, 4}, // Tet 1
	}

	EToE, EToF := tiConnect3D(EToV)

	// Find which face of Tet 0 connects to Tet 1
	sharedFace0 := -1
	for f := 0; f < 4; f++ {
		if EToE[0][f] == 1 {
			sharedFace0 = f
			break
		}
	}

	if sharedFace0 == -1 {
		t.Fatal("No shared face found between Tet 0 and Tet 1")
	}

	// The shared face should be face 2 of Tet 0 (vertices 1,2,3)
	// based on our face numbering convention
	expectedFace := 2
	if sharedFace0 != expectedFace {
		t.Errorf("Expected shared face to be face %d of Tet 0, got face %d",
			expectedFace, sharedFace0)
	}

	t.Logf("Tet 0 face %d connects to Tet 1 face %d",
		sharedFace0, EToF[0][sharedFace0])
}

// Helper function to test reciprocity property
func testReciprocity(t *testing.T, EToE, EToF [][]int) {
	K := len(EToE)
	for elem := 0; elem < K; elem++ {
		for face := 0; face < 4; face++ {
			neighbor := EToE[elem][face]
			if neighbor != elem {
				neighborFace := EToF[elem][face]

				// Check bounds
				if neighbor < 0 || neighbor >= K {
					t.Errorf("Invalid neighbor element %d for elem %d face %d",
						neighbor, elem, face)
					continue
				}

				if neighborFace < 0 || neighborFace >= 4 {
					t.Errorf("Invalid neighbor face %d for elem %d face %d",
						neighborFace, elem, face)
					continue
				}

				// Check reciprocity
				if EToE[neighbor][neighborFace] != elem {
					t.Errorf("Broken reciprocity: elem %d face %d -> elem %d, "+
						"but elem %d face %d -> elem %d",
						elem, face, neighbor,
						neighbor, neighborFace, EToE[neighbor][neighborFace])
				}

				if EToF[neighbor][neighborFace] != face {
					t.Errorf("Broken face reciprocity: elem %d face %d -> face %d, "+
						"but elem %d face %d -> face %d",
						elem, face, neighborFace,
						neighbor, neighborFace, EToF[neighbor][neighborFace])
				}
			}
		}
	}
}

// Helper function to test general connectivity properties
func testConnectivityProperties(t *testing.T, EToE, EToF [][]int) {
	K := len(EToE)

	// Property 1: Every element has exactly 4 faces
	for elem := 0; elem < K; elem++ {
		if len(EToE[elem]) != 4 || len(EToF[elem]) != 4 {
			t.Errorf("Element %d doesn'T have 4 faces", elem)
		}
	}

	// Property 2: Face indices are in valid range [0,3]
	for elem := 0; elem < K; elem++ {
		for face := 0; face < 4; face++ {
			if EToF[elem][face] < 0 || EToF[elem][face] >= 4 {
				t.Errorf("Invalid face index %d for elem %d face %d",
					EToF[elem][face], elem, face)
			}
		}
	}

	// Property 3: Reciprocity
	testReciprocity(t, EToE, EToF)

	// Property 4: No element connects to itself on interior faces
	for elem := 0; elem < K; elem++ {
		for face := 0; face < 4; face++ {
			if EToE[elem][face] == elem && EToF[elem][face] != face {
				t.Errorf("Element %d face %d self-connects to different face %d",
					elem, face, EToF[elem][face])
			}
		}
	}
}
