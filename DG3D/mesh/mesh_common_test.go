package mesh

import (
	"github.com/notargets/gocfd/utils"
	"testing"
)

// This test file verifies the face connectivity behavior in BuildConnectivity.
//
// The critical requirement is that EToF stores the neighbor's LOCAL face index,
// not a global face ID. This ensures reciprocal connectivity where:
// - If element A's face i connects to element B's face j
// - Then element B's face j connects to element A's face i
//
// The bug manifests as "reciprocal face connectivity failed" errors when
// EToF incorrectly stores global face IDs instead of local face indices.

// TestBuildConnectivity_BasicFunctionality tests the basic connectivity building
func TestBuildConnectivity_BasicFunctionality(t *testing.T) {
	// Create a simple two-tet mesh that shares a face
	mesh := &Mesh{
		Vertices: [][]float64{
			{0, 0, 0}, // 0
			{1, 0, 0}, // 1
			{0, 1, 0}, // 2
			{0, 0, 1}, // 3
			{1, 1, 1}, // 4
		},
		EtoV: [][]int{
			{0, 1, 2, 3}, // Tet 0
			{1, 2, 3, 4}, // Tet 1 - shares face {1,2,3} with Tet 0
		},
		ElementTypes: []utils.ElementType{
			utils.Tet,
			utils.Tet,
		},
		NumElements: 2,
		NumVertices: 5,
	}

	// Initialize required maps
	mesh.NodeIDMap = make(map[int]int)
	mesh.NodeArrayMap = make(map[int]int)
	mesh.ElementIDMap = make(map[int]int)
	mesh.FaceMap = make(map[string]int)

	for i := 0; i < 5; i++ {
		mesh.NodeIDMap[i] = i
		mesh.NodeArrayMap[i] = i
	}
	for i := 0; i < 2; i++ {
		mesh.ElementIDMap[i] = i
	}

	// Build connectivity
	mesh.BuildConnectivity()

	// Test 1: Check that arrays are initialized
	if mesh.EToE == nil {
		t.Fatal("EToE not initialized")
	}
	if mesh.EToF == nil {
		t.Fatal("EToF not initialized")
	}
	if len(mesh.EToE) != 2 {
		t.Errorf("Expected 2 elements in EToE, got %d", len(mesh.EToE))
	}
	if len(mesh.EToF) != 2 {
		t.Errorf("Expected 2 elements in EToF, got %d", len(mesh.EToF))
	}

	// Test 2: Each tet should have 4 faces
	for elem := 0; elem < 2; elem++ {
		if len(mesh.EToE[elem]) != 4 {
			t.Errorf("Element %d: expected 4 faces in EToE, got %d", elem, len(mesh.EToE[elem]))
		}
		if len(mesh.EToF[elem]) != 4 {
			t.Errorf("Element %d: expected 4 faces in EToF, got %d", elem, len(mesh.EToF[elem]))
		}
	}

	// Test 3: Find the shared face and verify reciprocal connectivity
	sharedFaceFound := false
	for f0 := 0; f0 < 4; f0++ {
		if mesh.EToE[0][f0] == 1 {
			// Element 0, face f0 connects to element 1
			f1 := mesh.EToF[0][f0]

			// Critical assertion: f1 must be a valid local face index for tet
			if f1 < 0 || f1 >= 4 {
				t.Errorf("Invalid face index: Element 0 face %d connects to Element 1 face %d (should be 0-3)", f0, f1)
			}

			// Verify reciprocal connectivity
			if mesh.EToE[1][f1] != 0 {
				t.Errorf("Reciprocal connectivity failed: Element 1 face %d should connect to Element 0, got Element %d",
					f1, mesh.EToE[1][f1])
			}
			if mesh.EToF[1][f1] != f0 {
				t.Errorf("Reciprocal face connectivity failed: Element 1 face %d should connect to Element 0 face %d, got face %d",
					f1, f0, mesh.EToF[1][f1])
			}

			sharedFaceFound = true
			t.Logf("Found shared face: Element 0 face %d <-> Element 1 face %d", f0, f1)
		}
	}

	if !sharedFaceFound {
		t.Error("No shared face found between the two tetrahedra")
	}

	// Test 4: Verify boundary faces
	boundaryCount := 0
	for elem := 0; elem < 2; elem++ {
		for face := 0; face < 4; face++ {
			if mesh.EToE[elem][face] == -1 {
				boundaryCount++
				if mesh.EToF[elem][face] != -1 {
					t.Errorf("Boundary face should have EToF = -1, got %d", mesh.EToF[elem][face])
				}
			}
		}
	}

	// Two tets sharing one face should have 6 boundary faces total
	expectedBoundaryFaces := 6
	if boundaryCount != expectedBoundaryFaces {
		t.Errorf("Expected %d boundary faces, got %d", expectedBoundaryFaces, boundaryCount)
	}
}

// TestBuildConnectivity_GetElementFaces tests the face vertex ordering
func TestBuildConnectivity_GetElementFaces(t *testing.T) {
	// Test tetrahedron face ordering
	tetVerts := []int{0, 1, 2, 3}
	tetFaces := GetElementFaces(utils.Tet, tetVerts)

	if len(tetFaces) != 4 {
		t.Errorf("Tet should have 4 faces, got %d", len(tetFaces))
	}

	// Each tet face should have 3 vertices
	for i, face := range tetFaces {
		if len(face) != 3 {
			t.Errorf("Tet face %d should have 3 vertices, got %d", i, len(face))
		}
	}

	// Test that the faces form valid triangles from the tet vertices
	expectedFaces := [][]int{
		{0, 2, 1}, // Face 0
		{0, 1, 3}, // Face 1
		{0, 3, 2}, // Face 2
		{1, 2, 3}, // Face 3
	}

	for i, expected := range expectedFaces {
		actual := tetFaces[i]
		if !equalSlices(actual, expected) {
			t.Errorf("Tet face %d: expected %v, got %v", i, expected, actual)
		}
	}
}

// Helper function to compare integer slices
func equalSlices(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
