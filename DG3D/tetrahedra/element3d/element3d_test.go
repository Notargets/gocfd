package element3d

import (
	"testing"
)

// TestConnect3D_Reciprocity tests that Connect3D creates reciprocal connectivity
func TestConnect3D_Reciprocity(t *testing.T) {
	// Create a simple two-tet mesh
	el := &Element3D{
		K: 2,
		EToV: [][]int{
			{0, 1, 2, 3}, // Tet 0
			{1, 2, 3, 4}, // Tet 1 - shares face {1,2,3} with Tet 0
		},
	}

	// Build connectivity
	ca := el.Connect3D()
	el.ConnectivityArrays = ca

	// Verify arrays are created
	if len(ca.EToE) != 2 || len(ca.EToF) != 2 {
		t.Fatalf("Expected 2 elements in connectivity arrays")
	}

	// Find the shared face
	sharedFaceFound := false
	for f0 := 0; f0 < 4; f0++ {
		if ca.EToE[0][f0] == 1 {
			f1 := ca.EToF[0][f0]

			// Critical checks
			if f1 < 0 || f1 >= 4 {
				t.Errorf("Invalid face index: Tet 0 face %d connects to Tet 1 face %d (should be 0-3)", f0, f1)
			}

			// Verify reciprocal connectivity
			if ca.EToE[1][f1] != 0 {
				t.Errorf("Reciprocal connectivity failed: Tet 1 face %d should connect to Tet 0, got %d",
					f1, ca.EToE[1][f1])
			}

			if ca.EToF[1][f1] != f0 {
				t.Errorf("Reciprocal face connectivity failed: Tet 1 face %d should connect to Tet 0 face %d, got %d",
					f1, f0, ca.EToF[1][f1])
			}

			sharedFaceFound = true
			t.Logf("Shared face found: Tet 0 face %d <-> Tet 1 face %d", f0, f1)
		}
	}

	if !sharedFaceFound {
		t.Error("No shared face found between the two tetrahedra")
	}

	// Count boundary faces
	boundaryCount := 0
	for elem := 0; elem < 2; elem++ {
		for face := 0; face < 4; face++ {
			if ca.EToE[elem][face] == elem && ca.EToF[elem][face] == face {
				boundaryCount++
			}
		}
	}

	// Two tets sharing one face should have 6 boundary faces
	if boundaryCount != 6 {
		t.Errorf("Expected 6 boundary faces, got %d", boundaryCount)
	}
}

// TestConnect3D_ComplexMesh tests connectivity with a larger mesh
func TestConnect3D_ComplexMesh(t *testing.T) {
	// Create a mesh with 4 tetrahedra
	el := &Element3D{
		K: 4,
		EToV: [][]int{
			{0, 1, 2, 3}, // Tet 0
			{1, 4, 2, 5}, // Tet 1
			{1, 2, 3, 5}, // Tet 2 - shares face {1,2,3} with Tet 0
			{2, 3, 5, 6}, // Tet 3 - shares face {2,3,5} with Tet 2
		},
	}

	// Build connectivity
	ca := el.Connect3D()

	// Verify all reciprocal connections
	failures := 0
	for elem := 0; elem < el.K; elem++ {
		for face := 0; face < 4; face++ {
			neighbor := ca.EToE[elem][face]

			// Skip boundary faces
			if neighbor == elem {
				continue
			}

			neighborFace := ca.EToF[elem][face]

			// Check bounds
			if neighborFace < 0 || neighborFace >= 4 {
				t.Errorf("Invalid face index: Element %d face %d -> Element %d face %d",
					elem, face, neighbor, neighborFace)
				failures++
				continue
			}

			// Check reciprocity
			if ca.EToE[neighbor][neighborFace] != elem {
				t.Errorf("Element %d face %d -> Element %d, but Element %d face %d -> Element %d",
					elem, face, neighbor, neighbor, neighborFace, ca.EToE[neighbor][neighborFace])
				failures++
			}

			if ca.EToF[neighbor][neighborFace] != face {
				t.Errorf("Element %d face %d: reciprocal face connectivity failed (expected face %d, got %d)",
					elem, face, face, ca.EToF[neighbor][neighborFace])
				failures++
			}
		}
	}

	if failures > 0 {
		t.Errorf("Found %d reciprocal connectivity failures", failures)
	} else {
		t.Log("All reciprocal connectivity verified successfully")
	}
}

// TestConnect3D_FaceVertexMapping verifies the face vertex ordering
func TestConnect3D_FaceVertexMapping(t *testing.T) {
	// The face vertex mapping should match what's used in the code
	expectedFaceVertices := [][]int{
		{0, 1, 2}, // Face 0
		{0, 1, 3}, // Face 1
		{1, 2, 3}, // Face 2
		{0, 2, 3}, // Face 3
	}

	// Create a single tet to verify face ordering
	el := &Element3D{
		K: 1,
		EToV: [][]int{
			{10, 11, 12, 13}, // Use distinct vertex IDs
		},
	}

	ca := el.Connect3D()

	// All faces should be boundary (self-connected)
	for f := 0; f < 4; f++ {
		if ca.EToE[0][f] != 0 || ca.EToF[0][f] != f {
			t.Errorf("Face %d should be boundary (self-connected)", f)
		}
	}

	// Log the expected face vertices for documentation
	t.Log("Face vertex mappings:")
	for f, verts := range expectedFaceVertices {
		actualVerts := make([]int, 3)
		for i := 0; i < 3; i++ {
			actualVerts[i] = el.EToV[0][verts[i]]
		}
		t.Logf("  Face %d: local vertices %v -> global vertices %v", f, verts, actualVerts)
	}
}
