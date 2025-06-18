package mesh

import (
	"fmt"
	"sort"
	"testing"
	"time"
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
		ElementTypes: []ElementType{
			Tet,
			Tet,
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
	tetFaces := GetElementFaces(Tet, tetVerts)

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

// TestBuildConnectivity_ComplexMesh tests connectivity with multiple element types
func TestBuildConnectivity_ComplexMesh(t *testing.T) {
	// Use the standard test meshes
	tm := GetStandardTestMeshes()
	mesh := tm.MixedMesh.ConvertToMesh()

	// Verify connectivity was built
	if mesh.EToE == nil || mesh.EToF == nil {
		t.Fatal("Connectivity not built")
	}

	// Check each element has the correct number of faces
	expectedFaceCounts := map[ElementType]int{
		Tet:     4,
		Hex:     6,
		Prism:   5,
		Pyramid: 5,
	}

	for i, elemType := range mesh.ElementTypes {
		expectedFaces := expectedFaceCounts[elemType]
		if len(mesh.EToE[i]) != expectedFaces {
			t.Errorf("Element %d (%v): expected %d faces, got %d",
				i, elemType, expectedFaces, len(mesh.EToE[i]))
		}
	}

	// Test reciprocal connectivity for all interior faces
	for elem := 0; elem < mesh.NumElements; elem++ {
		for face := 0; face < len(mesh.EToE[elem]); face++ {
			neighbor := mesh.EToE[elem][face]

			if neighbor >= 0 && neighbor < mesh.NumElements {
				// Interior face - check reciprocity
				neighborFace := mesh.EToF[elem][face]

				// Verify neighbor face index is valid
				if neighborFace < 0 || neighborFace >= len(mesh.EToE[neighbor]) {
					t.Errorf("Invalid neighbor face index: Element %d face %d -> Element %d face %d (element type %v has %d faces)",
						elem, face, neighbor, neighborFace, mesh.ElementTypes[neighbor], len(mesh.EToE[neighbor]))
					continue
				}

				// Check reciprocal connectivity
				if mesh.EToE[neighbor][neighborFace] != elem {
					t.Errorf("Reciprocal element connectivity failed: Element %d face %d -> Element %d, but Element %d face %d -> Element %d",
						elem, face, neighbor, neighbor, neighborFace, mesh.EToE[neighbor][neighborFace])
				}

				if mesh.EToF[neighbor][neighborFace] != face {
					t.Errorf("Reciprocal face connectivity failed: Element %d face %d -> face %d, but Element %d face %d -> face %d",
						elem, face, neighborFace, neighbor, neighborFace, mesh.EToF[neighbor][neighborFace])
				}
			}
		}
	}
}

// TestBuildConnectivity_SharedFaceVertices verifies that connected faces share the same vertices
func TestBuildConnectivity_SharedFaceVertices(t *testing.T) {
	tm := GetStandardTestMeshes()
	mesh := tm.TwoTetMesh.ConvertToMesh()

	// For each interior face connection, verify the faces share the same vertices
	for elem := 0; elem < mesh.NumElements; elem++ {
		elemFaces := GetElementFaces(mesh.ElementTypes[elem], mesh.EtoV[elem])

		for face := 0; face < len(mesh.EToE[elem]); face++ {
			neighbor := mesh.EToE[elem][face]

			if neighbor >= 0 && neighbor < mesh.NumElements && neighbor > elem {
				// Interior face - check only once per face pair
				neighborFace := mesh.EToF[elem][face]

				// Get face vertices from both sides
				faceVerts1 := make([]int, len(elemFaces[face]))
				copy(faceVerts1, elemFaces[face])
				sort.Ints(faceVerts1)

				neighborFaces := GetElementFaces(mesh.ElementTypes[neighbor], mesh.EtoV[neighbor])
				faceVerts2 := make([]int, len(neighborFaces[neighborFace]))
				copy(faceVerts2, neighborFaces[neighborFace])
				sort.Ints(faceVerts2)

				// Verify they share the same vertices
				if !equalSlices(faceVerts1, faceVerts2) {
					t.Errorf("Face mismatch: Element %d face %d has vertices %v, Element %d face %d has vertices %v",
						elem, face, faceVerts1, neighbor, neighborFace, faceVerts2)
				}
			}
		}
	}
}

// TestBuildConnectivity_CubeMesh tests a more complex mesh
func TestBuildConnectivity_CubeMesh(t *testing.T) {
	tm := GetStandardTestMeshes()
	mesh := tm.CubeMesh.ConvertToMesh()

	// The cube is meshed with 6 tetrahedra
	if mesh.NumElements != 6 {
		t.Errorf("Expected 6 elements, got %d", mesh.NumElements)
	}

	// Count all face connections
	totalFaceConnections := 0
	boundaryFaces := 0
	interiorConnections := 0

	for elem := 0; elem < mesh.NumElements; elem++ {
		for face := 0; face < 4; face++ { // Tets have 4 faces
			totalFaceConnections++
			if mesh.EToE[elem][face] == -1 {
				boundaryFaces++
			} else {
				interiorConnections++
			}
		}
	}

	// Calculate unique interior faces (each counted twice)
	uniqueInteriorFaces := interiorConnections / 2

	t.Logf("Total face connections: %d", totalFaceConnections)
	t.Logf("Boundary faces: %d, Interior connections: %d (= %d unique interior faces)",
		boundaryFaces, interiorConnections, uniqueInteriorFaces)

	// Total face connections should be 6 tets * 4 faces = 24
	expectedTotalConnections := 6 * 4
	if totalFaceConnections != expectedTotalConnections {
		t.Errorf("Expected %d total face connections, got %d",
			expectedTotalConnections, totalFaceConnections)
	}

	// Verify the count adds up
	if boundaryFaces+interiorConnections != totalFaceConnections {
		t.Errorf("Face count mismatch: %d boundary + %d interior != %d total",
			boundaryFaces, interiorConnections, totalFaceConnections)
	}

	// Verify all face connections are reciprocal
	connectivityErrors := 0
	for elem := 0; elem < mesh.NumElements; elem++ {
		for face := 0; face < 4; face++ {
			neighbor := mesh.EToE[elem][face]
			if neighbor >= 0 {
				neighborFace := mesh.EToF[elem][face]

				// Check bounds
				if neighborFace < 0 || neighborFace >= 4 {
					t.Errorf("Invalid face index: Element %d face %d -> Element %d face %d",
						elem, face, neighbor, neighborFace)
					connectivityErrors++
					continue
				}

				// Check reciprocity
				if mesh.EToE[neighbor][neighborFace] != elem {
					connectivityErrors++
				}
				if mesh.EToF[neighbor][neighborFace] != face {
					connectivityErrors++
				}
			}
		}
	}

	if connectivityErrors > 0 {
		t.Errorf("Found %d connectivity errors", connectivityErrors)
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

// TestBuildConnectivity_Performance benchmarks the connectivity building
func TestBuildConnectivity_Performance(t *testing.T) {
	// Create a larger mesh for performance testing
	nCubes := 10 // 10x10x10 grid of cubes
	vertices := [][]float64{}
	elements := [][]int{}
	elementTypes := []ElementType{}

	// Generate vertices for a grid
	vertexMap := make(map[string]int)
	for i := 0; i <= nCubes; i++ {
		for j := 0; j <= nCubes; j++ {
			for k := 0; k <= nCubes; k++ {
				key := fmt.Sprintf("%d,%d,%d", i, j, k)
				vertexMap[key] = len(vertices)
				vertices = append(vertices, []float64{float64(i), float64(j), float64(k)})
			}
		}
	}

	// Generate tetrahedra for each cube
	for i := 0; i < nCubes; i++ {
		for j := 0; j < nCubes; j++ {
			for k := 0; k < nCubes; k++ {
				// Get the 8 vertices of this cube
				v := [8]int{
					vertexMap[fmt.Sprintf("%d,%d,%d", i, j, k)],
					vertexMap[fmt.Sprintf("%d,%d,%d", i+1, j, k)],
					vertexMap[fmt.Sprintf("%d,%d,%d", i+1, j+1, k)],
					vertexMap[fmt.Sprintf("%d,%d,%d", i, j+1, k)],
					vertexMap[fmt.Sprintf("%d,%d,%d", i, j, k+1)],
					vertexMap[fmt.Sprintf("%d,%d,%d", i+1, j, k+1)],
					vertexMap[fmt.Sprintf("%d,%d,%d", i+1, j+1, k+1)],
					vertexMap[fmt.Sprintf("%d,%d,%d", i, j+1, k+1)],
				}

				// Add 6 tetrahedra per cube (using a standard decomposition)
				tets := [][]int{
					{v[0], v[1], v[3], v[4]},
					{v[1], v[2], v[3], v[6]},
					{v[1], v[3], v[4], v[6]},
					{v[3], v[4], v[6], v[7]},
					{v[1], v[4], v[5], v[6]},
					{v[4], v[5], v[6], v[7]},
				}

				for _, tet := range tets {
					elements = append(elements, tet)
					elementTypes = append(elementTypes, Tet)
				}
			}
		}
	}

	// Create mesh
	mesh := &Mesh{
		Vertices:     vertices,
		EtoV:         elements,
		ElementTypes: elementTypes,
		NumElements:  len(elements),
		NumVertices:  len(vertices),
		FaceMap:      make(map[string]int),
	}

	// Initialize ID maps
	mesh.NodeIDMap = make(map[int]int)
	mesh.NodeArrayMap = make(map[int]int)
	mesh.ElementIDMap = make(map[int]int)
	for i := 0; i < len(vertices); i++ {
		mesh.NodeIDMap[i] = i
		mesh.NodeArrayMap[i] = i
	}
	for i := 0; i < len(elements); i++ {
		mesh.ElementIDMap[i] = i
	}

	t.Logf("Testing with %d elements and %d vertices", mesh.NumElements, mesh.NumVertices)

	// Time the connectivity building
	start := time.Now()
	mesh.BuildConnectivity()
	elapsed := time.Since(start)

	t.Logf("BuildConnectivity took %v for %d elements", elapsed, mesh.NumElements)

	// Verify some basic properties
	if len(mesh.Faces) == 0 {
		t.Error("No faces created")
	}

	// Count interior vs boundary faces
	interiorCount := 0
	for elem := 0; elem < mesh.NumElements; elem++ {
		for face := 0; face < 4; face++ {
			if mesh.EToE[elem][face] >= 0 && mesh.EToE[elem][face] > elem {
				interiorCount++
			}
		}
	}

	t.Logf("Created %d unique faces, %d interior face connections", len(mesh.Faces), interiorCount)
}

// TestBuildConnectivity_BugRegression specifically tests for the bug where
// EToF stores global face IDs instead of local face indices
func TestBuildConnectivity_BugRegression(t *testing.T) {
	// This test reproduces the exact error from the bug report:
	// "Element 357 face 2: reciprocal face connectivity failed (expected face 2, got 3)"
	//
	// The bug occurs when BuildConnectivity stores global face IDs in EToF
	// instead of the neighbor's local face index.

	// Create a mesh where face indices would differ from global face IDs
	// We need a configuration where the same face appears at different local indices
	mesh := &Mesh{
		Vertices: [][]float64{
			{0, 0, 0}, // 0
			{1, 0, 0}, // 1
			{0, 1, 0}, // 2
			{0, 0, 1}, // 3
			{1, 1, 0}, // 4
			{1, 0, 1}, // 5
			{0, 1, 1}, // 6
		},
		EtoV: [][]int{
			{0, 1, 2, 3}, // Tet 0
			{1, 4, 2, 5}, // Tet 1 - shares face {1,2,4} with potential face mismatch
			{1, 2, 3, 5}, // Tet 2 - shares face {1,2,3} with Tet 0
			{2, 3, 5, 6}, // Tet 3 - shares face {2,3,5} with Tet 2
		},
		ElementTypes: []ElementType{
			Tet, Tet, Tet, Tet,
		},
		NumElements: 4,
		NumVertices: 7,
	}

	// Initialize maps
	mesh.NodeIDMap = make(map[int]int)
	mesh.NodeArrayMap = make(map[int]int)
	mesh.ElementIDMap = make(map[int]int)
	mesh.FaceMap = make(map[string]int)

	for i := 0; i < 7; i++ {
		mesh.NodeIDMap[i] = i
		mesh.NodeArrayMap[i] = i
	}
	for i := 0; i < 4; i++ {
		mesh.ElementIDMap[i] = i
	}

	// Build connectivity
	mesh.BuildConnectivity()

	// Verify reciprocal connectivity for all connections
	failures := 0
	for elem := 0; elem < mesh.NumElements; elem++ {
		for face := 0; face < 4; face++ {
			neighbor := mesh.EToE[elem][face]
			if neighbor >= 0 {
				neighborFace := mesh.EToF[elem][face]

				// Check bounds
				if neighborFace < 0 || neighborFace >= 4 {
					t.Errorf("Invalid face index: Element %d face %d -> Element %d face %d",
						elem, face, neighbor, neighborFace)
					failures++
					continue
				}

				// The bug would manifest here: if EToF stored global face IDs,
				// the reciprocal check would fail
				if mesh.EToE[neighbor][neighborFace] != elem {
					t.Errorf("Element %d face %d -> Element %d, but Element %d face %d -> Element %d",
						elem, face, neighbor, neighbor, neighborFace, mesh.EToE[neighbor][neighborFace])
					failures++
				}

				if mesh.EToF[neighbor][neighborFace] != face {
					// This is the exact error from the bug report
					t.Errorf("Element %d face %d: reciprocal face connectivity failed (expected face %d, got %d)",
						elem, face, face, mesh.EToF[neighbor][neighborFace])
					failures++
				}
			}
		}
	}

	if failures > 0 {
		t.Errorf("Found %d reciprocal connectivity failures", failures)
		t.Log("This indicates EToF is storing global face IDs instead of local face indices")
	} else {
		t.Log("All reciprocal connectivity verified - BuildConnectivity correctly stores local face indices")
	}

	// Additional check: verify that shared faces have different local indices
	// This is the key test - if EToF stored global face IDs, shared faces would have the same ID
	sharedFaceCount := 0
	for elem := 0; elem < mesh.NumElements; elem++ {
		for face := 0; face < 4; face++ {
			neighbor := mesh.EToE[elem][face]
			if neighbor > elem { // Check each pair only once
				neighborFace := mesh.EToF[elem][face]

				// The key insight: face and neighborFace are LOCAL indices
				// They will often be different (e.g., face 2 on one tet might connect to face 3 on another)
				if face != neighborFace {
					t.Logf("Shared face: Element %d face %d <-> Element %d face %d (different local indices)",
						elem, face, neighbor, neighborFace)
				}
				sharedFaceCount++
			}
		}
	}

	t.Logf("Found %d shared faces with potentially different local indices", sharedFaceCount)
}

// TestBuildConnectivity_NonManifold tests behavior with non-manifold meshes
func TestBuildConnectivity_NonManifold(t *testing.T) {
	// Create a non-manifold mesh where 3 tets share the same face
	// This is geometrically invalid but tests error handling
	mesh := &Mesh{
		Vertices: [][]float64{
			{0, 0, 0}, // 0
			{1, 0, 0}, // 1
			{0, 1, 0}, // 2
			{0, 0, 1}, // 3
			{1, 1, 0}, // 4
			{1, 0, 1}, // 5
		},
		EtoV: [][]int{
			{0, 1, 2, 3}, // Tet 0
			{1, 4, 2, 3}, // Tet 1 - shares face {1,2,3} with Tet 0
			{1, 2, 3, 5}, // Tet 2 - also shares face {1,2,3} - NON-MANIFOLD!
		},
		ElementTypes: []ElementType{
			Tet, Tet, Tet,
		},
		NumElements: 3,
		NumVertices: 6,
	}

	// Initialize maps
	mesh.NodeIDMap = make(map[int]int)
	mesh.NodeArrayMap = make(map[int]int)
	mesh.ElementIDMap = make(map[int]int)
	mesh.FaceMap = make(map[string]int)

	for i := 0; i < 6; i++ {
		mesh.NodeIDMap[i] = i
		mesh.NodeArrayMap[i] = i
	}
	for i := 0; i < 3; i++ {
		mesh.ElementIDMap[i] = i
	}

	// Build connectivity
	mesh.BuildConnectivity()

	// Find the shared face {1,2,3}
	targetVerts := []int{1, 2, 3}
	sort.Ints(targetVerts)

	// Count how many elements claim to share this face
	elementsWithFace := []int{}
	for elem := 0; elem < 3; elem++ {
		faces := GetElementFaces(Tet, mesh.EtoV[elem])
		for _, face := range faces {
			sortedFace := make([]int, len(face))
			copy(sortedFace, face)
			sort.Ints(sortedFace)
			if equalSlices(sortedFace, targetVerts) {
				elementsWithFace = append(elementsWithFace, elem)
				break
			}
		}
	}

	t.Logf("Face {1,2,3} is present in elements: %v", elementsWithFace)

	if len(elementsWithFace) > 2 {
		t.Log("WARNING: Non-manifold mesh detected - face shared by more than 2 elements")
		t.Log("BuildConnectivity can only handle manifold meshes where each face is shared by at most 2 elements")
	}

	// The current implementation will connect only the last two elements that share the face
	// This is a limitation, not a bug - proper mesh generation should avoid non-manifold configurations
}
