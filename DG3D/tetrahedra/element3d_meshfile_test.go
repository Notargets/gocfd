package tetrahedra

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"testing"

	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/DG3D/mesh/readers"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// getTestMeshPath returns the path to the test mesh file
func getTestMeshPath() string {
	// Try multiple possible paths
	possiblePaths := []string{
		"../mesh/cube-partitioned.neu",
		"../../mesh/cube-partitioned.neu",
		"DG3D/mesh/cube-partitioned.neu",
		"mesh/cube-partitioned.neu",
	}

	// Also try relative to the current file location
	_, filename, _, ok := runtime.Caller(0)
	if ok {
		dir := filepath.Dir(filename)
		possiblePaths = append(possiblePaths,
			filepath.Join(dir, "..", "mesh", "cube-partitioned.neu"),
			filepath.Join(dir, "..", "..", "mesh", "cube-partitioned.neu"),
		)
	}

	// Check each path
	for _, path := range possiblePaths {
		if _, err := os.Stat(path); err == nil {
			return path
		}
	}

	// If not found, return the expected path for error messages
	return "../mesh/cube-partitioned.neu"
}

func TestElement3D_MeshConnectivity(t *testing.T) {
	meshPath := getTestMeshPath()

	// First, let's debug the mesh connectivity before Element3D processes it
	t.Run("Debug_Mesh_Connectivity", func(t *testing.T) {
		// Load the mesh directly
		m, err := readers.ReadMeshFile(meshPath)
		require.NoError(t, err)

		// Build connectivity
		m.BuildConnectivity()

		// Check for invalid face indices
		invalidFound := false
		for elemID := 0; elemID < m.NumElements; elemID++ {
			if m.ElementTypes[elemID] == mesh.Tet {
				for faceID := 0; faceID < len(m.EToF[elemID]); faceID++ {
					neighborFace := m.EToF[elemID][faceID]
					neighborElem := m.EToE[elemID][faceID]

					// Skip boundary faces
					if neighborElem == -1 {
						continue
					}

					// Check if neighbor is also a tet
					if neighborElem < m.NumElements && m.ElementTypes[neighborElem] == mesh.Tet {
						if neighborFace >= 4 {
							t.Errorf("Invalid face index: Element %d (tet) face %d connects to Element %d (tet) face %d (>= 4)",
								elemID, faceID, neighborElem, neighborFace)
							invalidFound = true
						}
					}
				}
			}
		}

		if invalidFound {
			t.Log("Found invalid face indices in mesh connectivity")
		} else {
			t.Log("All face indices in mesh connectivity are valid")
		}
	})

	// Test with different polynomial orders
	orders := []int{1, 2, 3}

	for _, order := range orders {
		t.Run(fmt.Sprintf("Order_%d", order), func(t *testing.T) {
			// Create Element3D from mesh
			el, err := NewElement3D(order, meshPath)
			require.NoError(t, err, "Failed to create Element3D from mesh")
			require.NotNil(t, el)

			// Basic validation
			assert.Greater(t, el.K, 0, "Should have elements")
			assert.NotNil(t, el.ConnectivityArrays, "Should have connectivity")
			assert.NotNil(t, el.EToE, "Should have EToE")
			assert.NotNil(t, el.EToF, "Should have EToF")

			// Validate EToE and EToF dimensions
			assert.Equal(t, el.K, len(el.EToE), "EToE should have K elements")
			assert.Equal(t, el.K, len(el.EToF), "EToF should have K elements")

			// Validate each element's connectivity
			for k := 0; k < el.K; k++ {
				// Each tet should have exactly 4 faces
				assert.Equal(t, 4, len(el.EToE[k]), "Element %d: EToE should have 4 faces", k)
				assert.Equal(t, 4, len(el.EToF[k]), "Element %d: EToF should have 4 faces", k)

				// Validate face indices
				for f := 0; f < 4; f++ {
					neighborElem := el.EToE[k][f]
					neighborFace := el.EToF[k][f]

					// Check if boundary face
					if neighborElem == -1 {
						assert.Equal(t, -1, neighborFace,
							"Element %d face %d: boundary face should have neighborFace = -1", k, f)
					} else {
						// Interior face
						assert.GreaterOrEqual(t, neighborElem, 0,
							"Element %d face %d: neighbor element should be >= 0", k, f)
						assert.Less(t, neighborElem, el.K,
							"Element %d face %d: neighbor element %d should be < K=%d", k, f, neighborElem, el.K)

						// CRITICAL: Face index must be 0-3 for tetrahedra
						assert.GreaterOrEqual(t, neighborFace, 0,
							"Element %d face %d: neighbor face should be >= 0", k, f)
						assert.Less(t, neighborFace, 4,
							"Element %d face %d: neighbor face %d should be < 4 (tet has 4 faces)", k, f, neighborFace)

						// Verify reciprocal connectivity
						if neighborElem < el.K && neighborFace < 4 {
							reciprocalElem := el.EToE[neighborElem][neighborFace]
							reciprocalFace := el.EToF[neighborElem][neighborFace]

							assert.Equal(t, k, reciprocalElem,
								"Element %d face %d: reciprocal connectivity failed (expected elem %d, got %d)",
								k, f, k, reciprocalElem)
							assert.Equal(t, f, reciprocalFace,
								"Element %d face %d: reciprocal face connectivity failed (expected face %d, got %d)",
								k, f, f, reciprocalFace)
						}
					}
				}
			}
		})
	}
}

func TestElement3D_BuildMaps3D_WithRealMesh(t *testing.T) {
	meshPath := getTestMeshPath()

	// Create Element3D with order 1 (simplest case)
	el, err := NewElement3D(1, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el)

	// BuildMaps3D should have been called in constructor
	// Verify it completed without panic
	assert.NotNil(t, el.VmapM, "VmapM should be initialized")
	assert.NotNil(t, el.VmapP, "VmapP should be initialized")
	assert.NotNil(t, el.MapM, "MapM should be initialized")
	assert.NotNil(t, el.MapP, "MapP should be initialized")

	// Verify array sizes
	expectedSize := el.Nfp * 4 * el.K // Nfp nodes per face * 4 faces * K elements
	assert.Equal(t, expectedSize, len(el.VmapM), "VmapM size mismatch")
	assert.Equal(t, expectedSize, len(el.VmapP), "VmapP size mismatch")
	assert.Equal(t, expectedSize, len(el.MapM), "MapM size mismatch")
	assert.Equal(t, expectedSize, len(el.MapP), "MapP size mismatch")
}

func TestElement3D_FaceConnectivityConsistency(t *testing.T) {
	meshPath := getTestMeshPath()

	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)

	// Track face pairs to ensure each internal face is counted exactly twice
	facePairs := make(map[string]int)

	for k := 0; k < el.K; k++ {
		for f := 0; f < 4; f++ {
			neighbor := el.EToE[k][f]

			if neighbor != -1 && neighbor != k {
				// Create a unique key for this face pair
				var key string
				if k < neighbor {
					key = fmt.Sprintf("%d_%d", k, neighbor)
				} else {
					key = fmt.Sprintf("%d_%d", neighbor, k)
				}
				facePairs[key]++
			}
		}
	}

	// Each internal face should be counted exactly twice (once from each side)
	for key, count := range facePairs {
		assert.Equal(t, 2, count, "Face pair %s should be counted exactly twice", key)
	}
}

func TestElement3D_PartitionData(t *testing.T) {
	meshPath := getTestMeshPath()

	el, err := NewElement3D(1, meshPath)
	require.NoError(t, err)

	// If mesh has partition data, verify it's loaded
	if el.Mesh.EToP != nil {
		assert.Equal(t, el.Mesh.NumElements, len(el.Mesh.EToP),
			"EToP should have entry for each element")

		// Check partition IDs are reasonable
		maxPartID := -1
		for _, partID := range el.Mesh.EToP {
			assert.GreaterOrEqual(t, partID, 0, "Partition ID should be non-negative")
			if partID > maxPartID {
				maxPartID = partID
			}
		}

		t.Logf("Mesh has %d partitions (0 to %d)", maxPartID+1, maxPartID)
	} else {
		t.Log("Mesh has no partition data")
	}
}

func TestMeshConnectivityBug_Regression(t *testing.T) {
	// This is a regression test for the bug where BuildConnectivity stores global face IDs
	// in EToF instead of local face IDs. This test should FAIL with the buggy implementation
	// and PASS after the fix is applied.

	// Create a simple mesh with 2 tetrahedra sharing a face
	m := &mesh.Mesh{
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
		ElementTypes: []mesh.ElementType{
			mesh.Tet,
			mesh.Tet,
		},
		NumElements: 2,
		NumVertices: 5,
	}

	// Initialize maps
	m.NodeIDMap = make(map[int]int)
	m.NodeArrayMap = make(map[int]int)
	m.ElementIDMap = make(map[int]int)
	m.FaceMap = make(map[string]int)

	for i := 0; i < 5; i++ {
		m.NodeIDMap[i] = i
		m.NodeArrayMap[i] = i
	}
	for i := 0; i < 2; i++ {
		m.ElementIDMap[i] = i
	}

	// Build connectivity
	m.BuildConnectivity()

	// Check connectivity
	require.NotNil(t, m.EToE)
	require.NotNil(t, m.EToF)
	require.Equal(t, 2, len(m.EToE))
	require.Equal(t, 2, len(m.EToF))

	// Find the shared face between the two tets
	sharedFaceFound := false

	for f0 := 0; f0 < 4; f0++ {
		if m.EToE[0][f0] == 1 {
			// Element 0 face f0 connects to element 1
			f1 := m.EToF[0][f0]

			// CRITICAL ASSERTION: f1 must be a valid local face index (0-3) for tet 1
			assert.GreaterOrEqual(t, f1, 0,
				"Tet 0 face %d connects to Tet 1: neighbor face index must be >= 0", f0)
			assert.Less(t, f1, 4,
				"Tet 0 face %d connects to Tet 1: neighbor face index %d must be < 4", f0, f1)

			// Verify reciprocal connectivity
			assert.Equal(t, 0, m.EToE[1][f1],
				"Tet 1 face %d should connect back to Tet 0", f1)
			assert.Equal(t, f0, m.EToF[1][f1],
				"Tet 1 face %d should connect back to Tet 0 face %d", f1, f0)

			sharedFaceFound = true
			t.Logf("Shared face: Tet 0 face %d <-> Tet 1 face %d", f0, f1)
		}
	}

	assert.True(t, sharedFaceFound, "Should find the shared face between the two tets")
}

func TestElement3D_ReadsPartitionData(t *testing.T) {
	// Test that Element3D correctly reads partition data from cube-partitioned.neu
	meshPath := getTestMeshPath()

	// Create Element3D from partitioned mesh
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err, "Failed to create Element3D from partitioned mesh")
	require.NotNil(t, el)

	// The critical test: EToP must be populated
	require.NotNil(t, el.EToP, "Element3D.EToP MUST be populated from partitioned mesh file")
	assert.Equal(t, el.K, len(el.EToP), "EToP must have an entry for each element")

	// Verify the partitions are correct (p1=142, p2=141, p3=141, p4=141 elements)
	partitionCounts := make(map[int]int)
	for _, p := range el.EToP {
		partitionCounts[p]++
	}

	// Should have partitions 1,2,3,4
	assert.Equal(t, 4, len(partitionCounts))
	assert.Equal(t, 142, partitionCounts[1])
	assert.Equal(t, 141, partitionCounts[2])
	assert.Equal(t, 141, partitionCounts[3])
	assert.Equal(t, 141, partitionCounts[4])
}

// NEW TESTS FOR BUILDMAPS3D FIX

func TestBuildMaps3D_ValidVmapP(t *testing.T) {
	// Load the problematic mesh
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Check all VmapP values
	invalidCount := 0
	boundaryCount := 0
	interiorCount := 0

	for i := 0; i < len(el.VmapP); i++ {
		elemP := el.VmapP[i] / el.Np

		// Check if VmapP points to a valid element
		if elemP >= el.K {
			invalidCount++
			// Get details about this invalid mapping
			elem := i / (el.Nfp * 4)
			face := (i / el.Nfp) % 4
			point := i % el.Nfp

			t.Errorf("Invalid VmapP[%d] = %d points to element %d (max valid: %d) at elem=%d, face=%d, point=%d",
				i, el.VmapP[i], elemP, el.K-1, elem, face, point)
		}

		// Check if this is a boundary face (VmapP == VmapM)
		if el.VmapP[i] == el.VmapM[i] {
			boundaryCount++
		} else {
			interiorCount++
		}
	}

	// Report results
	t.Logf("Total face points: %d", len(el.VmapP))
	t.Logf("Boundary points: %d", boundaryCount)
	t.Logf("Interior points: %d", interiorCount)
	t.Logf("Invalid points: %d", invalidCount)

	if invalidCount > 0 {
		t.Fatalf("Found %d invalid VmapP entries", invalidCount)
	}

	// Additional validation: check boundary faces
	for k := 0; k < el.K; k++ {
		for f := 0; f < 4; f++ {
			// Check if this is a boundary face
			if el.EToE[k][f] == k && el.EToF[k][f] == f {
				// Verify all points on this face have VmapP = VmapM
				for p := 0; p < el.Nfp; p++ {
					idx := k*4*el.Nfp + f*el.Nfp + p
					if el.VmapP[idx] != el.VmapM[idx] {
						t.Errorf("Boundary face (elem=%d, face=%d) should have VmapP=VmapM at point %d", k, f, p)
					}
				}
			}
		}
	}
}

// TestBuildMaps3D_InteriorFaceConnectivity verifies interior face mappings
func TestBuildMaps3D_InteriorFaceConnectivity(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// For debugging, print first few elements' connectivity
	t.Logf("First few elements' connectivity:")
	for k := 0; k < min(3, el.K); k++ {
		t.Logf("Element %d: EToE=%v, EToF=%v", k, el.EToE[k], el.EToF[k])
	}

	// For each interior face, verify connectivity
	connectivityErrors := 0
	checkedFaces := 0

	for k1 := 0; k1 < el.K; k1++ {
		for f1 := 0; f1 < 4; f1++ {
			k2 := el.EToE[k1][f1]
			f2 := el.EToF[k1][f1]

			// Skip boundary faces
			if k2 == k1 && f2 == f1 {
				continue
			}

			// Skip invalid neighbors
			if k2 < 0 || k2 >= el.K {
				continue
			}

			checkedFaces++

			// For each point on the face
			hasError := false
			for p1 := 0; p1 < el.Nfp; p1++ {
				idx1 := k1*4*el.Nfp + f1*el.Nfp + p1

				// VmapP should point to a node in element k2
				vmapP1 := el.VmapP[idx1]
				elemP := vmapP1 / el.Np

				if elemP != k2 {
					if !hasError && connectivityErrors < 5 { // Limit debug output
						t.Errorf("Face (elem=%d,face=%d) -> (elem=%d,face=%d):", k1, f1, k2, f2)
						hasError = true
					}
					connectivityErrors++
					if connectivityErrors <= 20 {
						t.Errorf("  Point %d: VmapP=%d points to elem %d, expected %d",
							p1, vmapP1, elemP, k2)
					}
				}
			}
		}
	}

	t.Logf("Checked %d interior faces", checkedFaces)
	if connectivityErrors > 0 {
		t.Fatalf("Found %d connectivity errors", connectivityErrors)
	}
}

// TestBuildMaps3D_NodeMatching verifies that matching nodes have same coordinates
func TestBuildMaps3D_NodeMatching(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(1, meshPath) // Use order 1 for simpler debugging
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	NODETOL := 1e-7
	mismatchCount := 0

	// For each face point
	for i := 0; i < len(el.VmapM); i++ {
		// Skip boundary faces
		if el.VmapM[i] == el.VmapP[i] {
			continue
		}

		// Get coordinates of M and P nodes
		vmapM := el.VmapM[i]
		vmapP := el.VmapP[i]

		elemM := vmapM / el.Np
		nodeM := vmapM % el.Np
		elemP := vmapP / el.Np
		nodeP := vmapP % el.Np

		// Get coordinates
		xM := el.X.At(nodeM, elemM)
		yM := el.Y.At(nodeM, elemM)
		zM := el.Z.At(nodeM, elemM)

		xP := el.X.At(nodeP, elemP)
		yP := el.Y.At(nodeP, elemP)
		zP := el.Z.At(nodeP, elemP)

		// Check distance
		dx := xM - xP
		dy := yM - yP
		dz := zM - zP
		dist := math.Sqrt(dx*dx + dy*dy + dz*dz)

		if dist > NODETOL {
			mismatchCount++
			if mismatchCount <= 5 { // Limit output
				elem := i / (el.Nfp * 4)
				face := (i / el.Nfp) % 4
				point := i % el.Nfp

				t.Errorf("Node mismatch at face point %d (elem=%d,face=%d,point=%d):",
					i, elem, face, point)
				t.Errorf("  M: elem %d, node %d -> (%g, %g, %g)",
					elemM, nodeM, xM, yM, zM)
				t.Errorf("  P: elem %d, node %d -> (%g, %g, %g)",
					elemP, nodeP, xP, yP, zP)
				t.Errorf("  Distance: %g", dist)
			}
		}
	}

	if mismatchCount > 0 {
		t.Fatalf("Found %d node coordinate mismatches", mismatchCount)
	}
}

// DebugVmapP provides detailed debugging information
func DebugVmapP(el *Element3D) {
	fmt.Printf("=== VmapP Debug Information ===\n")
	fmt.Printf("Total elements (K): %d\n", el.K)
	fmt.Printf("Nodes per element (Np): %d\n", el.Np)
	fmt.Printf("Face points (Nfp): %d\n", el.Nfp)
	fmt.Printf("Total VmapP entries: %d\n", len(el.VmapP))

	// Count different types of mappings
	boundaryCount := 0
	interiorCount := 0
	invalidCount := 0

	// Track invalid entries for detailed reporting
	type invalidEntry struct {
		index, elem, face, point int
		vmapP, targetElem        int
	}
	var invalidEntries []invalidEntry

	for i := 0; i < len(el.VmapP); i++ {
		// Decode position
		elem := i / (el.Nfp * 4)
		face := (i / el.Nfp) % 4
		point := i % el.Nfp

		// Check where VmapP points
		targetElem := el.VmapP[i] / el.Np

		if targetElem >= el.K {
			invalidCount++
			if len(invalidEntries) < 10 { // Limit output
				invalidEntries = append(invalidEntries, invalidEntry{
					index: i, elem: elem, face: face, point: point,
					vmapP: el.VmapP[i], targetElem: targetElem,
				})
			}
		} else if el.VmapP[i] == el.VmapM[i] {
			boundaryCount++
		} else {
			interiorCount++
		}
	}

	fmt.Printf("\nMapping Statistics:\n")
	fmt.Printf("- Boundary face points: %d (%.1f%%)\n",
		boundaryCount, float64(boundaryCount)*100/float64(len(el.VmapP)))
	fmt.Printf("- Interior face points: %d (%.1f%%)\n",
		interiorCount, float64(interiorCount)*100/float64(len(el.VmapP)))
	fmt.Printf("- Invalid mappings: %d\n", invalidCount)

	if invalidCount > 0 {
		fmt.Printf("\n!!! FOUND %d INVALID VmapP ENTRIES !!!\n", invalidCount)
		fmt.Printf("First %d invalid entries:\n", len(invalidEntries))
		for _, inv := range invalidEntries {
			fmt.Printf("  VmapP[%d] = %d -> element %d (invalid!) at elem=%d, face=%d, point=%d\n",
				inv.index, inv.vmapP, inv.targetElem, inv.elem, inv.face, inv.point)
		}
	}

	// Check boundary face consistency
	fmt.Printf("\nBoundary Face Validation:\n")
	boundaryFaces := 0
	for k := 0; k < el.K && k < len(el.EToE); k++ {
		for f := 0; f < 4 && f < len(el.EToE[k]); f++ {
			if el.EToE[k][f] == k && el.EToF[k][f] == f {
				boundaryFaces++
				// Check all points on this boundary face
				allCorrect := true
				for p := 0; p < el.Nfp; p++ {
					idx := k*4*el.Nfp + f*el.Nfp + p
					if idx < len(el.VmapP) && el.VmapP[idx] != el.VmapM[idx] {
						allCorrect = false
					}
				}
				if !allCorrect {
					fmt.Printf("  WARNING: Boundary face (elem=%d, face=%d) has incorrect VmapP\n", k, f)
				}
			}
		}
	}
	fmt.Printf("Total boundary faces detected: %d\n", boundaryFaces)

	fmt.Printf("\n=== End Debug Information ===\n")
}
