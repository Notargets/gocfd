package tetrahedra

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"testing"

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

func TestElement3D_ValidateMeshStructure(t *testing.T) {
	meshPath := getTestMeshPath()

	el, err := NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Log mesh statistics
	t.Logf("Mesh statistics:")
	t.Logf("  Elements (K): %d", el.K)
	t.Logf("  Vertices: %d", el.Mesh.NumVertices)
	t.Logf("  Element type: %v", el.Mesh.ElementTypes[0])

	// Count boundary vs interior faces
	boundaryFaces := 0
	interiorFaces := 0

	for k := 0; k < el.K; k++ {
		for f := 0; f < 4; f++ {
			if el.EToE[k][f] == -1 {
				boundaryFaces++
			} else if el.EToE[k][f] != k {
				interiorFaces++
			}
		}
	}

	t.Logf("  Boundary faces: %d", boundaryFaces)
	t.Logf("  Interior faces: %d", interiorFaces)

	// Sanity check: interior faces should be even (each counted from both sides)
	assert.Equal(t, 0, interiorFaces%2, "Interior faces should be even (counted from both sides)")
}
