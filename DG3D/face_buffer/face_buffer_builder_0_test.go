package facebuffer

// Level 0 Tests: The simplest possible face buffer builder tests
//
// These tests form the foundation of a hierarchical testing approach:
// - Level 0: Single element and simple non-partitioned meshes (this file)
// - Level 1: More complex non-partitioned meshes
// - Level 2: Simple partitioned meshes
// - Level 3: Complex partitioned meshes
// - Level 4: Real mesh files with full features
//
// The goal is to isolate and test basic functionality before adding complexity.

import (
	"github.com/notargets/gocfd/DG3D/tetrahedra/tetelement"
	"testing"

	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// Level 0: Single element tests - the absolute simplest cases
// These tests use synthetic meshes with known, verifiable properties

func _TestLevel0_SingleTetAllBoundary(t *testing.T) {
	// Create a single tetrahedron mesh using test helpers
	tm := mesh.GetStandardTestMeshes()
	singleTetMesh := &mesh.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []mesh.ElementSet{tm.SingleTet},
		Dimension: 3,
		BoundingBox: [2][3]float64{
			{0, 0, 0},
			{1, 1, 1},
		},
	}

	// Convert to mesh and create Element3D
	meshObj := singleTetMesh.ConvertToMesh()

	// Create Element3D from mesh - order 1 for simplicity
	el, err := createElementFromMesh(meshObj, 1)
	require.NoError(t, err, "Failed to create Element3D")
	require.NotNil(t, el)

	// Basic validations
	assert.Equal(t, 1, el.K, "Should have exactly 1 element")
	assert.Nil(t, el.EToP, "Should not have partition data")
	assert.Nil(t, el.SplitElement3D, "Should not be split")

	// Verify connectivity arrays
	require.NotNil(t, el.ConnectivityArrays, "Should have connectivity arrays")
	require.NotNil(t, el.EToE, "Should have EToE")
	require.NotNil(t, el.EToF, "Should have EToF")
	assert.Equal(t, 1, len(el.EToE), "Should have 1 element in EToE")
	assert.Equal(t, 4, len(el.EToE[0]), "Element should have 4 faces in EToE")

	// For a single tet, all neighbors should be -1 (boundary)
	for f := 0; f < 4; f++ {
		assert.Equal(t, -1, el.EToE[0][f], "Face %d should be boundary", f)
		assert.Equal(t, -1, el.EToF[0][f], "Face %d should have no neighbor", f)
	}

	// Create face buffer builder
	builder := NewFaceBufferBuilder(el, 3) // 3 equations

	// Verify builder dimensions
	assert.Equal(t, uint32(4), builder.Nface, "Should have 4 faces per tet")
	assert.Equal(t, uint32(1), builder.K, "Should have 1 element")
	assert.Equal(t, uint32(1), builder.NPart, "Non-partitioned mesh has 1 partition")
	assert.Equal(t, uint32(0), builder.MyPartID, "Default partition ID is 0")

	// Build face buffer
	err = builder.BuildFromElement3D(el)
	require.NoError(t, err, "Should build successfully")

	// Get statistics
	stats := builder.GetBuildStatistics()

	// For a single tet, all faces should be boundary
	expectedBoundaryPoints := uint32(4) * builder.Nfp // 4 faces * points per face
	assert.Equal(t, expectedBoundaryPoints, stats["boundary_points"],
		"All face points should be boundary")
	assert.Equal(t, uint32(0), stats["interior_points"],
		"Should have no interior points")
	assert.Equal(t, uint32(0), stats["remote_points"],
		"Should have no remote points")

	// Verify all connections are boundary type
	for i, conn := range builder.connections {
		assert.Equal(t, BoundaryFace, conn.Type,
			"Connection %d should be boundary", i)
		assert.Equal(t, uint32(i), conn.MPos,
			"Connection %d should have correct M position", i)
		assert.Equal(t, uint32(0), conn.PPos,
			"Boundary faces should have PPos = 0")
	}

	// Validate the build
	err = builder.ValidateBuild()
	assert.NoError(t, err, "Single tet should validate successfully")
}

func _TestLevel0_TwoTetsSharedFace(t *testing.T) {
	// Create two tetrahedra sharing one face
	tm := mesh.GetStandardTestMeshes()
	twoTetMesh := tm.TwoTetMesh.ConvertToMesh()

	// Create Element3D
	el, err := createElementFromMesh(twoTetMesh, 1)
	require.NoError(t, err)
	require.NotNil(t, el)

	// Basic validations
	assert.Equal(t, 2, el.K, "Should have 2 elements")
	assert.Equal(t, 5, el.Mesh.NumVertices, "Should have 5 vertices")

	// Verify connectivity
	require.NotNil(t, el.EToE, "Should have EToE")
	require.NotNil(t, el.EToF, "Should have EToF")
	assert.Equal(t, 2, len(el.EToE), "Should have 2 elements in connectivity")

	// Find the shared face - exactly one face per tet should connect to the other
	sharedFaceCount := 0
	for elem := 0; elem < 2; elem++ {
		for face := 0; face < 4; face++ {
			if el.EToE[elem][face] == (1 - elem) { // Points to the other element
				sharedFaceCount++
			}
		}
	}
	assert.Equal(t, 2, sharedFaceCount, "Should have exactly 2 shared face connections (one from each side)")

	// Create and build face buffer
	builder := NewFaceBufferBuilder(el, 3)
	err = builder.BuildFromElement3D(el)
	require.NoError(t, err)

	stats := builder.GetBuildStatistics()

	// Two tets sharing one face should have:
	// - 6 boundary faces (3 per tet * 2 tets - 1 shared on each side)
	// - 2 interior faces (the shared face counted from each side)
	expectedBoundaryPoints := uint32(6) * builder.Nfp
	expectedInteriorPoints := uint32(2) * builder.Nfp

	assert.Equal(t, expectedBoundaryPoints, stats["boundary_points"],
		"Should have correct number of boundary points")
	assert.Equal(t, expectedInteriorPoints, stats["interior_points"],
		"Should have correct number of interior points")
	assert.Equal(t, uint32(0), stats["remote_points"],
		"Should have no remote points in non-partitioned mesh")

	// Find and verify the shared face connections
	sharedFaceCount = 0
	for _, conn := range builder.connections {
		if conn.Type == LocalNeighbor {
			sharedFaceCount++

			// Verify reciprocity: the neighbor's neighbor should point back
			neighborMPos := conn.PPos
			var recipConn *BuildTimeConnection
			for j := range builder.connections {
				if builder.connections[j].MPos == neighborMPos {
					recipConn = &builder.connections[j]
					break
				}
			}

			require.NotNil(t, recipConn,
				"Should find reciprocal connection for shared face")
			assert.Equal(t, LocalNeighbor, recipConn.Type,
				"Reciprocal connection should also be LocalNeighbor")
			assert.Equal(t, conn.MPos, recipConn.PPos,
				"Reciprocal connection should point back")
		}
	}

	assert.Equal(t, int(expectedInteriorPoints), sharedFaceCount,
		"Should have correct number of shared face connections")

	// Validate the build
	err = builder.ValidateBuild()
	assert.NoError(t, err, "Two tet mesh should validate successfully")
}

func _TestLevel0_CubeMeshNonPartitioned(t *testing.T) {
	// Use the standard cube mesh (6 tets)
	tm := mesh.GetStandardTestMeshes()
	cubeMesh := tm.CubeMesh.ConvertToMesh()

	// Create Element3D
	el, err := createElementFromMesh(cubeMesh, 1)
	require.NoError(t, err)
	require.NotNil(t, el)

	assert.Equal(t, 6, el.K, "Cube should have 6 tetrahedra")
	assert.Equal(t, 9, el.Mesh.NumVertices, "Cube should have 9 vertices") // 8 corners + 1 center

	// Verify connectivity exists
	require.NotNil(t, el.EToE, "Should have EToE connectivity")
	require.NotNil(t, el.EToF, "Should have EToF connectivity")
	assert.Equal(t, 6, len(el.EToE), "Should have connectivity for 6 elements")

	// Build face buffer
	builder := NewFaceBufferBuilder(el, 3)
	err = builder.BuildFromElement3D(el)
	require.NoError(t, err)

	stats := builder.GetBuildStatistics()

	// Cube mesh should have:
	// - 12 exterior faces (2 per cube face * 6 cube faces)
	// - Interior faces connecting the tets
	assert.Greater(t, stats["boundary_points"], uint32(0),
		"Should have boundary faces on cube exterior")
	assert.Greater(t, stats["interior_points"], uint32(0),
		"Should have interior faces between tets")
	assert.Equal(t, uint32(0), stats["remote_points"],
		"Should have no remote points")

	// Verify total face points are classified
	totalPoints := builder.K * builder.Nface * builder.Nfp
	totalClassified := stats["boundary_points"] + stats["interior_points"] + stats["remote_points"]
	assert.Equal(t, totalPoints, totalClassified,
		"All face points should be classified")

	// Validate
	err = builder.ValidateBuild()
	assert.NoError(t, err, "Cube mesh should validate successfully")
}

// Helper function to create Element3D from a mesh
func createElementFromMesh(meshObj *mesh.Mesh,
	order int) (*tetelement.Element3D, error) {
	// Use the actual Element3D constructor
	return tetelement.NewElement3DFromMesh(order, meshObj)
}
