package facebuffer

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"testing"

	"github.com/notargets/gocfd/DG3D/tetrahedra"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// getTestMeshPath returns the path to the test mesh file
// It handles different working directory scenarios
func getTestMeshPath() string {
	// Try multiple possible paths
	possiblePaths := []string{
		"../mesh/cube-partitioned.neu",                    // From facebuffer directory
		"DG3D/mesh/cube-partitioned.neu",                  // From project root
		"mesh/cube-partitioned.neu",                       // From DG3D directory
		"../../DG3D/mesh/cube-partitioned.neu",            // From nested test directory
		filepath.Join("testdata", "cube-partitioned.neu"), // Local test data
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

func TestMeshLoading_Basic(t *testing.T) {
	// Test basic mesh loading first
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Logf("Failed to load mesh from: %s", meshPath)
		t.Logf("Error: %v", err)
		t.Logf("Working directory: %s", getWorkingDirectory())
		t.Skip("Cannot load test mesh file - skipping tests that require it")
	}

	require.NotNil(t, el)
	t.Logf("Successfully loaded mesh with %d elements", el.K)
	if el.EToP != nil {
		t.Logf("Mesh has partition data")
	} else {
		t.Logf("Mesh has no partition data")
	}
}

func getWorkingDirectory() string {
	wd, err := os.Getwd()
	if err != nil {
		return "unknown"
	}
	return wd
}

func TestNewFaceBufferBuilder_WithPartitionedMesh(t *testing.T) {
	// Load real partitioned mesh
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file %s: %v", meshPath, err)
	}
	require.NotNil(t, el)

	// Verify mesh has basic structure
	require.Greater(t, el.K, 0, "Mesh should have elements")

	// Skip if no partition data (mesh might not be partitioned)
	if el.EToP == nil {
		t.Skip("Test mesh is not partitioned, skipping partition-specific test")
	}

	builder := NewFaceBufferBuilder(el, 5) // 5 equations

	// Verify builder extracted correct dimensions
	assert.Equal(t, uint32(4), builder.Nface) // Tetrahedra have 4 faces
	assert.Equal(t, uint32(el.Nfp), builder.Nfp)
	assert.Equal(t, uint32(el.K), builder.K)
	assert.Equal(t, uint32(5), builder.Neq)

	// Verify partition info was extracted
	assert.Greater(t, builder.NPart, uint32(1), "Partitioned mesh should have multiple partitions")
	assert.Equal(t, uint32(el.EToP[0]), builder.MyPartID, "Should use first element's partition as local partition")

	expectedTotalFacePoints := uint32(4) * uint32(el.Nfp) * uint32(el.K)
	assert.Equal(t, expectedTotalFacePoints, builder.totalFacePoints)
}

func TestNewFaceBufferBuilder_NonPartitioned(t *testing.T) {
	// Load mesh and clear partition data to simulate non-partitioned
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file %s: %v", meshPath, err)
	}

	el.EToP = nil // Clear partition data

	builder := NewFaceBufferBuilder(el, 3)

	// Should default to single partition
	assert.Equal(t, uint32(1), builder.NPart)
	assert.Equal(t, uint32(0), builder.MyPartID)
}

// Add these updated tests to face_buffer_builder_test.go
// These tests properly use the split Element3D structures

func TestBuildFromElement3D_ProcessesRealConnectivity_WithSplitMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Mesh should be split")

	// Test each partition separately
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]
		t.Run(fmt.Sprintf("Partition_%d", partID), func(t *testing.T) {
			builder := NewFaceBufferBuilder(partEl, 3)

			err = builder.BuildFromElement3D(partEl)
			require.NoError(t, err, "Should successfully process partition %d connectivity", partID)

			stats := builder.GetBuildStatistics()

			// Verify all face points were processed
			expectedTotalPoints := uint32(partEl.K) * 4 * uint32(partEl.Nfp)
			assert.Equal(t, expectedTotalPoints, stats["total_face_points"])
			assert.Equal(t, stats["total_face_points"], stats["total_connections"])

			// Should have interior connections (faces between elements in same partition)
			assert.Greater(t, stats["interior_points"], uint32(0),
				"Partition %d should have interior faces", partID)

			// Should have boundary faces (domain boundaries + partition boundaries)
			assert.Greater(t, stats["boundary_points"], uint32(0),
				"Partition %d should have boundary faces", partID)

			// Should NOT have remote connections (partition boundaries marked as boundary)
			assert.Equal(t, uint32(0), stats["remote_points"],
				"Split mesh should not have remote points (partition boundaries are marked as boundary)")
		})
	}
}

func TestBuildFromElement3D_SinglePartitionProcessing(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(2, meshPath)
	require.NoError(t, err)

	// Test partition 4 specifically (141 elements)
	partEl, err := el.GetPartition(4)
	require.NoError(t, err)

	builder := NewFaceBufferBuilder(partEl, 3)

	// Verify builder uses correct partition info
	assert.Equal(t, uint32(4), builder.MyPartID)
	assert.Equal(t, uint32(141), builder.K)

	err = builder.BuildFromElement3D(partEl)
	require.NoError(t, err)

	stats := builder.GetBuildStatistics()

	// All connections should be classified
	totalClassified := stats["interior_points"] + stats["boundary_points"] + stats["remote_points"]
	assert.Equal(t, stats["total_connections"], totalClassified)

	// Log statistics for debugging
	t.Logf("Partition 4 statistics:")
	t.Logf("  Total face points: %d", stats["total_face_points"])
	t.Logf("  Interior points: %d", stats["interior_points"])
	t.Logf("  Boundary points: %d", stats["boundary_points"])
	t.Logf("  Remote points: %d", stats["remote_points"])
}

func TestBuild_CreatesValidRuntime_WithSplitMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Use partition 1 (142 elements)
	partEl, err := el.GetPartition(1)
	require.NoError(t, err)

	builder := NewFaceBufferBuilder(partEl, 3)
	runtime, err := builder.Build(partEl)
	require.NoError(t, err)
	require.NotNil(t, runtime)

	// Verify runtime dimensions match partition
	assert.Equal(t, uint32(partEl.K), runtime.K)
	assert.Equal(t, uint32(142), runtime.K) // Partition 1 has 142 elements
	assert.Equal(t, uint32(partEl.Nfp), runtime.Nfp)
	assert.Equal(t, builder.totalFacePoints, runtime.TotalFacePoints)

	// Verify arrays are properly sized
	assert.Equal(t, int(runtime.TotalFacePoints), len(runtime.FaceTypes))
	assert.Equal(t, int(runtime.TotalFacePoints), len(runtime.BCTypes))
	assert.Equal(t, int(runtime.TotalFacePoints), len(runtime.PartitionIDs))

	expectedMArraySize := int(runtime.TotalFacePoints * runtime.Neq)
	assert.Equal(t, expectedMArraySize, len(runtime.MArray))

	// Verify local indices count matches interior connections
	stats := builder.GetBuildStatistics()
	assert.Equal(t, int(stats["interior_points"]), len(runtime.LocalPIndices))
}

func TestPartitionBoundaryHandling_WithSplitMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Test that partition boundaries are handled correctly
	// In the split mesh, connections to other partitions should be marked as boundary
	for partIdx, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		builder := NewFaceBufferBuilder(partEl, 3)
		err = builder.BuildFromElement3D(partEl)
		require.NoError(t, err)

		// Count boundary faces that were originally interior in the global mesh
		partitionBoundaryCount := 0

		for k := 0; k < partEl.K; k++ {
			for f := 0; f < 4; f++ {
				if partEl.EToE[k][f] == -1 {
					// This is a boundary in the split mesh
					// Check if it was interior in the original mesh
					globalIdx, err := partEl.GetOriginalElementIndex(k)
					require.NoError(t, err)

					if globalIdx < el.K && f < len(el.EToE[globalIdx]) {
						globalNeighbor := el.EToE[globalIdx][f]
						if globalNeighbor != -1 {
							// Was interior in global mesh, now boundary
							partitionBoundaryCount++
						}
					}
				}
			}
		}

		t.Logf("Partition %d (idx %d): %d faces are partition boundaries",
			partID, partIdx, partitionBoundaryCount)

		// Each partition should have some partition boundary faces
		assert.Greater(t, partitionBoundaryCount, 0,
			"Partition %d should have partition boundary faces", partID)
	}
}

func TestValidateBuild_WithSplitMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Test validation for each partition
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		builder := NewFaceBufferBuilder(partEl, 3)
		err = builder.BuildFromElement3D(partEl)
		require.NoError(t, err)

		// Each partition should pass validation
		err = builder.ValidateBuild()
		assert.NoError(t, err, "Partition %d connectivity should validate successfully", partID)
	}
}

func TestNonPartitionedMesh_StillWorks(t *testing.T) {
	// Load the mesh but simulate non-partitioned by clearing partition data
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Clear partition data to simulate non-partitioned mesh
	el.EToP = nil
	el.SplitElement3D = nil

	// Should not have split meshes
	assert.Nil(t, el.SplitElement3D, "Non-partitioned mesh should not be split")
	assert.Nil(t, el.EToP, "Should not have partition data")

	// Face buffer should still work
	builder := NewFaceBufferBuilder(el, 3)

	// Verify it treats as single partition
	assert.Equal(t, uint32(1), builder.NPart)
	assert.Equal(t, uint32(0), builder.MyPartID)

	err = builder.BuildFromElement3D(el)
	require.NoError(t, err)

	stats := builder.GetBuildStatistics()

	// Should have both interior and boundary faces
	assert.Greater(t, stats["boundary_points"], uint32(0), "Should have boundary faces")
	assert.Greater(t, stats["interior_points"], uint32(0), "Should have interior faces")
	assert.Equal(t, uint32(0), stats["remote_points"], "Should have no remote points")

	// All face points should be classified
	totalClassified := stats["interior_points"] + stats["boundary_points"] + stats["remote_points"]
	assert.Equal(t, stats["total_connections"], totalClassified)
}

func TestConsistencyBetweenPartitions(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(2, meshPath)
	require.NoError(t, err)

	// Build face buffers for all partitions
	builders := make([]*FaceBufferBuilder, len(el.SplitElement3D))
	totalInterior := uint32(0)
	totalBoundary := uint32(0)
	totalFacePoints := uint32(0)

	for i, partEl := range el.SplitElement3D {
		builders[i] = NewFaceBufferBuilder(partEl, 3)
		err = builders[i].BuildFromElement3D(partEl)
		require.NoError(t, err)

		stats := builders[i].GetBuildStatistics()
		totalInterior += stats["interior_points"]
		totalBoundary += stats["boundary_points"]
		totalFacePoints += stats["total_face_points"]
	}

	// Total face points should match expected for entire mesh
	expectedTotal := uint32(el.K) * 4 * uint32(el.Nfp)
	assert.Equal(t, expectedTotal, totalFacePoints,
		"Sum of partition face points should equal total mesh face points")

	t.Logf("Mesh totals across all partitions:")
	t.Logf("  Total interior points: %d", totalInterior)
	t.Logf("  Total boundary points: %d", totalBoundary)
	t.Logf("  Total face points: %d", totalFacePoints)
}

// Helper function to count faces between partitions in the original mesh
func countInterPartitionFaces(el *tetrahedra.Element3D) map[string]int {
	// Key format: "partA-partB" where partA < partB
	interPartitionFaces := make(map[string]int)

	for k := 0; k < el.K; k++ {
		for f := 0; f < 4; f++ {
			neighbor := el.EToE[k][f]
			if neighbor != -1 && neighbor != k {
				partA := el.EToP[k]
				partB := el.EToP[neighbor]

				if partA != partB {
					// Create consistent key
					var key string
					if partA < partB {
						key = fmt.Sprintf("%d-%d", partA, partB)
					} else {
						key = fmt.Sprintf("%d-%d", partB, partA)
					}
					interPartitionFaces[key]++
				}
			}
		}
	}

	// Each face is counted twice (once from each side)
	for key := range interPartitionFaces {
		interPartitionFaces[key] /= 2
	}

	return interPartitionFaces
}

func TestInterPartitionFaceCount(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)

	if el.EToP == nil {
		t.Skip("Need partitioned mesh for this test")
	}

	// Count inter-partition faces in original mesh
	interPartitionFaces := countInterPartitionFaces(el)

	totalInterPartitionFaces := 0
	for key, count := range interPartitionFaces {
		t.Logf("Faces between partitions %s: %d", key, count)
		totalInterPartitionFaces += count
	}

	t.Logf("Total inter-partition faces: %d", totalInterPartitionFaces)

	// This information helps verify that partition boundaries are handled correctly
	assert.Greater(t, totalInterPartitionFaces, 0,
		"Partitioned mesh should have faces between partitions")
}
