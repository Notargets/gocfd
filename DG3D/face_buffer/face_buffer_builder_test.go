package facebuffer

import (
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

func TestBuildFromElement3D_ProcessesRealConnectivity(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file: %v", err)
	}

	builder := NewFaceBufferBuilder(el, 3)

	err = builder.BuildFromElement3D(el)
	require.NoError(t, err, "Should successfully process real mesh connectivity")

	stats := builder.GetBuildStatistics()

	// Verify all face points were processed
	expectedTotalPoints := uint32(el.K) * 4 * uint32(el.Nfp)
	assert.Equal(t, expectedTotalPoints, stats["total_face_points"])
	assert.Equal(t, stats["total_face_points"], stats["total_connections"])

	// Should have some interior connections (mesh has internal faces)
	assert.Greater(t, stats["interior_points"], uint32(0), "Real mesh should have interior faces")

	// Should have boundary faces
	assert.Greater(t, stats["boundary_points"], uint32(0), "Real mesh should have boundary faces")

	// Should have remote connections if partitioned
	if el.EToP != nil {
		assert.Greater(t, stats["remote_partitions"], uint32(0), "Partitioned mesh should have remote connections")
		assert.Greater(t, stats["remote_points"], uint32(0), "Should have remote face points")
	}
}

func TestBuildFromElement3D_HandlesPartitionBoundaries(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file: %v", err)
	}

	if el.EToP == nil {
		t.Skip("Need partitioned mesh for this test")
	}

	builder := NewFaceBufferBuilder(el, 4)
	err = builder.BuildFromElement3D(el)
	require.NoError(t, err)

	stats := builder.GetBuildStatistics()

	// Verify partition classification is working
	totalConnections := stats["interior_points"] + stats["boundary_points"] + stats["remote_points"]
	assert.Equal(t, stats["total_connections"], totalConnections, "All connections should be classified")

	// Check remote partition tracking
	assert.Equal(t, uint32(len(builder.remotePartitions)), stats["remote_partitions"])

	// Verify each remote partition has connections
	for partID, rpi := range builder.remotePartitions {
		assert.Greater(t, len(rpi.Connections), 0, "Remote partition %d should have connections", partID)
		assert.Equal(t, partID, rpi.PartitionID)
	}
}

func TestBuild_CreatesValidRuntime(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file: %v", err)
	}

	builder := NewFaceBufferBuilder(el, 3)
	runtime, err := builder.Build(el)
	require.NoError(t, err)
	require.NotNil(t, runtime)

	// Verify runtime dimensions match source
	assert.Equal(t, uint32(el.K), runtime.K)
	assert.Equal(t, uint32(el.Nfp), runtime.Nfp)
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

	// Verify remote buffers match remote partitions
	assert.Equal(t, len(builder.remotePartitions), len(runtime.RemoteBuffers))
}

func TestBuild_ValidatesConnectivity(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file: %v", err)
	}

	builder := NewFaceBufferBuilder(el, 3)

	// Should auto-build and validate successfully
	runtime, err := builder.Build(el)
	require.NoError(t, err)

	// Validate the runtime structure
	err = runtime.ValidateRuntime()
	assert.NoError(t, err, "Runtime structure should be valid")
}

func TestValidateBuild_WithRealMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file: %v", err)
	}

	builder := NewFaceBufferBuilder(el, 3)
	err = builder.BuildFromElement3D(el)
	require.NoError(t, err)

	// Real mesh should pass validation
	err = builder.ValidateBuild()
	assert.NoError(t, err, "Real mesh connectivity should validate successfully")
}

func TestPartitionClassification_WithRealMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file: %v", err)
	}

	if el.EToP == nil {
		t.Skip("Need partitioned mesh for this test")
	}

	builder := NewFaceBufferBuilder(el, 3)
	err = builder.BuildFromElement3D(el)
	require.NoError(t, err)

	// Verify connections are properly classified based on real EToE/EToP data
	for _, conn := range builder.connections {
		elemID := conn.MPos / (builder.Nface * builder.Nfp)
		faceID := (conn.MPos % (builder.Nface * builder.Nfp)) / builder.Nfp

		if int(elemID) < len(el.EToE) && int(faceID) < len(el.EToE[elemID]) {
			neighbor := el.EToE[elemID][faceID]

			if neighbor == -1 {
				// Should be boundary
				assert.Equal(t, BoundaryFace, conn.Type, "Boundary face should be classified as boundary")
			} else if el.EToP[neighbor] == el.EToP[elemID] {
				// Should be local neighbor
				assert.Equal(t, LocalNeighbor, conn.Type, "Local neighbor should be classified as local")
			} else {
				// Should be remote neighbor
				assert.Equal(t, RemoteNeighbor, conn.Type, "Remote neighbor should be classified as remote")
				assert.Equal(t, uint32(el.EToP[neighbor]), conn.PartitionID, "Remote partition ID should match EToP")
			}
		}
	}
}

func TestGetBuildStatistics_WithRealMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file: %v", err)
	}

	builder := NewFaceBufferBuilder(el, 4)
	err = builder.BuildFromElement3D(el)
	require.NoError(t, err)

	stats := builder.GetBuildStatistics()

	// Check all expected statistics are present
	expectedKeys := []string{
		"total_face_points", "interior_points", "total_connections",
		"remote_partitions", "boundary_points", "remote_points",
	}

	for _, key := range expectedKeys {
		_, exists := stats[key]
		assert.True(t, exists, "Missing statistic: %s", key)
		assert.GreaterOrEqual(t, stats[key], uint32(0), "Statistic %s should be non-negative", key)
	}

	// Verify accounting is correct
	assert.Equal(t, stats["total_face_points"], stats["total_connections"])
	totalClassified := stats["interior_points"] + stats["boundary_points"] + stats["remote_points"]
	assert.Equal(t, stats["total_connections"], totalClassified, "All connections should be classified")
}

func TestRemotePartitionTracking_WithRealMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	if err != nil {
		t.Skipf("Skipping test: failed to load mesh file: %v", err)
	}

	if el.EToP == nil {
		t.Skip("Need partitioned mesh for this test")
	}

	builder := NewFaceBufferBuilder(el, 3)
	err = builder.BuildFromElement3D(el)
	require.NoError(t, err)

	localPartID := builder.MyPartID

	// Find all unique remote partitions from EToP
	remotePartitions := make(map[int]bool)
	for elemID := 0; elemID < el.K; elemID++ {
		for faceID := 0; faceID < 4; faceID++ {
			if faceID < len(el.EToE[elemID]) {
				neighbor := el.EToE[elemID][faceID]
				if neighbor >= 0 && neighbor < len(el.EToP) {
					neighborPart := el.EToP[neighbor]
					if uint32(neighborPart) != localPartID {
						remotePartitions[neighborPart] = true
					}
				}
			}
		}
	}

	// Verify builder found all remote partitions
	assert.Equal(t, len(remotePartitions), len(builder.remotePartitions),
		"Should track all remote partitions found in mesh")

	for partID := range remotePartitions {
		_, exists := builder.remotePartitions[uint32(partID)]
		assert.True(t, exists, "Should track remote partition %d", partID)
	}
}
