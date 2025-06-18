package facebuffer

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/DG3D/tetrahedra"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestPartitionBoundaryBCTypes verifies partition boundaries are correctly marked with BCType=255
func TestPartitionBoundaryBCTypes(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	// Check each partition for partition boundary BC types
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		partitionBoundaryCount := 0
		domainBoundaryCount := 0

		for k := 0; k < partEl.K; k++ {
			for f := 0; f < 4; f++ {
				if partEl.EToE[k][f] == -1 {
					bcIdx := k*4 + f
					bcType := partEl.BCType[bcIdx]

					if bcType == tetrahedra.BCPartitionBoundary {
						partitionBoundaryCount++
					} else {
						domainBoundaryCount++
					}
				}
			}
		}

		t.Logf("Partition %d: %d partition boundaries, %d domain boundaries",
			partID, partitionBoundaryCount, domainBoundaryCount)

		// Each partition should have some partition boundaries (except possibly corner cases)
		// Most partitions should have partition boundaries
		assert.Greater(t, partitionBoundaryCount, 0,
			"Partition %d should have partition boundary faces marked with BCType=255", partID)
	}
}

// TestBuildStatisticsWithPartitionBoundaries verifies statistics distinguish partition vs domain boundaries
func TestBuildStatisticsWithPartitionBoundaries(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]
		builder := NewFaceBufferBuilder(partEl, 3)
		err = builder.BuildFromElement3D(partEl)
		require.NoError(t, err)

		stats := builder.GetBuildStatistics()

		t.Logf("Partition %d statistics:", partID)
		t.Logf("  Total boundary points: %d", stats["boundary_points"])
		t.Logf("  Domain boundary points: %d", stats["domain_boundary_points"])
		t.Logf("  Partition boundary points: %d", stats["partition_boundary_points"])

		// Verify that partition boundaries are being counted
		assert.Greater(t, stats["partition_boundary_points"], uint32(0),
			"Partition %d should have partition boundary points", partID)

		// Total boundary points should equal domain + partition boundaries
		totalBoundaries := stats["domain_boundary_points"] + stats["partition_boundary_points"]
		assert.Equal(t, stats["boundary_points"], totalBoundaries,
			"Total boundary points should equal domain + partition boundaries")
	}
}

// TestPartitionBoundaryProcessing verifies partition boundaries are processed correctly
func TestPartitionBoundaryProcessing(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(2, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	// Test partition 4 specifically
	partEl, err := el.GetPartition(4)
	require.NoError(t, err)

	builder := NewFaceBufferBuilder(partEl, 3)
	err = builder.BuildFromElement3D(partEl)
	require.NoError(t, err)

	runtime, err := builder.Build(partEl)
	require.NoError(t, err)

	// Count different BC types
	bcTypeCounts := make(map[uint32]int)
	for i, faceType := range runtime.FaceTypes {
		if faceType == BoundaryFace {
			bcType := runtime.BCTypes[i]
			bcTypeCounts[bcType]++
		}
	}

	t.Logf("Partition 4 BC type distribution:")
	for bcType, count := range bcTypeCounts {
		bcName := tetrahedra.GetBCTypeName(int(bcType))
		t.Logf("  BC Type %d (%s): %d faces", bcType, bcName, count)
	}

	// Should have partition boundaries
	assert.Greater(t, bcTypeCounts[uint32(tetrahedra.BCPartitionBoundary)], 0,
		"Partition 4 should have partition boundary faces (BCType=255)")
}

// TestNoPartitionBoundariesInGlobalMesh verifies global mesh doesn't have partition boundaries
func TestNoPartitionBoundariesInGlobalMesh(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Check the global mesh before splitting
	partitionBoundaryCount := 0
	for i := 0; i < len(el.BCType); i++ {
		if el.BCType[i] == tetrahedra.BCPartitionBoundary {
			partitionBoundaryCount++
		}
	}

	assert.Equal(t, 0, partitionBoundaryCount,
		"Global mesh should not have partition boundaries before splitting")
}

// TestPartitionBoundarySymmetry verifies partition boundaries are symmetric between partitions
func TestPartitionBoundarySymmetry(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	// Count inter-partition faces in the original mesh
	interPartitionFaces := make(map[string]int)
	for k := 0; k < el.K; k++ {
		for f := 0; f < 4; f++ {
			neighbor := el.EToE[k][f]
			if neighbor != -1 && neighbor != k {
				partA := el.EToP[k]
				partB := el.EToP[neighbor]

				if partA != partB {
					key := fmt.Sprintf("%d-%d", min(partA, partB), max(partA, partB))
					interPartitionFaces[key]++
				}
			}
		}
	}

	// Each face is counted twice (once from each side)
	for key := range interPartitionFaces {
		interPartitionFaces[key] /= 2
	}

	t.Logf("Inter-partition face counts:")
	totalInterPartitionFaces := 0
	for key, count := range interPartitionFaces {
		t.Logf("  %s: %d faces", key, count)
		totalInterPartitionFaces += count
	}

	// Now count partition boundaries in split meshes
	totalPartitionBoundaries := 0
	for _, partEl := range el.SplitElement3D {
		for k := 0; k < partEl.K; k++ {
			for f := 0; f < 4; f++ {
				if partEl.EToE[k][f] == -1 {
					bcIdx := k*4 + f
					if partEl.BCType[bcIdx] == tetrahedra.BCPartitionBoundary {
						totalPartitionBoundaries++
					}
				}
			}
		}
	}

	t.Logf("Total inter-partition faces in original: %d", totalInterPartitionFaces)
	t.Logf("Total partition boundaries in split meshes: %d", totalPartitionBoundaries)

	// Each inter-partition face should appear as a partition boundary on both sides
	assert.Equal(t, totalInterPartitionFaces*2, totalPartitionBoundaries,
		"Each inter-partition face should create two partition boundaries (one on each side)")
}

// TestPartition4HasBoundaryFaces verifies the specific failing test case
func TestPartition4HasBoundaryFaces(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := tetrahedra.NewElement3D(1, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Mesh should be split")

	// Test partition 4 specifically which was failing
	var part4Found bool
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]
		if partID == 4 {
			part4Found = true

			builder := NewFaceBufferBuilder(partEl, 3)
			err = builder.BuildFromElement3D(partEl)
			require.NoError(t, err, "Should successfully process partition 4 connectivity")

			stats := builder.GetBuildStatistics()

			// The original failing assertion
			assert.Greater(t, stats["boundary_points"], uint32(0),
				"Partition 4 should have boundary faces")

			// Additional diagnostics
			t.Logf("Partition 4 detailed statistics:")
			t.Logf("  Elements: %d", partEl.K)
			t.Logf("  Total face points: %d", stats["total_face_points"])
			t.Logf("  Interior points: %d", stats["interior_points"])
			t.Logf("  Total boundary points: %d", stats["boundary_points"])
			t.Logf("  Domain boundary points: %d", stats["domain_boundary_points"])
			t.Logf("  Partition boundary points: %d", stats["partition_boundary_points"])

			// Count faces with EToE = -1
			boundaryFaceCount := 0
			for k := 0; k < partEl.K; k++ {
				for f := 0; f < 4; f++ {
					if partEl.EToE[k][f] == -1 {
						boundaryFaceCount++
					}
				}
			}
			t.Logf("  Faces with EToE = -1: %d", boundaryFaceCount)
		}
	}

	require.True(t, part4Found, "Partition 4 should exist in split mesh")
}
