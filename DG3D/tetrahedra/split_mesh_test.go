package tetrahedra

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestSplitByPartition_BasicFunctionality(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el)

	// Verify mesh is partitioned
	require.NotNil(t, el.EToP, "Test requires partitioned mesh")

	// Expected partition distribution from debug output
	expectedPartitions := map[int]int{
		1: 142,
		2: 141,
		3: 141,
		4: 141,
	}

	// Verify split was performed
	require.NotNil(t, el.SplitElement3D, "SplitElement3D should be populated")
	assert.Equal(t, 4, len(el.SplitElement3D), "Should have 4 partitions")

	// Check each partition
	actualPartitions := make(map[int]int)
	for _, partEl := range el.SplitElement3D {
		require.NotNil(t, partEl.EToP)
		require.Greater(t, len(partEl.EToP), 0)

		partID := partEl.EToP[0]
		actualPartitions[partID] = partEl.K

		// All elements in partition should have same partition ID
		for _, p := range partEl.EToP {
			assert.Equal(t, partID, p, "All elements in partition should have same ID")
		}
	}

	// Verify partition sizes match expected
	for partID, expectedK := range expectedPartitions {
		assert.Equal(t, expectedK, actualPartitions[partID],
			"Partition %d should have %d elements", partID, expectedK)
	}
}

func TestSplitByPartition_GeometricDataIntegrity(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)

	// Test each partition
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		// Check dimensions
		assert.Equal(t, el.Np, partEl.Np, "Np should match parent")
		assert.Equal(t, el.Nfp, partEl.Nfp, "Nfp should match parent")

		// Check vertex coordinates are properly sized
		assert.Equal(t, partEl.K*4, partEl.VX.Len(),
			"VX size incorrect for partition %d", partID)
		assert.Equal(t, partEl.K*4, partEl.VY.Len(),
			"VY size incorrect for partition %d", partID)
		assert.Equal(t, partEl.K*4, partEl.VZ.Len(),
			"VZ size incorrect for partition %d", partID)

		// Check physical coordinates matrices
		rows, cols := partEl.X.Dims()
		assert.Equal(t, el.Np, rows, "X rows should be Np")
		assert.Equal(t, partEl.K, cols, "X cols should be local K")

		// Check geometric factors if present
		if partEl.GeometricFactors != nil {
			rows, cols := partEl.J.Dims()
			assert.Equal(t, el.Np, rows, "J rows should be Np")
			assert.Equal(t, partEl.K, cols, "J cols should be local K")

			// Verify Jacobian is positive
			for k := 0; k < partEl.K; k++ {
				for n := 0; n < partEl.Np; n++ {
					J := partEl.J.At(n, k)
					assert.Greater(t, J, 0.0,
						"Jacobian must be positive at node %d, elem %d, partition %d", n, k, partID)
				}
			}
		}

		// Check face geometric factors if present
		if partEl.FaceGeometricFactors != nil {
			rows, cols := partEl.SJ.Dims()
			assert.Equal(t, el.Nfp*4, rows, "SJ rows should be Nfp*4")
			assert.Equal(t, partEl.K, cols, "SJ cols should be local K")
		}
	}
}

func TestSplitByPartition_ConnectivityRemapping(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)

	// For each partition, verify connectivity is properly remapped
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		require.NotNil(t, partEl.ConnectivityArrays,
			"Partition %d should have connectivity arrays", partID)

		// Check EToE and EToF dimensions
		assert.Equal(t, partEl.K, len(partEl.EToE),
			"EToE length should match K for partition %d", partID)
		assert.Equal(t, partEl.K, len(partEl.EToF),
			"EToF length should match K for partition %d", partID)

		// For each element in partition
		for k := 0; k < partEl.K; k++ {
			assert.Equal(t, 4, len(partEl.EToE[k]),
				"Each tet should have 4 faces")
			assert.Equal(t, 4, len(partEl.EToF[k]),
				"Each tet should have 4 faces")

			// Check that local neighbors are in valid range
			for f := 0; f < 4; f++ {
				neighbor := partEl.EToE[k][f]
				if neighbor != -1 {
					assert.GreaterOrEqual(t, neighbor, 0,
						"Local neighbor should be >= 0")
					assert.Less(t, neighbor, partEl.K,
						"Local neighbor should be < K")

					// Verify reciprocal connectivity for local neighbors
					nf := partEl.EToF[k][f]
					if nf >= 0 && nf < 4 {
						recipNeighbor := partEl.EToE[neighbor][nf]
						recipFace := partEl.EToF[neighbor][nf]

						// For local connections, should point back
						if recipNeighbor != -1 {
							assert.Equal(t, k, recipNeighbor,
								"Reciprocal neighbor should point back")
							assert.Equal(t, f, recipFace,
								"Reciprocal face should match")
						}
					}
				}
			}
		}
	}
}

func TestSplitByPartition_BoundaryConditions(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)

	// Skip if no boundary conditions in original mesh
	if el.BCType == nil {
		t.Skip("No boundary conditions in test mesh")
	}

	// For each partition, check BC data
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		if partEl.BCType != nil {
			// BC array should be sized for local elements
			expectedSize := partEl.K * 4 // 4 faces per tet
			assert.Equal(t, expectedSize, len(partEl.BCType),
				"BCType size incorrect for partition %d", partID)
		}
	}
}

func TestSplitByPartition_OriginalIndexMapping(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)

	// Build reverse mapping to verify
	globalToPartition := make(map[int]struct{ partIdx, localIdx int })

	for partIdx, partEl := range el.SplitElement3D {
		require.NotNil(t, partEl.tetIndices, "tetIndices should be populated")
		assert.Equal(t, partEl.K, len(partEl.tetIndices),
			"tetIndices length should match K")

		// Each local index should map to unique global index
		for localIdx := 0; localIdx < partEl.K; localIdx++ {
			globalIdx := partEl.tetIndices[localIdx]

			// Verify global index is in valid range
			assert.GreaterOrEqual(t, globalIdx, 0)
			assert.Less(t, globalIdx, el.K)

			// Check for duplicates
			if prev, exists := globalToPartition[globalIdx]; exists {
				t.Errorf("Global element %d mapped to multiple partitions: "+
					"partition %d (local %d) and partition %d (local %d)",
					globalIdx, partIdx, localIdx, prev.partIdx, prev.localIdx)
			}
			globalToPartition[globalIdx] = struct{ partIdx, localIdx int }{partIdx, localIdx}

			// Test GetOriginalElementIndex
			origIdx, err := partEl.GetOriginalElementIndex(localIdx)
			assert.NoError(t, err)
			assert.Equal(t, globalIdx, origIdx)
		}
	}

	// Verify all global elements are accounted for
	assert.Equal(t, el.K, len(globalToPartition),
		"All global elements should be mapped to partitions")
}

func TestSplitByPartition_GetPartition(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)

	// Test retrieving each partition
	for partID := 1; partID <= 4; partID++ {
		partEl, err := el.GetPartition(partID)
		require.NoError(t, err, "Should find partition %d", partID)
		require.NotNil(t, partEl)

		// Verify it's the correct partition
		assert.Equal(t, partID, partEl.EToP[0])
	}

	// Test non-existent partition
	_, err = el.GetPartition(99)
	assert.Error(t, err, "Should error for non-existent partition")
}

func TestSplitByPartition_NonPartitionedMesh(t *testing.T) {
	// Create a simple non-partitioned mesh
	meshPath := getTestMeshPath()
	el, err := NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Clear partition data to simulate non-partitioned mesh
	el.EToP = nil
	el.SplitElement3D = nil // Clear any existing split

	// Split should succeed but create no partitions
	err = el.SplitByPartition()
	assert.NoError(t, err)
	assert.Nil(t, el.SplitElement3D, "Non-partitioned mesh should have nil SplitElement3D")
}

func TestSplitByPartition_DataConsistency(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)

	// For a few elements, verify data matches between original and split
	testCases := []struct {
		globalElem int
		partition  int
	}{
		{0, 4},
		{142, 1}, // First element of partition 1
		{283, 2}, // Somewhere in partition 2
	}

	for _, tc := range testCases {
		if tc.globalElem >= el.K {
			continue // Skip if test element exceeds mesh size
		}

		partEl, err := el.GetPartition(tc.partition)
		require.NoError(t, err)

		// Find local index
		localIdx := -1
		for i, globalIdx := range partEl.tetIndices {
			if globalIdx == tc.globalElem {
				localIdx = i
				break
			}
		}

		if localIdx == -1 {
			// Element not in this partition, skip
			continue
		}

		// Compare vertex coordinates
		for v := 0; v < 4; v++ {
			origVX := el.VX.At(tc.globalElem*4 + v)
			splitVX := partEl.VX.At(localIdx*4 + v)
			assert.InDelta(t, origVX, splitVX, 1e-10,
				"VX mismatch for elem %d vertex %d", tc.globalElem, v)

			origVY := el.VY.At(tc.globalElem*4 + v)
			splitVY := partEl.VY.At(localIdx*4 + v)
			assert.InDelta(t, origVY, splitVY, 1e-10,
				"VY mismatch for elem %d vertex %d", tc.globalElem, v)

			origVZ := el.VZ.At(tc.globalElem*4 + v)
			splitVZ := partEl.VZ.At(localIdx*4 + v)
			assert.InDelta(t, origVZ, splitVZ, 1e-10,
				"VZ mismatch for elem %d vertex %d", tc.globalElem, v)
		}

		// Compare physical coordinates at a few nodes
		for n := 0; n < min(5, el.Np); n++ {
			origX := el.X.At(n, tc.globalElem)
			splitX := partEl.X.At(n, localIdx)
			assert.InDelta(t, origX, splitX, 1e-10,
				"X coordinate mismatch at node %d", n)

			origY := el.Y.At(n, tc.globalElem)
			splitY := partEl.Y.At(n, localIdx)
			assert.InDelta(t, origY, splitY, 1e-10,
				"Y coordinate mismatch at node %d", n)

			origZ := el.Z.At(n, tc.globalElem)
			splitZ := partEl.Z.At(n, localIdx)
			assert.InDelta(t, origZ, splitZ, 1e-10,
				"Z coordinate mismatch at node %d", n)
		}

		// Compare Jacobian if available
		if el.GeometricFactors != nil && partEl.GeometricFactors != nil {
			for n := 0; n < min(3, el.Np); n++ {
				origJ := el.J.At(n, tc.globalElem)
				splitJ := partEl.J.At(n, localIdx)
				assert.InDelta(t, origJ, splitJ, 1e-10,
					"Jacobian mismatch at node %d", n)
			}
		}
	}
}

func TestSplitByPartition_VolumeConservation(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(1, meshPath) // Order 1 for simpler volume calculation
	require.NoError(t, err)

	// Calculate total volume from original mesh
	originalVolume := 0.0
	for k := 0; k < el.K; k++ {
		// Get vertices
		v0 := [3]float64{el.VX.At(k*4 + 0), el.VY.At(k*4 + 0), el.VZ.At(k*4 + 0)}
		v1 := [3]float64{el.VX.At(k*4 + 1), el.VY.At(k*4 + 1), el.VZ.At(k*4 + 1)}
		v2 := [3]float64{el.VX.At(k*4 + 2), el.VY.At(k*4 + 2), el.VZ.At(k*4 + 2)}
		v3 := [3]float64{el.VX.At(k*4 + 3), el.VY.At(k*4 + 3), el.VZ.At(k*4 + 3)}

		vol := computeTetVolume(v0, v1, v2, v3)
		originalVolume += vol
	}

	// Calculate total volume from split meshes
	splitVolume := 0.0
	for _, partEl := range el.SplitElement3D {
		for k := 0; k < partEl.K; k++ {
			v0 := [3]float64{partEl.VX.At(k*4 + 0), partEl.VY.At(k*4 + 0), partEl.VZ.At(k*4 + 0)}
			v1 := [3]float64{partEl.VX.At(k*4 + 1), partEl.VY.At(k*4 + 1), partEl.VZ.At(k*4 + 1)}
			v2 := [3]float64{partEl.VX.At(k*4 + 2), partEl.VY.At(k*4 + 2), partEl.VZ.At(k*4 + 2)}
			v3 := [3]float64{partEl.VX.At(k*4 + 3), partEl.VY.At(k*4 + 3), partEl.VZ.At(k*4 + 3)}

			vol := computeTetVolume(v0, v1, v2, v3)
			splitVolume += vol
		}
	}

	// Volumes should match
	assert.InDelta(t, originalVolume, splitVolume, 1e-10*originalVolume,
		"Total volume should be conserved after splitting")
}
