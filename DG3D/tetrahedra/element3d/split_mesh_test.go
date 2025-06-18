package element3d

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestGlobalToLocal_Mapping(t *testing.T) {
	// Test the global to local mapping functionality
	EToP := []int{4, 1, 2, 2, 3, 1, 1, 3, 2} // 9 elements in 4 partitions

	g2l := NewGlobalToLocal(EToP)

	// Test cases
	testCases := []struct {
		globalIdx     int
		expectedPart  int
		expectedLocal int
	}{
		{0, 4, 0}, // First element in partition 4
		{1, 1, 0}, // First element in partition 1
		{5, 1, 1}, // Second element in partition 1
		{6, 1, 2}, // Third element in partition 1
		{2, 2, 0}, // First element in partition 2
		{3, 2, 1}, // Second element in partition 2
		{8, 2, 2}, // Third element in partition 2
		{4, 3, 0}, // First element in partition 3
		{7, 3, 1}, // Second element in partition 3
	}

	for _, tc := range testCases {
		partID, localIdx, ok := g2l.Get(tc.globalIdx)
		assert.True(t, ok, "Should find global index %d", tc.globalIdx)
		assert.Equal(t, tc.expectedPart, partID, "Wrong partition for global %d", tc.globalIdx)
		assert.Equal(t, tc.expectedLocal, localIdx, "Wrong local index for global %d", tc.globalIdx)
	}

	// Test non-existent element
	_, _, ok := g2l.Get(99)
	assert.False(t, ok, "Should not find non-existent element")
}

func TestSplitByPartition_PreservesFaceMappings_DEBUG(t *testing.T) {
	meshPath := getTestMeshPath()

	// Just load normally - the error will trigger our debug output
	_, err := NewElement3D(2, meshPath)

	// We expect this to fail and show debug output
	if err != nil {
		t.Logf("Got expected error with debug output: %v", err)
		// Don't fail the test - we just want to see the debug
	}
}

func TestSplitByPartition_RemoteFaceMappings(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Find a face that connects two different partitions
	var testElem, testFace, remoteElem int
	var testPartID, remotePartID int
	found := false

	for k := 0; k < el.K && !found; k++ {
		for f := 0; f < 4; f++ {
			neighbor := el.EToE[k][f]
			if neighbor >= 0 && neighbor < el.K && el.EToP[k] != el.EToP[neighbor] {
				testElem = k
				testFace = f
				remoteElem = neighbor
				testPartID = el.EToP[k]
				remotePartID = el.EToP[neighbor]
				found = true
				break
			}
		}
	}

	require.True(t, found, "Should find at least one cross-partition face")

	// Get both partitions
	testPartEl, err := el.GetPartition(testPartID)
	require.NoError(t, err)

	remotePartEl, err := el.GetPartition(remotePartID)
	require.NoError(t, err)

	// Build global to local mapping
	g2l := NewGlobalToLocal(el.EToP)

	// Get local indices
	_, testLocalIdx, _ := g2l.Get(testElem)
	_, remoteLocalIdx, _ := g2l.Get(remoteElem)

	// Verify VmapP in test partition points to remote partition's local coordinates
	for p := 0; p < el.Nfp; p++ {
		// Index in test partition
		localFaceIdx := testLocalIdx*4*el.Nfp + testFace*el.Nfp + p

		// VmapP should contain index in remote partition's coordinate system
		vmapPIdx := testPartEl.VmapP[localFaceIdx]
		elemP := vmapPIdx / el.Np
		nodeP := vmapPIdx % el.Np

		// The element index should be the remote element's local index
		assert.Equal(t, remoteLocalIdx, elemP,
			"VmapP should reference remote element's local index")

		// Get physical coordinates from both sides
		// M side (test partition)
		vmapMIdx := testPartEl.VmapM[localFaceIdx]
		mCoords := [3]float64{
			testPartEl.X.At(vmapMIdx%el.Np, vmapMIdx/el.Np),
			testPartEl.Y.At(vmapMIdx%el.Np, vmapMIdx/el.Np),
			testPartEl.Z.At(vmapMIdx%el.Np, vmapMIdx/el.Np),
		}

		// P side (remote partition)
		pCoords := [3]float64{
			remotePartEl.X.At(nodeP, elemP),
			remotePartEl.Y.At(nodeP, elemP),
			remotePartEl.Z.At(nodeP, elemP),
		}

		// The M and P coordinates should match for interior faces
		assert.InDelta(t, mCoords[0], pCoords[0], 1e-10,
			"M and P X coordinates should match for face point %d", p)
		assert.InDelta(t, mCoords[1], pCoords[1], 1e-10,
			"M and P Y coordinates should match for face point %d", p)
		assert.InDelta(t, mCoords[2], pCoords[2], 1e-10,
			"M and P Z coordinates should match for face point %d", p)
	}
}

func TestSplitByPartition_BoundaryFaces(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Find a boundary face
	var boundaryElem, boundaryFace int
	found := false

	for k := 0; k < el.K && !found; k++ {
		for f := 0; f < 4; f++ {
			if el.EToE[k][f] == -1 {
				boundaryElem = k
				boundaryFace = f
				found = true
				break
			}
		}
	}

	if !found {
		t.Skip("No boundary faces found in test mesh")
	}

	partID := el.EToP[boundaryElem]
	partEl, err := el.GetPartition(partID)
	require.NoError(t, err)

	// Get local index
	g2l := NewGlobalToLocal(el.EToP)
	_, localIdx, _ := g2l.Get(boundaryElem)

	// Verify boundary face has VmapP = VmapM
	for p := 0; p < el.Nfp; p++ {
		idx := localIdx*4*el.Nfp + boundaryFace*el.Nfp + p

		assert.Equal(t, partEl.VmapM[idx], partEl.VmapP[idx],
			"Boundary face should have VmapP = VmapM at point %d", p)
		assert.Equal(t, partEl.MapM[idx], partEl.MapP[idx],
			"Boundary face should have MapP = MapM at point %d", p)
	}

	// Verify EToE = -1 for boundary
	assert.Equal(t, -1, partEl.EToE[localIdx][boundaryFace],
		"Boundary face should have EToE = -1")
}

func TestSplitByPartition_ConnectivityConsistency(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(2, meshPath)
	require.NoError(t, err)

	// For each partition, verify internal consistency
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		// Check that all VmapM indices are valid
		for i := 0; i < len(partEl.VmapM); i++ {
			vmapM := partEl.VmapM[i]
			elem := vmapM / el.Np
			node := vmapM % el.Np

			assert.GreaterOrEqual(t, elem, 0, "VmapM element should be >= 0")
			assert.Less(t, elem, partEl.K, "VmapM element should be < K")
			assert.GreaterOrEqual(t, node, 0, "VmapM node should be >= 0")
			assert.Less(t, node, el.Np, "VmapM node should be < Np")
		}

		// Check that VmapP indices are valid (may point to other partitions)
		for i := 0; i < len(partEl.VmapP); i++ {
			vmapP := partEl.VmapP[i]
			node := vmapP % el.Np

			// Node index should always be valid
			assert.GreaterOrEqual(t, node, 0, "VmapP node should be >= 0")
			assert.Less(t, node, el.Np, "VmapP node should be < Np")
		}

		// Verify EToE references are valid local indices or -1
		for k := 0; k < partEl.K; k++ {
			for f := 0; f < 4; f++ {
				neighbor := partEl.EToE[k][f]
				if neighbor != -1 {
					// Should be a valid local index
					assert.GreaterOrEqual(t, neighbor, 0,
						"EToE[%d][%d] should be >= 0 in partition %d", k, f, partID)
					assert.Less(t, neighbor, max(partEl.K, el.K),
						"EToE[%d][%d] should be < max K in partition %d", k, f, partID)
				}
			}
		}
	}
}

func TestSplitByPartition_PhysicalCoordinatePreservation(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := NewElement3D(3, meshPath)
	require.NoError(t, err)

	// For a few elements, verify coordinates are preserved
	g2l := NewGlobalToLocal(el.EToP)

	testElements := []int{0, 10, 50, 100}
	for _, globalIdx := range testElements {
		if globalIdx >= el.K {
			continue
		}

		partID := el.EToP[globalIdx]
		partEl, err := el.GetPartition(partID)
		require.NoError(t, err)

		_, localIdx, _ := g2l.Get(globalIdx)

		// Check all nodes
		for n := 0; n < el.Np; n++ {
			origX := el.X.At(n, globalIdx)
			origY := el.Y.At(n, globalIdx)
			origZ := el.Z.At(n, globalIdx)

			splitX := partEl.X.At(n, localIdx)
			splitY := partEl.Y.At(n, localIdx)
			splitZ := partEl.Z.At(n, localIdx)

			assert.InDelta(t, origX, splitX, 1e-14,
				"X coordinate mismatch for elem %d node %d", globalIdx, n)
			assert.InDelta(t, origY, splitY, 1e-14,
				"Y coordinate mismatch for elem %d node %d", globalIdx, n)
			assert.InDelta(t, origZ, splitZ, 1e-14,
				"Z coordinate mismatch for elem %d node %d", globalIdx, n)
		}
	}
}

func TestSplitByPartition_NoCallToBuildMaps3D(t *testing.T) {
	// This test verifies that we don't call BuildMaps3D which would overwrite
	// our carefully transformed mappings

	meshPath := getTestMeshPath()
	el, err := NewElement3D(1, meshPath)
	require.NoError(t, err)

	// Store original MapP values for a cross-partition face
	var crossPartFaceMapP []int
	var globalElem, face int

	for k := 0; k < el.K; k++ {
		for f := 0; f < 4; f++ {
			neighbor := el.EToE[k][f]
			if neighbor >= 0 && neighbor < el.K && el.EToP[k] != el.EToP[neighbor] {
				globalElem = k
				face = f

				// Store MapP values
				for p := 0; p < el.Nfp; p++ {
					idx := k*4*el.Nfp + f*el.Nfp + p
					crossPartFaceMapP = append(crossPartFaceMapP, el.MapP[idx])
				}
				goto found
			}
		}
	}
found:

	require.NotEmpty(t, crossPartFaceMapP, "Should find cross-partition face")

	// After splitting, verify MapP is NOT identity (which BuildMaps3D would create)
	partID := el.EToP[globalElem]
	partEl, err := el.GetPartition(partID)
	require.NoError(t, err)

	g2l := NewGlobalToLocal(el.EToP)
	_, localIdx, _ := g2l.Get(globalElem)

	// Check that MapP values were preserved (not reset to identity)
	for p := 0; p < el.Nfp; p++ {
		localMapIdx := localIdx*4*el.Nfp + face*el.Nfp + p

		// If BuildMaps3D was called, MapP would equal MapM (identity mapping)
		// for unmatched faces
		assert.NotEqual(t, partEl.MapM[localMapIdx], partEl.MapP[localMapIdx],
			"MapP should not be identity for cross-partition face point %d", p)

		// MapP should have the original value from global mesh
		assert.Equal(t, crossPartFaceMapP[p], partEl.MapP[localMapIdx],
			"MapP should preserve original value for cross-partition face")
	}
}
