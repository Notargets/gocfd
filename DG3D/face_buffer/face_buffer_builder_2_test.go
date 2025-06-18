package facebuffer

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/tetrahedra/element3d"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// TestLocalNeighborReciprocity verifies that M-P connections are reciprocal for local neighbors
func _TestLocalNeighborReciprocity(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := element3d.NewElement3D(2, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	// Test each partition
	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]
		t.Run(fmt.Sprintf("Partition_%d", partID), func(t *testing.T) {
			builder := NewFaceBufferBuilder(partEl, 3)
			err = builder.BuildFromElement3D(partEl)
			require.NoError(t, err)

			// Build a map from M position to connection for quick lookup
			mPosToConn := make(map[uint32]*BuildTimeConnection)
			for i := range builder.connections {
				conn := &builder.connections[i]
				mPosToConn[conn.MPos] = conn
			}

			// Track verified pairs to avoid double-checking
			verifiedPairs := make(map[string]bool)
			reciprocityFailures := 0

			// Check each local neighbor connection
			for _, conn := range builder.connections {
				if conn.Type != LocalNeighbor {
					continue
				}

				// Create unique key for this pair
				var pairKey string
				if conn.MPos < conn.PPos {
					pairKey = fmt.Sprintf("%d-%d", conn.MPos, conn.PPos)
				} else {
					pairKey = fmt.Sprintf("%d-%d", conn.PPos, conn.MPos)
				}

				// Skip if already verified
				if verifiedPairs[pairKey] {
					continue
				}
				verifiedPairs[pairKey] = true

				// Find the reciprocal connection
				recipConn, exists := mPosToConn[conn.PPos]
				if !exists {
					t.Errorf("No connection found at P position %d (reciprocal of M position %d)",
						conn.PPos, conn.MPos)
					reciprocityFailures++
					continue
				}

				// Verify reciprocity
				if recipConn.Type != LocalNeighbor {
					t.Errorf("Reciprocal connection at pos %d is not LocalNeighbor (is %v)",
						conn.PPos, recipConn.Type)
					reciprocityFailures++
					continue
				}

				if recipConn.PPos != conn.MPos {
					t.Errorf("Reciprocity failure: M[%d]->P[%d], but M[%d]->P[%d] (expected P[%d])",
						conn.MPos, conn.PPos, recipConn.MPos, recipConn.PPos, conn.MPos)
					reciprocityFailures++
				}
			}

			assert.Equal(t, 0, reciprocityFailures,
				"Found %d reciprocity failures in partition %d", reciprocityFailures, partID)
		})
	}
}

// TestLocalNeighborPhysicalCoordinates verifies M and P sides have matching physical coordinates
func _TestLocalNeighborPhysicalCoordinates(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := element3d.NewElement3D(1, meshPath) // Order 1 for simpler testing
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	const tolerance = 1e-10

	// Test first partition for detailed verification
	partEl := el.SplitElement3D[0]
	partID := partEl.EToP[0]

	builder := NewFaceBufferBuilder(partEl, 3)
	err = builder.BuildFromElement3D(partEl)
	require.NoError(t, err)

	// Build runtime to get the actual mapping
	runtime, err := builder.Build(partEl)
	require.NoError(t, err)

	// For each local neighbor connection, verify physical coordinates match
	coordinateMismatches := 0
	maxError := 0.0

	for i, faceType := range runtime.FaceTypes {
		if faceType != LocalNeighbor {
			continue
		}

		// Get M-side coordinates
		mElem, mFace, mPoint := unflattenFaceIndex(uint32(i), builder.Nface, builder.Nfp)
		mX, mY, mZ := getFacePointCoordinates(partEl, mElem, mFace, mPoint)

		// Get P-side coordinates
		pIdx := runtime.LocalPIndices[countLocalNeighborsBeforeIndex(runtime.FaceTypes, i)]
		pElem, pFace, pPoint := unflattenFaceIndex(pIdx, builder.Nface, builder.Nfp)
		pX, pY, pZ := getFacePointCoordinates(partEl, pElem, pFace, pPoint)

		// Check coordinate match
		dx := math.Abs(mX - pX)
		dy := math.Abs(mY - pY)
		dz := math.Abs(mZ - pZ)
		dist := math.Sqrt(dx*dx + dy*dy + dz*dz)

		if dist > tolerance {
			coordinateMismatches++
			if dist > maxError {
				maxError = dist
			}

			// Log first few mismatches for debugging
			if coordinateMismatches <= 5 {
				t.Logf("Coordinate mismatch at face point %d:", i)
				t.Logf("  M-side: elem %d, face %d, point %d -> (%.6e, %.6e, %.6e)",
					mElem, mFace, mPoint, mX, mY, mZ)
				t.Logf("  P-side: elem %d, face %d, point %d -> (%.6e, %.6e, %.6e)",
					pElem, pFace, pPoint, pX, pY, pZ)
				t.Logf("  Distance: %.6e", dist)
			}
		}
	}

	assert.Equal(t, 0, coordinateMismatches,
		"Found %d coordinate mismatches in partition %d (max error: %.6e)",
		coordinateMismatches, partID, maxError)
}

// TestFacePointOrdering verifies that face points are correctly matched despite orientation differences
func TestFacePointOrdering(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := element3d.NewElement3D(2, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	// Focus on one partition
	partEl := el.SplitElement3D[0]
	builder := NewFaceBufferBuilder(partEl, 3)
	err = builder.BuildFromElement3D(partEl)
	require.NoError(t, err)

	// For each element, verify face connectivity makes sense
	nfaces := uint32(4)
	elementsChecked := 0

	for k := uint32(0); k < builder.K; k++ {
		hasInteriorFace := false

		for f := uint32(0); f < nfaces; f++ {
			// Check all points on this face
			allPointsBoundary := true
			allPointsInterior := true

			for fp := uint32(0); fp < builder.Nfp; fp++ {
				mPos := k*nfaces*builder.Nfp + f*builder.Nfp + fp

				// Find this connection
				var conn *BuildTimeConnection
				for i := range builder.connections {
					if builder.connections[i].MPos == mPos {
						conn = &builder.connections[i]
						break
					}
				}

				require.NotNil(t, conn, "Should find connection for position %d", mPos)

				if conn.Type == BoundaryFace {
					allPointsInterior = false
				} else {
					allPointsBoundary = false
					hasInteriorFace = true
				}
			}

			// All points on a face should have the same type
			assert.True(t, allPointsBoundary || allPointsInterior,
				"Face %d of element %d has mixed boundary/interior points", f, k)
		}

		if hasInteriorFace {
			elementsChecked++
		}
	}

	t.Logf("Checked %d elements with interior faces", elementsChecked)
	assert.Greater(t, elementsChecked, 0, "Should have checked some elements with interior faces")
}

// TestBoundaryFaceClassification verifies boundary faces are correctly identified
func TestBoundaryFaceClassification(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := element3d.NewElement3D(1, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		builder := NewFaceBufferBuilder(partEl, 3)
		err = builder.BuildFromElement3D(partEl)
		require.NoError(t, err)

		// Check that boundary faces correspond to EToE = -1
		for k := 0; k < partEl.K; k++ {
			for f := 0; f < 4; f++ {
				isBoundaryInMesh := partEl.EToE[k][f] == -1

				// Check all points on this face
				for fp := 0; fp < int(builder.Nfp); fp++ {
					mPos := uint32(k)*builder.Nface*builder.Nfp +
						uint32(f)*builder.Nfp + uint32(fp)

					// Find connection
					var conn *BuildTimeConnection
					for i := range builder.connections {
						if builder.connections[i].MPos == mPos {
							conn = &builder.connections[i]
							break
						}
					}

					require.NotNil(t, conn)

					if isBoundaryInMesh {
						assert.Equal(t, BoundaryFace, conn.Type,
							"Face %d of element %d should be boundary", f, k)
					} else {
						assert.Equal(t, LocalNeighbor, conn.Type,
							"Face %d of element %d should be interior", f, k)
					}
				}
			}
		}

		t.Logf("Partition %d: Boundary face classification verified", partID)
	}
}

// TestIndexRangeValidity verifies all indices are within valid ranges
func TestIndexRangeValidity(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := element3d.NewElement3D(3, meshPath) // Higher order for more face points
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	partEl := el.SplitElement3D[0]
	builder := NewFaceBufferBuilder(partEl, 4) // 4 equations
	err = builder.BuildFromElement3D(partEl)
	require.NoError(t, err)

	runtime, err := builder.Build(partEl)
	require.NoError(t, err)

	// Verify all indices are in valid ranges
	maxMIndex := runtime.K * 4 * runtime.Nfp
	localPCount := 0

	for i, faceType := range runtime.FaceTypes {
		assert.Less(t, uint32(i), maxMIndex, "Face type index out of range")

		switch faceType {
		case LocalNeighbor:
			assert.Less(t, localPCount, len(runtime.LocalPIndices),
				"LocalPIndices access out of range")
			pIdx := runtime.LocalPIndices[localPCount]
			assert.Less(t, pIdx, maxMIndex,
				"P index %d out of range (max %d)", pIdx, maxMIndex)
			localPCount++

		case BoundaryFace:
			assert.Less(t, uint32(i), uint32(len(runtime.BCTypes)),
				"BCTypes access out of range")

		case RemoteNeighbor:
			assert.Less(t, uint32(i), uint32(len(runtime.PartitionIDs)),
				"PartitionIDs access out of range")
		}
	}

	assert.Equal(t, len(runtime.LocalPIndices), localPCount,
		"LocalPIndices count mismatch")
}

// Helper functions

// unflattenFaceIndex converts a flattened face point index back to (elem, face, point)
func unflattenFaceIndex(idx, nface, nfp uint32) (elem, face, point uint32) {
	elem = idx / (nface * nfp)
	remainder := idx % (nface * nfp)
	face = remainder / nfp
	point = remainder % nfp
	return
}

// getFacePointCoordinates gets the physical coordinates of a face point
func getFacePointCoordinates(el *element3d.Element3D, elem, face, point uint32) (x, y, z float64) {
	// Get the volume node index for this face point
	faceNodeIdx := el.Fmask[face][point]

	// Get physical coordinates
	x = el.X.At(int(faceNodeIdx), int(elem))
	y = el.Y.At(int(faceNodeIdx), int(elem))
	z = el.Z.At(int(faceNodeIdx), int(elem))
	return
}

// countLocalNeighborsBeforeIndex counts how many LocalNeighbor entries appear before index i
func countLocalNeighborsBeforeIndex(faceTypes []FaceType, idx int) int {
	count := 0
	for i := 0; i < idx; i++ {
		if faceTypes[i] == LocalNeighbor {
			count++
		}
	}
	return count
}

// TestNaturalTraversalOrder verifies the face buffer follows the expected traversal order
func TestNaturalTraversalOrder(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := element3d.NewElement3D(1, meshPath)
	require.NoError(t, err)

	partEl := el.SplitElement3D[0]
	builder := NewFaceBufferBuilder(partEl, 3)
	err = builder.BuildFromElement3D(partEl)
	require.NoError(t, err)

	// Verify connections are in natural traversal order
	for i, conn := range builder.connections {
		expectedMPos := uint32(i)
		assert.Equal(t, expectedMPos, conn.MPos,
			"Connection %d has wrong M position", i)

		// Verify the position matches element/face/point traversal
		elem, face, point := unflattenFaceIndex(conn.MPos, builder.Nface, builder.Nfp)
		reconstructedPos := elem*builder.Nface*builder.Nfp + face*builder.Nfp + point
		assert.Equal(t, conn.MPos, reconstructedPos,
			"Position reconstruction failed for connection %d", i)
	}
}

// TestPartitionBoundaryAsRemote would test remote neighbors, but in the split mesh
// architecture, partition boundaries are marked as domain boundaries, not remote.
// This test is included for completeness but expects no remote connections.
func TestPartitionBoundaryAsRemote(t *testing.T) {
	meshPath := getTestMeshPath()
	el, err := element3d.NewElement3D(1, meshPath)
	require.NoError(t, err)
	require.NotNil(t, el.SplitElement3D, "Need split mesh")

	for _, partEl := range el.SplitElement3D {
		partID := partEl.EToP[0]

		builder := NewFaceBufferBuilder(partEl, 3)
		err = builder.BuildFromElement3D(partEl)
		require.NoError(t, err)

		// Count remote connections
		remoteCount := 0
		for _, conn := range builder.connections {
			if conn.Type == RemoteNeighbor {
				remoteCount++
			}
		}

		assert.Equal(t, 0, remoteCount,
			"Partition %d should have no remote neighbors in split mesh (found %d)",
			partID, remoteCount)
	}
}
