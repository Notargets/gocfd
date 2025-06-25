// Package facebuffer implements efficient face data indexing for DG methods
package facebuffer

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/tetrahedra/tetelement"
)

// Face index encoding constants
const (
	BoundaryPlaceholder int32 = -999
	RemoteFace          int32 = -9999
)

// FaceBuffer holds the runtime face indexing data
type FaceBuffer struct {
	// Dimensions
	Nfp    uint32 // Face points per face
	Nfaces uint32 // Faces per element
	K      uint32 // Number of elements in this partition

	// Face connectivity (face-level indexing)
	FaceIndex []int32 // Size: Nfaces × K, encodes all face types

	// Remote partition send indices
	RemoteSendIndices map[uint32][]uint32 // [partitionID] -> indices into M buffer
}

// BuildFaceBuffer creates face buffer from Element3D connectivity
// Uses EToE and EToP to determine face connections:
// - EToE[elem][face] == elem: Boundary face
// - EToP[EToE[elem][face]] != EToP[elem]: Remote face
// - Otherwise: Interior face
func BuildFaceBuffer(el *tetelement.Element3D) (*FaceBuffer, error) {
	// Extract dimensions
	nfp := uint32(el.Nfp)
	nfaces := uint32(4) // Tetrahedra have 4 faces
	k := uint32(el.K)

	// Initialize result
	fb := &FaceBuffer{
		Nfp:               nfp,
		Nfaces:            nfaces,
		K:                 k,
		FaceIndex:         make([]int32, nfaces*k),
		RemoteSendIndices: make(map[uint32][]uint32),
	}

	// Determine if this is a parallel run
	isParallel := el.EToP != nil
	var myPartID int
	if isParallel && len(el.EToP) > 0 {
		// Find this partition's ID from first element
		myPartID = el.EToP[0]
	}

	// Process each face
	for elem := uint32(0); elem < k; elem++ {
		for face := uint32(0); face < nfaces; face++ {
			faceIdx := face + elem*nfaces
			neighborElem := el.EToE[elem][face]

			// Check if boundary (self-reference)
			if neighborElem == int(elem) {
				fb.FaceIndex[faceIdx] = BoundaryPlaceholder
				continue
			}

			// Check if remote (different partition)
			if isParallel && el.EToP[neighborElem] != myPartID {
				fb.FaceIndex[faceIdx] = RemoteFace

				// Track which M values need to be sent to neighbor partition
				neighborPartID := uint32(el.EToP[neighborElem])
				// Add all face points for this face to send list
				for p := uint32(0); p < nfp; p++ {
					mIdx := elem*nfaces*nfp + face*nfp + p
					fb.RemoteSendIndices[neighborPartID] = append(
						fb.RemoteSendIndices[neighborPartID], mIdx)
				}
				continue
			}

			// Interior face - find P position in neighbor's M buffer
			pFaceStart := findNeighborFaceStart(el, elem, face, uint32(neighborElem))
			fb.FaceIndex[faceIdx] = int32(pFaceStart)
		}
	}

	return fb, nil
}

// findNeighborFaceStart finds the start position of the matching face in neighbor's M buffer
func findNeighborFaceStart(el *tetelement.Element3D, elem, face, neighborElem uint32) uint32 {
	nfp := uint32(el.Nfp)
	nfaces := uint32(4)

	// Get a reference face point to match
	refMIdx := elem*nfaces*nfp + face*nfp // First point of this face
	refVolNode := el.VmapP[refMIdx]       // Where it connects to

	// Search neighbor's faces for matching connection
	for nFace := uint32(0); nFace < nfaces; nFace++ {
		// Check first point of neighbor's face
		nIdx := neighborElem*nfaces*nfp + nFace*nfp
		if el.VmapM[nIdx] == refVolNode {
			// Found the matching face - return its start position
			return nIdx
		}
	}

	// This should never happen with valid connectivity
	panic(fmt.Sprintf("No matching face found in neighbor element %d for element %d face %d",
		neighborElem, elem, face))
}

// GetStats returns face buffer statistics for validation
func (fb *FaceBuffer) GetStats() map[string]int {
	stats := map[string]int{
		"total_faces":       int(fb.K * fb.Nfaces),
		"boundary_faces":    0,
		"interior_faces":    0,
		"remote_faces":      0,
		"remote_partitions": len(fb.RemoteSendIndices),
		"total_face_points": int(fb.K * fb.Nfaces * fb.Nfp),
	}

	for _, faceCode := range fb.FaceIndex {
		switch {
		case faceCode == BoundaryPlaceholder:
			stats["boundary_faces"]++
		case faceCode == RemoteFace:
			stats["remote_faces"]++
		case faceCode >= 0:
			stats["interior_faces"]++
		}
	}

	// Convert face counts to point counts
	stats["boundary_points"] = stats["boundary_faces"] * int(fb.Nfp)
	stats["interior_points"] = stats["interior_faces"] * int(fb.Nfp)
	stats["remote_points"] = stats["remote_faces"] * int(fb.Nfp)

	return stats
}

// ValidateFaceBuffer performs consistency checks
func (fb *FaceBuffer) ValidateFaceBuffer() error {
	maxMIdx := fb.K * fb.Nfaces * fb.Nfp

	// Check all face indices
	for i, faceCode := range fb.FaceIndex {
		if faceCode >= 0 {
			// Interior face - check P index is in range
			// The index points to start of face, so check if all points would be valid
			if uint32(faceCode)+fb.Nfp-1 >= maxMIdx {
				return fmt.Errorf("FaceIndex[%d] = %d would exceed M buffer size", i, faceCode)
			}
		} else if faceCode != BoundaryPlaceholder && faceCode != RemoteFace {
			// Also accept BC codes in range [-998, -1]
			if faceCode < -998 || faceCode > -1 {
				return fmt.Errorf("FaceIndex[%d] = %d is not a valid face code", i, faceCode)
			}
		}
	}

	// Check remote send indices are in range
	for partID, indices := range fb.RemoteSendIndices {
		for i, idx := range indices {
			if idx >= maxMIdx {
				return fmt.Errorf("RemoteSendIndices[%d][%d] = %d out of range",
					partID, i, idx)
			}
		}
	}

	return nil
}

// ApplyBoundaryConditions replaces boundary placeholders with specific BC codes
// This is Phase 2 of the face buffer construction
func (fb *FaceBuffer) ApplyBoundaryConditions(bcData map[int32]int32) error {
	for i, faceCode := range fb.FaceIndex {
		if faceCode == BoundaryPlaceholder {
			elem := uint32(i) / fb.Nfaces
			face := uint32(i) % fb.Nfaces

			// Look up BC type for this element/face
			// This would come from mesh BC data
			bcKey := int32(elem*fb.Nfaces + face)
			if bcType, ok := bcData[bcKey]; ok {
				if bcType < -998 || bcType > -1 {
					return fmt.Errorf("BC type %d out of range [-998, -1]", bcType)
				}
				fb.FaceIndex[i] = bcType
			} else {
				// Default to wall BC if not specified
				fb.FaceIndex[i] = -1
			}
		}
	}
	return nil
}
