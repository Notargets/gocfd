// Package facebuffer implements efficient face data indexing for DG methods
package facebuffer

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/tetrahedra/tetelement"
)

// FaceType identifies the type of face connection
type FaceType uint8

const (
	BoundaryFace FaceType = iota
	InteriorFace
	RemoteFace
)

// FaceBuffer holds the runtime face indexing data
type FaceBuffer struct {
	// Dimensions
	Nfp             uint32 // Face points per face
	Nfaces          uint32 // Faces per element (4 for tetrahedra)
	K               uint32 // Number of elements in this partition
	TotalFacePoints uint32 // Total M buffer size = Nfp * Nfaces * K

	// Face classification
	FaceTypes []FaceType // Type of each face point connection

	// Local interior P indices
	LocalPIndices []uint32 // P position in M buffer for each interior point

	// Remote partition send indices
	RemoteSendIndices map[uint32][]uint32 // [partitionID] -> indices into M buffer
}

// BuildFaceBuffer creates face buffer from Element3D connectivity
// Uses ONLY VmapM/VmapP/EToE/EToP to determine face connections:
// - VmapM/VmapP: Identify boundary faces and find connecting volume nodes
// - EToE: Determine neighbor elements
// - EToP: Identify remote partitions (if present)
func BuildFaceBuffer(el *tetelement.Element3D) (*FaceBuffer, error) {
	// Extract dimensions
	nfp := uint32(el.Nfp)
	nfaces := uint32(4) // Tetrahedra
	k := uint32(el.K)
	totalFacePoints := nfp * nfaces * k

	// Initialize result
	fb := &FaceBuffer{
		Nfp:               nfp,
		Nfaces:            nfaces,
		K:                 k,
		TotalFacePoints:   totalFacePoints,
		FaceTypes:         make([]FaceType, totalFacePoints),
		LocalPIndices:     make([]uint32, 0, totalFacePoints/2), // Estimate ~half are interior
		RemoteSendIndices: make(map[uint32][]uint32),
	}

	// Determine if this is a parallel run
	isParallel := el.EToP != nil
	var myPartID int
	if isParallel {
		// Find this partition's ID (all local elements should have same partition ID)
		myPartID = el.EToP[0]
	}

	// Process each face point in natural traversal order
	for idx := uint32(0); idx < totalFacePoints; idx++ {
		// VmapM and VmapP tell us the volume nodes
		volM := el.VmapM[idx]
		volP := el.VmapP[idx]

		// Check if boundary (self-reference)
		if volM == volP {
			fb.FaceTypes[idx] = BoundaryFace
			continue
		}

		// Interior or remote connection
		// Find which element and face we're on
		elem := idx / (nfaces * nfp)
		face := (idx % (nfaces * nfp)) / nfp

		// Get neighbor element from EToE
		neighborElem := el.EToE[elem][face]

		// Check if this is actually a self-reference (boundary)
		if neighborElem == int(elem) {
			fb.FaceTypes[idx] = BoundaryFace
			continue
		}

		// Check if remote connection
		if isParallel {
			// Check if neighbor is outside our partition's element range
			if uint32(neighborElem) >= k {
				// Neighbor is in a different partition
				fb.FaceTypes[idx] = RemoteFace
				continue
			}

			// Check if neighbor is in a different partition ID
			if neighborElem < len(el.EToP) && el.EToP[neighborElem] != myPartID {
				fb.FaceTypes[idx] = RemoteFace
				continue
			}
		}

		// Interior connection within this partition
		fb.FaceTypes[idx] = InteriorFace

		// Find P position in local M array
		pPos := findPPositionInNeighbor(el, idx, uint32(neighborElem))
		fb.LocalPIndices = append(fb.LocalPIndices, pPos)
	}

	return fb, nil
}

// findPPositionInNeighbor finds where the P point is located in neighbor's M array
// Note: This only works for neighbors within the same partition
func findPPositionInNeighbor(el *tetelement.Element3D, mIdx uint32, neighborElem uint32) uint32 {
	// The P volume node we're looking for
	targetVolNode := el.VmapP[mIdx]

	// Search through neighbor element's face points
	nfp := uint32(el.Nfp)
	nfaces := uint32(4)
	k := uint32(el.K)

	// Make sure neighbor element is within our partition
	if neighborElem >= k {
		panic(fmt.Sprintf("Neighbor element %d is outside partition range [0,%d)", neighborElem, k))
	}

	// Search all face points of neighbor element
	for face := uint32(0); face < nfaces; face++ {
		for point := uint32(0); point < nfp; point++ {
			neighborIdx := neighborElem*nfaces*nfp + face*nfp + point
			if neighborIdx < uint32(len(el.VmapM)) && el.VmapM[neighborIdx] == targetVolNode {
				return neighborIdx
			}
		}
	}

	// This should never happen with valid connectivity
	panic(fmt.Sprintf("No matching P position found for M index %d", mIdx))
}

// GetStats returns face buffer statistics for validation
func (fb *FaceBuffer) GetStats() map[string]int {
	stats := map[string]int{
		"total_face_points": int(fb.TotalFacePoints),
		"boundary_points":   0,
		"interior_points":   0,
		"remote_points":     0,
		"remote_partitions": len(fb.RemoteSendIndices),
	}

	for _, ft := range fb.FaceTypes {
		switch ft {
		case BoundaryFace:
			stats["boundary_points"]++
		case InteriorFace:
			stats["interior_points"]++
		case RemoteFace:
			stats["remote_points"]++
		}
	}

	return stats
}

// ValidateFaceBuffer performs consistency checks
func (fb *FaceBuffer) ValidateFaceBuffer() error {
	// Check that interior points equals length of LocalPIndices
	interiorCount := 0
	for _, ft := range fb.FaceTypes {
		if ft == InteriorFace {
			interiorCount++
		}
	}

	if interiorCount != len(fb.LocalPIndices) {
		return fmt.Errorf("interior point count %d doesn't match LocalPIndices length %d",
			interiorCount, len(fb.LocalPIndices))
	}

	// Check all local P indices are in range
	for i, pIdx := range fb.LocalPIndices {
		if pIdx >= fb.TotalFacePoints {
			return fmt.Errorf("LocalPIndices[%d] = %d out of range", i, pIdx)
		}
	}

	// Check remote indices are in range
	for partID, indices := range fb.RemoteSendIndices {
		for i, idx := range indices {
			if idx >= fb.TotalFacePoints {
				return fmt.Errorf("RemoteSendIndices[%d][%d] = %d out of range",
					partID, i, idx)
			}
		}
	}

	return nil
}
