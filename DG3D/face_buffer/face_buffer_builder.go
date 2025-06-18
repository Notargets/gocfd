// Package facebuffer implements efficient face data indexing for DG methods
// Separates build-time complexity from runtime execution structures
package facebuffer

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/tetrahedra/element3d"
	"sort"
)

// FaceType identifies the type of face connection
type FaceType uint8

const (
	BoundaryFace FaceType = iota
	LocalNeighbor
	RemoteNeighbor
)

// BuildTimeConnection represents a face connection during build phase
type BuildTimeConnection struct {
	MPos        uint32   // M-side position in natural traversal order
	Type        FaceType // Type of connection
	BCType      uint32   // Boundary condition type (if boundary)
	PartitionID uint32   // Remote partition ID (if remote)
	PPos        uint32   // P-side position (local index or remote index)
}

// RemotePartitionInfo holds build-time info for one remote partition
type RemotePartitionInfo struct {
	PartitionID   uint32
	Connections   []BuildTimeConnection
	SendIndices   []uint32 // Will become runtime SendIndices
	ExpectedRecvs uint32   // Number of values expected from this partition
}

// FaceBufferBuilder handles the complex build-time logic
type FaceBufferBuilder struct {
	// Problem dimensions
	Nface    uint32 // Number of faces per element
	Nfp      uint32 // Number of face points per face
	K        uint32 // Number of elements in this partition
	Neq      uint32 // Number of equations
	NPart    uint32 // Total number of partitions
	MyPartID uint32 // This partition's ID

	// Build-time data structures
	connections        []BuildTimeConnection
	remotePartitions   map[uint32]*RemotePartitionInfo
	totalFacePoints    uint32
	interiorPointCount uint32
}

// NewFaceBufferBuilder creates a new builder from Element3D
func NewFaceBufferBuilder(el *element3d.Element3D, neq uint32) *FaceBufferBuilder {
	// Extract dimensions from Element3D
	nface := uint32(4) // Tetrahedra have 4 faces
	nfp := uint32(el.Nfp)
	k := uint32(el.K)

	// Determine partition info
	var npart, myPartID uint32
	if el.EToP != nil {
		// Partitioned mesh - find max partition ID and determine local partition
		maxPartID := uint32(0)
		myPartID = uint32(el.EToP[0]) // Assume all local elements have same partition ID
		for _, partID := range el.EToP {
			if uint32(partID) > maxPartID {
				maxPartID = uint32(partID)
			}
		}
		npart = maxPartID + 1
	} else {
		// Non-partitioned mesh - single partition
		npart = 1
		myPartID = 0
	}

	return &FaceBufferBuilder{
		Nface:            nface,
		Nfp:              nfp,
		K:                k,
		Neq:              neq,
		NPart:            npart,
		MyPartID:         myPartID,
		connections:      make([]BuildTimeConnection, 0),
		remotePartitions: make(map[uint32]*RemotePartitionInfo),
		totalFacePoints:  nface * nfp * k,
	}
}

// AddBoundaryFace adds a boundary condition face point
func (fb *FaceBufferBuilder) AddBoundaryFace(elem, face, point uint32, bcType uint32) error {
	mPos := fb.computeMPos(elem, face, point)
	if mPos >= fb.totalFacePoints {
		return fmt.Errorf("M position %d exceeds total face points %d", mPos, fb.totalFacePoints)
	}

	conn := BuildTimeConnection{
		MPos:   mPos,
		Type:   BoundaryFace,
		BCType: bcType,
	}
	fb.connections = append(fb.connections, conn)
	return nil
}

// AddLocalNeighbor adds a local interior face connection
func (fb *FaceBufferBuilder) AddLocalNeighbor(elemM, faceM, pointM, elemP, faceP, pointP uint32) error {
	mPos := fb.computeMPos(elemM, faceM, pointM)
	pPos := fb.computeMPos(elemP, faceP, pointP) // P position in same M array

	if mPos >= fb.totalFacePoints || pPos >= fb.totalFacePoints {
		return fmt.Errorf("face positions out of range: M=%d, P=%d, total=%d",
			mPos, pPos, fb.totalFacePoints)
	}

	conn := BuildTimeConnection{
		MPos: mPos,
		Type: LocalNeighbor,
		PPos: pPos,
	}
	fb.connections = append(fb.connections, conn)
	fb.interiorPointCount++
	return nil
}

// AddRemoteNeighbor adds a remote face connection
func (fb *FaceBufferBuilder) AddRemoteNeighbor(elemM, faceM, pointM uint32,
	remotePartID, remoteElemP, remoteFaceP, remotePointP uint32) error {

	mPos := fb.computeMPos(elemM, faceM, pointM)
	if mPos >= fb.totalFacePoints {
		return fmt.Errorf("M position %d exceeds total face points %d", mPos, fb.totalFacePoints)
	}

	if remotePartID == fb.MyPartID {
		return fmt.Errorf("remote partition ID cannot be same as local partition %d", fb.MyPartID)
	}

	// Compute P position in remote partition's M array indexing
	remotePPos := fb.computeMPos(remoteElemP, remoteFaceP, remotePointP)

	conn := BuildTimeConnection{
		MPos:        mPos,
		Type:        RemoteNeighbor,
		PartitionID: remotePartID,
		PPos:        remotePPos,
	}
	fb.connections = append(fb.connections, conn)

	// Track remote partition info
	if _, exists := fb.remotePartitions[remotePartID]; !exists {
		fb.remotePartitions[remotePartID] = &RemotePartitionInfo{
			PartitionID: remotePartID,
			Connections: make([]BuildTimeConnection, 0),
			SendIndices: make([]uint32, 0),
		}
	}
	fb.remotePartitions[remotePartID].Connections = append(
		fb.remotePartitions[remotePartID].Connections, conn)

	return nil
}

// ProcessRemoteIndices processes the remote indices that other partitions will send us
// This is called when we receive the send indices from remote partitions
func (fb *FaceBufferBuilder) ProcessRemoteIndices(remotePartID uint32, theirSendIndices []uint32) error {
	rpi, exists := fb.remotePartitions[remotePartID]
	if !exists {
		return fmt.Errorf("unknown remote partition %d", remotePartID)
	}

	// Store the indices they will use to send us data
	rpi.SendIndices = make([]uint32, len(theirSendIndices))
	copy(rpi.SendIndices, theirSendIndices)
	rpi.ExpectedRecvs = uint32(len(theirSendIndices))

	return nil
}

// BuildFromElement3D automatically constructs face buffer from Element3D connectivity
func (fb *FaceBufferBuilder) BuildFromElement3D(el *element3d.Element3D) error {
	// Ensure connectivity is built
	if el.EToE == nil || el.EToF == nil {
		if el.ConnectivityArrays == nil {
			el.ConnectivityArrays = el.Connect3D()
		}
	}

	nfaces := 4 // Tetrahedra

	// Process all faces in natural traversal order
	for k := 0; k < el.K; k++ {
		for f := 0; f < nfaces; f++ {
			for fp := 0; fp < el.Nfp; fp++ {
				neighbor := el.EToE[k][f]

				if neighbor == -1 {
					// Boundary face - determine BC type from mesh data
					bcType := uint32(1) // Default wall BC
					if len(el.BCType) > k*nfaces+f {
						bcType = uint32(el.BCType[k*nfaces+f])
					}

					// Check if this is a partition boundary
					if bcType == uint32(element3d.BCPartitionBoundary) {
						// For now, treat partition boundaries as regular boundaries
						// In a distributed implementation, these would trigger special handling
					}

					err := fb.AddBoundaryFace(uint32(k), uint32(f), uint32(fp), bcType)
					if err != nil {
						return fmt.Errorf("error adding boundary face: %v", err)
					}

				} else if el.EToP == nil || el.EToP[neighbor] == el.EToP[k] {
					// Local interior neighbor (same partition or non-partitioned)
					neighborFace := el.EToF[k][f]

					// CRITICAL FIX: Use VmapP to get the correct face point mapping
					// VmapP contains the proper face-to-face node correspondence
					faceIdx := k*nfaces*el.Nfp + f*el.Nfp + fp
					vmapP := el.VmapP[faceIdx]

					// Extract element and node from VmapP
					neighborElem := vmapP / el.Np
					neighborNode := vmapP % el.Np

					// Verify the neighbor element matches
					if neighborElem != neighbor {
						return fmt.Errorf("VmapP inconsistency: expected element %d, got %d at k=%d, f=%d, fp=%d",
							neighbor, neighborElem, k, f, fp)
					}

					// Find which face point on the neighbor face has this node
					var neighborFP int = -1
					for nfp := 0; nfp < el.Nfp; nfp++ {
						if el.Fmask[neighborFace][nfp] == neighborNode {
							neighborFP = nfp
							break
						}
					}

					if neighborFP == -1 {
						return fmt.Errorf("could not find face point for node %d on face %d of element %d",
							neighborNode, neighborFace, neighbor)
					}

					// Now we have the correct face point correspondence
					err := fb.AddLocalNeighbor(
						uint32(k), uint32(f), uint32(fp),
						uint32(neighbor), uint32(neighborFace), uint32(neighborFP))
					if err != nil {
						return fmt.Errorf("error adding local neighbor: %v", err)
					}

				} else {
					// Remote neighbor (different partition)
					remotePartID := uint32(el.EToP[neighbor])
					neighborFace := el.EToF[k][f]

					// For remote neighbors, we also need to use VmapP to get the correct mapping
					faceIdx := k*nfaces*el.Nfp + f*el.Nfp + fp
					vmapP := el.VmapP[faceIdx]

					// Extract node from VmapP
					neighborNode := vmapP % el.Np

					// Find which face point on the neighbor face corresponds to this node
					var neighborFP int = -1
					for nfp := 0; nfp < el.Nfp; nfp++ {
						if el.Fmask[neighborFace][nfp] == neighborNode {
							neighborFP = nfp
							break
						}
					}

					if neighborFP == -1 {
						// For remote neighbors, the face point might not be directly findable
						// due to different element numbering in remote partition
						// Use the same face point index as a fallback
						neighborFP = fp
					}

					err := fb.AddRemoteNeighbor(
						uint32(k), uint32(f), uint32(fp),
						remotePartID, uint32(neighbor), uint32(neighborFace), uint32(neighborFP))
					if err != nil {
						return fmt.Errorf("error adding remote neighbor: %v", err)
					}
				}
			}
		}
	}

	return fb.ValidateBuild()
}

// Build creates the runtime structure from Element3D
func (fb *FaceBufferBuilder) Build(el *element3d.Element3D) (*FaceBufferRuntime, error) {
	// Auto-build connections from Element3D if no manual connections added
	if len(fb.connections) == 0 {
		err := fb.BuildFromElement3D(el)
		if err != nil {
			return nil, fmt.Errorf("failed to build from Element3D: %v", err)
		}
	}

	// Sort connections by M position for natural traversal order
	sort.Slice(fb.connections, func(i, j int) bool {
		return fb.connections[i].MPos < fb.connections[j].MPos
	})

	// Validate we have connections for all face points
	if uint32(len(fb.connections)) != fb.totalFacePoints {
		return nil, fmt.Errorf("connection count %d does not match total face points %d",
			len(fb.connections), fb.totalFacePoints)
	}

	// Build runtime arrays
	faceTypes := make([]FaceType, fb.totalFacePoints)
	bcTypes := make([]uint32, fb.totalFacePoints)
	partitionIDs := make([]uint32, fb.totalFacePoints)
	localPIndices := make([]uint32, 0, fb.interiorPointCount)

	// Process connections in M traversal order
	for i, conn := range fb.connections {
		if uint32(i) != conn.MPos {
			return nil, fmt.Errorf("connection order mismatch at position %d", i)
		}

		faceTypes[i] = conn.Type

		switch conn.Type {
		case BoundaryFace:
			bcTypes[i] = conn.BCType

		case LocalNeighbor:
			localPIndices = append(localPIndices, conn.PPos)

		case RemoteNeighbor:
			partitionIDs[i] = conn.PartitionID
		}
	}

	// Build remote partition data
	remoteBuffers := make(map[uint32]*RemoteBufferData)
	for partID, rpi := range fb.remotePartitions {
		remoteBuffers[partID] = &RemoteBufferData{
			PartitionID:   partID,
			SendIndices:   make([]uint32, len(rpi.SendIndices)),
			ExpectedRecvs: rpi.ExpectedRecvs,
		}
		copy(remoteBuffers[partID].SendIndices, rpi.SendIndices)
	}

	runtime := &FaceBufferRuntime{
		// Dimensions
		Nface:           fb.Nface,
		Nfp:             fb.Nfp,
		K:               fb.K,
		Neq:             fb.Neq,
		NPart:           fb.NPart,
		MyPartID:        fb.MyPartID,
		TotalFacePoints: fb.totalFacePoints,

		// Runtime arrays
		FaceTypes:     faceTypes,
		BCTypes:       bcTypes,
		PartitionIDs:  partitionIDs,
		LocalPIndices: localPIndices,
		RemoteBuffers: remoteBuffers,
	}

	return runtime, runtime.Initialize("golang")
}

// computeMPos calculates M position in natural traversal order
func (fb *FaceBufferBuilder) computeMPos(elem, face, point uint32) uint32 {
	return elem*fb.Nface*fb.Nfp + face*fb.Nfp + point
}

// GetBuildStatistics returns build-time statistics
func (fb *FaceBufferBuilder) GetBuildStatistics() map[string]uint32 {
	stats := map[string]uint32{
		"total_face_points": fb.totalFacePoints,
		"interior_points":   fb.interiorPointCount,
		"total_connections": uint32(len(fb.connections)),
		"remote_partitions": uint32(len(fb.remotePartitions)),
	}

	var boundaryPoints, remotePoints, partitionBoundaryPoints, domainBoundaryPoints uint32
	for _, conn := range fb.connections {
		switch conn.Type {
		case BoundaryFace:
			boundaryPoints++
			if conn.BCType == uint32(element3d.BCPartitionBoundary) {
				partitionBoundaryPoints++
			} else {
				domainBoundaryPoints++
			}
		case RemoteNeighbor:
			remotePoints++
		}
	}

	stats["boundary_points"] = boundaryPoints
	stats["remote_points"] = remotePoints
	stats["partition_boundary_points"] = partitionBoundaryPoints
	stats["domain_boundary_points"] = domainBoundaryPoints

	return stats
}

// ValidateBuild performs validation checks
func (fb *FaceBufferBuilder) ValidateBuild() error {
	// Check all positions are covered
	positionsSeen := make(map[uint32]bool)
	for _, conn := range fb.connections {
		if positionsSeen[conn.MPos] {
			return fmt.Errorf("duplicate M position %d", conn.MPos)
		}
		positionsSeen[conn.MPos] = true
	}

	// Check for gaps
	for i := uint32(0); i < fb.totalFacePoints; i++ {
		if !positionsSeen[i] {
			return fmt.Errorf("missing connection for M position %d", i)
		}
	}

	// Validate local P indices are in range
	for _, conn := range fb.connections {
		if conn.Type == LocalNeighbor && conn.PPos >= fb.totalFacePoints {
			return fmt.Errorf("local P position %d out of range", conn.PPos)
		}
	}

	return nil
}
