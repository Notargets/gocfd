// Package facebuffer implements efficient face data indexing for DG methods
// Separates build-time complexity from runtime execution structures
package facebuffer

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/tetrahedra/tetelement"
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
func NewFaceBufferBuilder(el *tetelement.Element3D, neq uint32) *FaceBufferBuilder {
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
func (fb *FaceBufferBuilder) BuildFromElement3D(el *tetelement.Element3D) error {
	// TODO: Implement this
	return fb.ValidateBuild()
}

// Build creates the runtime structure from Element3D
func (fb *FaceBufferBuilder) Build(el *tetelement.Element3D) (*FaceBufferRuntime, error) {
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
			// if conn.BCType == uint32(tetelement.BCPartitionBoundary) {
			// 	partitionBoundaryPoints++
			// } else {
			// 	domainBoundaryPoints++
			// }
			// TODO: Implement
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
