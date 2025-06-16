package facebuffer

import (
	"fmt"
)

// RemoteBufferData holds runtime data for one remote partition
type RemoteBufferData struct {
	PartitionID   uint32
	SendIndices   []uint32  // Indices into local M array for send buffer assembly
	SendBuffer    []float32 // Buffer to send to this partition
	ReceiveBuffer []float32 // Buffer received from this partition
	ExpectedRecvs uint32    // Expected number of receives from this partition
}

// FaceBufferRuntime contains only the data needed for execution
// All arrays are flat and kernel-compatible
type FaceBufferRuntime struct {
	// Problem dimensions
	Nface           uint32
	Nfp             uint32
	K               uint32
	Neq             uint32
	NPart           uint32
	MyPartID        uint32
	TotalFacePoints uint32

	// Main data array: all face points in natural traversal order
	MArray []float32

	// Face classification arrays (indexed by M position)
	FaceTypes    []FaceType // Type of each face point
	BCTypes      []uint32   // BC type (if boundary)
	PartitionIDs []uint32   // Remote partition ID (if remote)

	// Local neighbor indexing
	LocalPIndices []uint32 // Maps M position → P position in same MArray

	// Remote partition data
	RemoteBuffers map[uint32]*RemoteBufferData

	// Runtime traversal counters
	localCounter   uint32
	remoteCounters map[uint32]uint32
	currentMPos    uint32
}

// Initialize allocates runtime arrays - supports golang or occa allocation
func (rt *FaceBufferRuntime) Initialize(allocType string) error {
	// Allocate main data array
	totalElements := rt.TotalFacePoints * rt.Neq

	switch allocType {
	case "golang":
		rt.MArray = make([]float32, totalElements)
	case "occa":
		// Placeholder for Occa memory allocation
		// rt.MArray = OccaMemAlloc(totalElements * sizeof(float32))
		rt.MArray = make([]float32, totalElements) // Fallback for now
	default:
		return fmt.Errorf("unknown allocation type: %s", allocType)
	}

	// Initialize remote buffers
	for _, rb := range rt.RemoteBuffers {
		rb.SendBuffer = make([]float32, len(rb.SendIndices))
		rb.ReceiveBuffer = make([]float32, rb.ExpectedRecvs)
	}

	// Initialize counters
	rt.remoteCounters = make(map[uint32]uint32)
	for partID := range rt.RemoteBuffers {
		rt.remoteCounters[partID] = 0
	}

	return nil
}

// ResetCounters resets all traversal counters for new solve iteration
func (rt *FaceBufferRuntime) ResetCounters() {
	rt.localCounter = 0
	rt.currentMPos = 0
	for partID := range rt.remoteCounters {
		rt.remoteCounters[partID] = 0
	}
}

// GetPValue returns P-side value for current M position during traversal
func (rt *FaceBufferRuntime) GetPValue(equation uint32) float32 {
	if rt.currentMPos >= rt.TotalFacePoints {
		return 0.0
	}

	faceType := rt.FaceTypes[rt.currentMPos]

	switch faceType {
	case LocalNeighbor:
		// Get P value from same M array using precomputed index
		if rt.localCounter >= uint32(len(rt.LocalPIndices)) {
			return 0.0
		}
		pIndex := rt.LocalPIndices[rt.localCounter]
		rt.localCounter++
		dataIndex := pIndex*rt.Neq + equation
		if dataIndex < uint32(len(rt.MArray)) {
			return rt.MArray[dataIndex]
		}
		return 0.0

	case RemoteNeighbor:
		// Get P value from remote partition's receive buffer
		partID := rt.PartitionIDs[rt.currentMPos]
		rb, exists := rt.RemoteBuffers[partID]
		if !exists {
			return 0.0
		}

		counter := rt.remoteCounters[partID]
		rt.remoteCounters[partID]++

		if counter < uint32(len(rb.ReceiveBuffer)) {
			return rb.ReceiveBuffer[counter]
		}
		return 0.0

	case BoundaryFace:
		// Compute boundary condition value
		bcType := rt.BCTypes[rt.currentMPos]
		return rt.computeBoundaryValue(equation, bcType)

	default:
		return 0.0
	}
}

// AdvanceMPos advances to next M position (call after GetPValue)
func (rt *FaceBufferRuntime) AdvanceMPos() {
	rt.currentMPos++
}

// GetMValue returns current M-side value
func (rt *FaceBufferRuntime) GetMValue(equation uint32) float32 {
	if rt.currentMPos >= rt.TotalFacePoints {
		return 0.0
	}
	dataIndex := rt.currentMPos*rt.Neq + equation
	if dataIndex < uint32(len(rt.MArray)) {
		return rt.MArray[dataIndex]
	}
	return 0.0
}

// SetMValue sets value in M array at specific location
func (rt *FaceBufferRuntime) SetMValue(elem, face, point, equation uint32, value float32) {
	mPos := elem*rt.Nface*rt.Nfp + face*rt.Nfp + point
	dataIndex := mPos*rt.Neq + equation
	if dataIndex < uint32(len(rt.MArray)) {
		rt.MArray[dataIndex] = value
	}
}

// AssembleSendBuffers prepares send buffers for all remote partitions
func (rt *FaceBufferRuntime) AssembleSendBuffers(equation uint32) {
	for _, rb := range rt.RemoteBuffers {
		// Assemble send buffer using indirect access via SendIndices
		for i, sendIndex := range rb.SendIndices {
			dataIndex := sendIndex*rt.Neq + equation
			if dataIndex < uint32(len(rt.MArray)) && i < len(rb.SendBuffer) {
				rb.SendBuffer[i] = rt.MArray[dataIndex]
			}
		}
	}
}

// SetReceiveBuffer stores received buffer from a remote partition
func (rt *FaceBufferRuntime) SetReceiveBuffer(partitionID uint32, buffer []float32) error {
	rb, exists := rt.RemoteBuffers[partitionID]
	if !exists {
		return fmt.Errorf("unknown remote partition %d", partitionID)
	}

	if len(buffer) != len(rb.ReceiveBuffer) {
		return fmt.Errorf("buffer size mismatch: expected %d, got %d",
			len(rb.ReceiveBuffer), len(buffer))
	}

	copy(rb.ReceiveBuffer, buffer)
	return nil
}

// GetSendBuffer returns send buffer for MPI communication
func (rt *FaceBufferRuntime) GetSendBuffer(partitionID uint32) []float32 {
	if rb, exists := rt.RemoteBuffers[partitionID]; exists {
		return rb.SendBuffer
	}
	return nil
}

// TraverseMBuffer demonstrates full M buffer traversal with P value gathering
func (rt *FaceBufferRuntime) TraverseMBuffer(equation uint32, riemannFunc func(float32, float32) float32) []float32 {
	rt.ResetCounters()
	results := make([]float32, rt.TotalFacePoints)

	for rt.currentMPos = 0; rt.currentMPos < rt.TotalFacePoints; rt.currentMPos++ {
		mValue := rt.GetMValue(equation)
		pValue := rt.GetPValue(equation)

		// Apply Riemann solver
		result := riemannFunc(mValue, pValue)
		results[rt.currentMPos] = result
	}

	rt.ResetCounters() // Clean up for next use
	return results
}

// computeBoundaryValue computes boundary condition value
func (rt *FaceBufferRuntime) computeBoundaryValue(equation, bcType uint32) float32 {
	// Simple boundary conditions for testing
	switch bcType {
	case 1: // Wall - zero value
		return 0.0
	case 2: // Inflow - fixed value
		return 1.0
	case 3: // Outflow - use interior value
		return rt.GetMValue(equation)
	default:
		return 0.0
	}
}

// GetRuntimeStatistics returns runtime memory usage statistics
func (rt *FaceBufferRuntime) GetRuntimeStatistics() map[string]uint32 {
	var totalSendBuffer, totalRecvBuffer uint32
	for _, rb := range rt.RemoteBuffers {
		totalSendBuffer += uint32(len(rb.SendBuffer))
		totalRecvBuffer += uint32(len(rb.ReceiveBuffer))
	}

	return map[string]uint32{
		"m_array_size":        uint32(len(rt.MArray)),
		"local_indices_count": uint32(len(rt.LocalPIndices)),
		"remote_partitions":   uint32(len(rt.RemoteBuffers)),
		"total_send_buffer":   totalSendBuffer,
		"total_recv_buffer":   totalRecvBuffer,
		"face_types_size":     uint32(len(rt.FaceTypes)),
		"bc_types_size":       uint32(len(rt.BCTypes)),
		"partition_ids_size":  uint32(len(rt.PartitionIDs)),
	}
}

// ValidateRuntime performs runtime structure validation
func (rt *FaceBufferRuntime) ValidateRuntime() error {
	// Check array sizes match
	if uint32(len(rt.FaceTypes)) != rt.TotalFacePoints {
		return fmt.Errorf("FaceTypes size mismatch: %d vs %d",
			len(rt.FaceTypes), rt.TotalFacePoints)
	}

	if uint32(len(rt.BCTypes)) != rt.TotalFacePoints {
		return fmt.Errorf("BCTypes size mismatch: %d vs %d",
			len(rt.BCTypes), rt.TotalFacePoints)
	}

	if uint32(len(rt.PartitionIDs)) != rt.TotalFacePoints {
		return fmt.Errorf("PartitionIDs size mismatch: %d vs %d",
			len(rt.PartitionIDs), rt.TotalFacePoints)
	}

	// Check M array size
	expectedMArraySize := rt.TotalFacePoints * rt.Neq
	if uint32(len(rt.MArray)) != expectedMArraySize {
		return fmt.Errorf("MArray size mismatch: %d vs %d",
			len(rt.MArray), expectedMArraySize)
	}

	// Validate local indices are in range
	for i, pIndex := range rt.LocalPIndices {
		if pIndex >= rt.TotalFacePoints {
			return fmt.Errorf("local P index %d at position %d out of range", pIndex, i)
		}
	}

	// Validate remote buffer consistency
	for partID, rb := range rt.RemoteBuffers {
		if rb.PartitionID != partID {
			return fmt.Errorf("remote buffer partition ID mismatch: %d vs %d",
				rb.PartitionID, partID)
		}

		if uint32(len(rb.ReceiveBuffer)) != rb.ExpectedRecvs {
			return fmt.Errorf("receive buffer size mismatch for partition %d: %d vs %d",
				partID, len(rb.ReceiveBuffer), rb.ExpectedRecvs)
		}
	}

	return nil
}

// Simple Riemann solver for testing
func SimpleRiemannSolver(mValue, pValue float32) float32 {
	return (mValue + pValue) * 0.5
}
