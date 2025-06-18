package tetelement

import (
	"fmt"
)

// GlobalToLocal represents the mapping from global element ID to partition-local coordinates
type GlobalToLocal struct {
	mapping map[int]struct {
		PartitionID int
		LocalIndex  int
	}
}

// NewGlobalToLocal builds the mapping from the partition assignments
func NewGlobalToLocal(EToP []int) *GlobalToLocal {
	g2l := &GlobalToLocal{
		mapping: make(map[int]struct {
			PartitionID int
			LocalIndex  int
		}),
	}

	// First pass: count elements per partition to assign local indices
	partitionCounts := make(map[int]int)

	// Second pass: build the mapping
	for globalIdx, partID := range EToP {
		localIdx := partitionCounts[partID]
		g2l.mapping[globalIdx] = struct {
			PartitionID int
			LocalIndex  int
		}{partID, localIdx}
		partitionCounts[partID]++
	}

	return g2l
}

// Get returns the partition ID and local index for a global element
func (g2l *GlobalToLocal) Get(globalIdx int) (partitionID, localIndex int, ok bool) {
	if info, exists := g2l.mapping[globalIdx]; exists {
		return info.PartitionID, info.LocalIndex, true
	}
	return -1, -1, false
}

// GetLocalIndex returns just the local index for a global element
func (g2l *GlobalToLocal) GetLocalIndex(globalIdx int) int {
	if info, exists := g2l.mapping[globalIdx]; exists {
		return info.LocalIndex
	}
	return -1
}

// SplitByPartition splits the Element3D into separate Element3D instances for each partition
// This should be called after EToP is populated and before face_buffer or other parallel processing
func (el *Element3D) SplitByPartition() error {
	if el.EToP == nil {
		// Nothing to split - non-partitioned mesh
		return nil
	}

	// Build global to local mapping for ALL partitions
	g2l := NewGlobalToLocal(el.EToP)

	// Count elements per partition and create partition mapping
	partitionCount := make(map[int]int)
	partitionElements := make(map[int][]int) // partition -> list of element indices

	for elemID, partID := range el.EToP {
		partitionCount[partID]++
		partitionElements[partID] = append(partitionElements[partID], elemID)
	}

	// Create split Element3D for each partition
	el.SplitElement3D = make([]*Element3D, 0, len(partitionCount))

	for partID, elemList := range partitionElements {
		splitEl, err := el.createPartitionElement(partID, elemList, g2l)
		if err != nil {
			return fmt.Errorf("failed to create partition %d: %v", partID, err)
		}
		el.SplitElement3D = append(el.SplitElement3D, splitEl)
	}

	return nil
}

// createPartitionElement creates a new Element3D for a specific partition
func (el *Element3D) createPartitionElement(partitionID int, elemIndices []int, g2l *GlobalToLocal) (*Element3D, error) {
	localK := len(elemIndices)

	// Create global to local element mapping for this partition
	globalToLocal := make(map[int]int)
	for localIdx, globalIdx := range elemIndices {
		globalToLocal[globalIdx] = localIdx
	}

	// Create new Element3D for this partition
	partEl := &Element3D{
		K:    localK,
		Mesh: el.Mesh, // Reference to original mesh
	}

	return partEl, nil
}

// GetPartition returns the Element3D for a specific partition ID
func (el *Element3D) GetPartition(partitionID int) (*Element3D, error) {
	if el.SplitElement3D == nil {
		return nil, fmt.Errorf("mesh has not been split by partition")
	}

	// Find partition with matching ID
	for _, partEl := range el.SplitElement3D {
		if partEl.EToP != nil && len(partEl.EToP) > 0 && partEl.EToP[0] == partitionID {
			return partEl, nil
		}
	}

	return nil, fmt.Errorf("partition %d not found", partitionID)
}
