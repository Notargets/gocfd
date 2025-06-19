package tetelement

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/DG3D/tetrahedra/gonudg"
	"github.com/notargets/gocfd/utils"
)

// MeshSplitter handles splitting a global mesh into partition-local Element3D structures
type MeshSplitter struct {
	// Inputs
	GlobalMesh *mesh.Mesh
	VX, VY, VZ []float64 // Global vertex coordinates
	EToP       []int     // Element to partition mapping
	Order      int       // Polynomial order for DG3D

	// Outputs
	PartitionElements []*Element3D  // Split elements with DG3D initialized
	PEToE             map[int][]int // Partition element to global element mapping
	EToP_0based       []int         // Normalized to 0-based indexing
}

// NewMeshSplitter creates a new mesh splitter
func NewMeshSplitter(order int, globalMesh *mesh.Mesh, vx, vy, vz []float64, eToP []int) *MeshSplitter {
	return &MeshSplitter{
		Order:      order,
		GlobalMesh: globalMesh,
		VX:         vx,
		VY:         vy,
		VZ:         vz,
		EToP:       eToP,
		PEToE:      make(map[int][]int),
	}
}

// SplitMesh splits the global mesh into partition Element3D structures
func (ms *MeshSplitter) SplitMesh() ([]*Element3D, map[int][]int, error) {
	// Normalize EToP to 0-based
	ms.EToP_0based = normalizeToZeroBased(ms.EToP)

	// Find number of partitions
	numPartitions := 0
	for _, partID := range ms.EToP_0based {
		if partID >= numPartitions {
			numPartitions = partID + 1
		}
	}

	if numPartitions == 0 {
		return nil, nil, fmt.Errorf("no partitions found")
	}

	// Initialize partition data structures
	ms.PEToE = make(map[int][]int)
	partEToV := make(map[int][][]int)
	partElementTypes := make(map[int][]utils.ElementType)
	partBCs := make(map[int]map[string][]mesh.BoundaryElement)

	// Build PEToE and partition element data
	for elemID := 0; elemID < ms.GlobalMesh.NumElements; elemID++ {
		partID := ms.EToP_0based[elemID]

		// Track element mapping
		ms.PEToE[partID] = append(ms.PEToE[partID], elemID)

		// Copy element connectivity (keeping global vertex IDs!)
		partEToV[partID] = append(partEToV[partID], ms.GlobalMesh.EtoV[elemID])

		// Copy element type
		partElementTypes[partID] = append(partElementTypes[partID], ms.GlobalMesh.ElementTypes[elemID])
	}

	// Transform BCs to local element indices
	for partID := 0; partID < numPartitions; partID++ {
		partBCs[partID] = make(map[string][]mesh.BoundaryElement)

		// Create global to local element mapping for this partition
		globalToLocal := make(map[int]int)
		for localIdx, globalIdx := range ms.PEToE[partID] {
			globalToLocal[globalIdx] = localIdx
		}

		// Transform each BC
		for bcTag, bcs := range ms.GlobalMesh.BoundaryElements {
			for _, bc := range bcs {
				// Check if this BC's parent element belongs to this partition
				if localElemID, exists := globalToLocal[bc.ParentElement]; exists {
					// Create new BC with local element index
					localBC := mesh.BoundaryElement{
						ElementType:   bc.ElementType,
						Nodes:         bc.Nodes, // Keep global node IDs
						ParentElement: localElemID,
						ParentFace:    bc.ParentFace, // Face ID unchanged
					}
					partBCs[partID][bcTag] = append(partBCs[partID][bcTag], localBC)
				}
			}
		}
	}

	// Create partition Element3D structures
	ms.PartitionElements = make([]*Element3D, numPartitions)
	for partID := 0; partID < numPartitions; partID++ {
		// Create DG3D for this partition
		dg3d, err := gonudg.NewDG3D(ms.Order, ms.VX, ms.VY, ms.VZ, partEToV[partID])
		if err != nil {
			return nil, nil, fmt.Errorf("failed to create DG3D for partition %d: %v", partID, err)
		}

		// Create Element3D
		el := &Element3D{
			K:    len(partEToV[partID]),
			DG3D: dg3d,
			Mesh: nil, // As requested
			EToP: nil, // Not needed for individual partitions
		}

		// Initialize BC maps
		el.BCMaps = &BCFaceMap{
			NodeBC: make(map[int]utils.BCType),
			FaceBC: make(map[FaceKey]utils.BCType),
		}

		// Populate FaceBC from transformed boundary elements
		for bcName, boundaryElems := range partBCs[partID] {
			bcType := utils.ParseBCName(bcName)

			for _, be := range boundaryElems {
				faceKey := FaceKey{
					Element: be.ParentElement,
					Face:    be.ParentFace,
				}
				el.BCMaps.FaceBC[faceKey] = bcType
			}
		}

		// Map nodes to BC types (similar to mapNodesToBC in bc_mapping.go)
		if err := mapNodesToBCForPartition(el); err != nil {
			return nil, nil, fmt.Errorf("failed to map nodes to BC for partition %d: %v", partID, err)
		}

		ms.PartitionElements[partID] = el
	}

	return ms.PartitionElements, ms.PEToE, nil
}

// mapNodesToBCForPartition maps boundary node indices to their BC types for a partition
func mapNodesToBCForPartition(el *Element3D) error {
	// Similar to mapNodesToBC in bc_mapping.go
	for i, mapBIdx := range el.DG3D.MapB {
		// Calculate element and face from the mapBIdx
		nodesPerFace := el.DG3D.Nfp
		facesPerElem := el.DG3D.Nfaces
		nodesPerElem := nodesPerFace * facesPerElem

		elem := mapBIdx / nodesPerElem
		localIdx := mapBIdx % nodesPerElem
		face := localIdx / nodesPerFace

		// Get BC type for this face
		faceKey := FaceKey{Element: elem, Face: face}
		bcType, exists := el.BCMaps.FaceBC[faceKey]
		if !exists {
			// Default BC type if not specified
			bcType = utils.BCWall
		}

		// Store BC type for this boundary node
		el.BCMaps.NodeBC[i] = bcType
	}

	return nil
}

// normalizeToZeroBased converts arbitrary partition IDs to 0-based sequential IDs
func normalizeToZeroBased(eToP []int) []int {
	if len(eToP) == 0 {
		return []int{}
	}

	// Find unique partition IDs
	uniqueIDs := make(map[int]bool)
	for _, id := range eToP {
		uniqueIDs[id] = true
	}

	// Create sorted list of unique IDs
	var sortedIDs []int
	for id := range uniqueIDs {
		sortedIDs = append(sortedIDs, id)
	}

	// Simple sort
	for i := 0; i < len(sortedIDs); i++ {
		for j := i + 1; j < len(sortedIDs); j++ {
			if sortedIDs[j] < sortedIDs[i] {
				sortedIDs[i], sortedIDs[j] = sortedIDs[j], sortedIDs[i]
			}
		}
	}

	// Create mapping from old to new IDs
	mapping := make(map[int]int)
	for newID, oldID := range sortedIDs {
		mapping[oldID] = newID
	}

	// Apply mapping
	result := make([]int, len(eToP))
	for i, oldID := range eToP {
		result[i] = mapping[oldID]
	}

	return result
}
