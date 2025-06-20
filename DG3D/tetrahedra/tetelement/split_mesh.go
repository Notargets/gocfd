package tetelement

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/DG3D/tetrahedra/gonudg"
	"github.com/notargets/gocfd/utils"
)

// SplitMesh splits the global mesh into partition Element3D structures
func (el *Element3D) SplitMesh() (err error) {
	// Normalize EToP to 0-based
	el.EToP = normalizeToZeroBased(el.EToP)

	// Find number of partitions
	numPartitions := 0
	for _, partID := range el.EToP {
		if partID >= numPartitions {
			numPartitions = partID + 1
		}
	}

	if numPartitions == 0 {
		err = fmt.Errorf("no partitions found")
		return
	}

	// Initialize partition data structures
	el.PEToE = make(map[int][]int)
	partEToV := make(map[int][][]int)
	partBCs := make(map[int]map[string][]mesh.BoundaryElement)

	// Build PEToE and partition element data
	for elemID := 0; elemID < el.Mesh.NumElements; elemID++ {
		partID := el.EToP[elemID]

		// Track element mapping
		el.PEToE[partID] = append(el.PEToE[partID], elemID)

		// Copy element connectivity (keeping global vertex IDs!)
		partEToV[partID] = append(partEToV[partID], el.Mesh.EtoV[elemID])
	}

	// Transform BCs to local element indices
	for partID := 0; partID < numPartitions; partID++ {
		partBCs[partID] = make(map[string][]mesh.BoundaryElement)

		// Create global to local element mapping for this partition
		globalToLocal := make(map[int]int)
		for localIdx, globalIdx := range el.PEToE[partID] {
			globalToLocal[globalIdx] = localIdx
		}

		// Transform each BC
		for bcTag, bcs := range el.Mesh.BoundaryElements {
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
	el.Split = make([]*Element3D, numPartitions)
	for partID := 0; partID < numPartitions; partID++ {
		// Create DG3D for this partition
		dg3d, err2 := gonudg.NewDG3D(el.N, el.VX, el.VY, el.VZ, partEToV[partID])
		if err2 != nil {
			err = fmt.Errorf("failed to create DG3D for partition %d: %v",
				partID, err2)
			return
		}

		// Create Element3D
		elNew := &Element3D{
			K:    len(partEToV[partID]),
			DG3D: dg3d,
			Mesh: nil, // As requested
			EToP: nil, // Not needed for individual partitions
		}

		// Initialize BC maps
		elNew.BCMaps = &BCFaceMap{
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
				elNew.BCMaps.FaceBC[faceKey] = bcType
			}
		}

		// Map nodes to BC types (similar to mapNodesToBC in bc_mapping.go)
		if err = mapNodesToBCForPartition(elNew); err != nil {
			err = fmt.Errorf("failed to map nodes to BC for partition %d: %v",
				partID, err)
			return
		}

		el.Split[partID] = elNew
	}

	return
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
