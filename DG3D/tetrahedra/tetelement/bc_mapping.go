package tetelement

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
)

// BCFaceMap maps boundary face nodes to their BC types
type BCFaceMap struct {
	NodeBC map[int]utils.BCType     // Maps MapB indices (0 to len(MapB)-1) to BC types
	FaceBC map[FaceKey]utils.BCType // Maps element/face pairs to BC types
}

// FaceKey uniquely identifies a face in the mesh
type FaceKey struct {
	Element int
	Face    int
}

// BuildBCMaps processes the mesh boundary information and creates BC mappings
// This should be called after mesh construction and BuildMaps3D
func (el *Element3D) BuildBCMaps() error {
	if el.Mesh == nil {
		return fmt.Errorf("mesh not initialized")
	}

	if el.DG3D == nil {
		return fmt.Errorf("DG3D not initialized")
	}

	// Initialize BC maps
	el.BCMaps = &BCFaceMap{
		NodeBC: make(map[int]utils.BCType),
		FaceBC: make(map[FaceKey]utils.BCType),
	}

	// First, build element face to BC type mapping from mesh boundary elements
	if err := el.buildFaceBCMap(); err != nil {
		return fmt.Errorf("failed to build face BC map: %v", err)
	}

	// Then, map boundary nodes to BC types
	if err := el.mapNodesToBC(); err != nil {
		return fmt.Errorf("failed to map nodes to BC: %v", err)
	}

	return nil
}

// buildFaceBCMap creates a mapping from element/face pairs to BC types
func (el *Element3D) buildFaceBCMap() error {
	if el.Mesh.BoundaryElements == nil {
		// No boundary conditions defined
		return nil
	}

	// Process each boundary condition
	for bcName, boundaryElems := range el.Mesh.BoundaryElements {
		// Convert BC name to BCType
		bcType := utils.ParseBCName(bcName)

		// Map each boundary element to its parent element/face
		for _, be := range boundaryElems {
			// Create face key
			faceKey := FaceKey{
				Element: be.ParentElement,
				Face:    be.ParentFace,
			}

			// Store BC type for this face
			el.BCMaps.FaceBC[faceKey] = bcType
		}
	}

	return nil
}

// mapNodesToBC maps boundary node indices to their BC types
func (el *Element3D) mapNodesToBC() error {
	// MapB contains indices into the face arrays where boundary nodes are located
	// VmapB contains the corresponding volume node indices

	// For each boundary node index in MapB
	for i, mapBIdx := range el.DG3D.MapB {
		// MapB[i] gives us the position in the face arrays
		// We need to figure out which element and face this belongs to

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

// GetBCType returns the BC type for a boundary node
// boundaryNodeIdx is the index into MapB/VmapB arrays (not the volume node index)
func (el *Element3D) GetBCType(boundaryNodeIdx int) utils.BCType {
	if el.BCMaps == nil {
		return utils.BCNone
	}

	if bcType, ok := el.BCMaps.NodeBC[boundaryNodeIdx]; ok {
		return bcType
	}

	// Default to wall if not found
	return utils.BCWall
}

// GetFaceBCType returns the BC type for an element face
func (el *Element3D) GetFaceBCType(element, face int) utils.BCType {
	if el.BCMaps == nil {
		return utils.BCNone
	}

	faceKey := FaceKey{Element: element, Face: face}
	if bcType, ok := el.BCMaps.FaceBC[faceKey]; ok {
		return bcType
	}

	return utils.BCNone
}

// GetBoundaryStatistics returns statistics about boundary conditions
func (el *Element3D) GetBoundaryStatistics() map[utils.BCType]int {
	stats := make(map[utils.BCType]int)

	if el.BCMaps == nil {
		return stats
	}

	// Count nodes by BC type
	for _, bcType := range el.BCMaps.NodeBC {
		stats[bcType]++
	}

	return stats
}

// ValidateBCMapping performs validation checks on the BC mapping
func (el *Element3D) ValidateBCMapping() error {
	if el.BCMaps == nil {
		return fmt.Errorf("BC maps not initialized")
	}

	// Check that all boundary nodes have a BC type
	for i := range el.DG3D.MapB {
		if _, ok := el.BCMaps.NodeBC[i]; !ok {
			return fmt.Errorf("boundary node index %d (MapB position) has no BC type assigned", i)
		}
	}

	// Check that BC types are valid
	for i, bcType := range el.BCMaps.NodeBC {
		if bcType >= utils.BCUserDefined5 {
			return fmt.Errorf("invalid BC type %d at boundary node index %d", bcType, i)
		}
	}

	return nil
}
