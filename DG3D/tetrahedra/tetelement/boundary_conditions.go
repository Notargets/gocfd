package tetelement

import (
	"fmt"
	"sort"
	"strings"

	"github.com/notargets/gocfd/DG3D/mesh"
)

// Boundary condition type constants
const (
	BCInterior          = 0
	BCWall              = 1
	BCInlet             = 2
	BCOutlet            = 3
	BCFarfield          = 4
	BCSymmetry          = 5
	BCPeriodic          = 6
	BCPartitionBoundary = 255 // Special marker for partition interfaces
)

// BoundaryFaceInfo stores information about a boundary face
type BoundaryFaceInfo struct {
	ElementIndex int
	FaceIndex    int
	BCType       int
}

// IsPartitionBoundary checks if a BC type is a partition boundary
func IsPartitionBoundary(bcType int) bool {
	return bcType == BCPartitionBoundary
}

// extractBoundaryConditions extracts boundary condition types from mesh boundary elements
// and populates the BCType array for each face of each element
func (el *Element3D) extractBoundaryConditions() error {
	// Initialize BCType array - 4 faces per tet
	el.BCType = make([]int, el.K*4)

	// Default all faces to interior (0) or wall (-1)
	for i := range el.BCType {
		el.BCType[i] = 0 // 0 typically means interior face
	}

	// Create a map from boundary tag names to BC type integers
	// This mapping should be configurable based on your solver's BC types
	bcTypeMap := map[string]int{
		"wall":     1,
		"inlet":    2,
		"outlet":   3,
		"farfield": 4,
		"symmetry": 5,
		"periodic": 6,
		// Add more as needed
	}

	// If no boundary elements exist, return early
	if len(el.Mesh.BoundaryElements) == 0 {
		return nil
	}

	// For each boundary tag and its elements
	for tagName, boundaryElems := range el.Mesh.BoundaryElements {
		// Get BC type for this tag
		bcType, exists := bcTypeMap[strings.ToLower(tagName)]
		if !exists {
			// If tag name not in map, try to parse as a number or use default
			if strings.Contains(strings.ToLower(tagName), "wall") {
				bcType = 1
			} else if strings.Contains(strings.ToLower(tagName), "inlet") {
				bcType = 2
			} else if strings.Contains(strings.ToLower(tagName), "outlet") {
				bcType = 3
			} else {
				bcType = 1 // Default to wall
			}
		}

		// Process each boundary element
		for _, belem := range boundaryElems {
			// Skip if not a triangular face (tets have triangular faces)
			if belem.ElementType != mesh.Triangle {
				continue
			}

			// Find which element and face this boundary corresponds to
			if belem.ParentElement >= 0 && belem.ParentFace >= 0 {
				// If parent info is available (from some mesh formats)
				if elemIdx := el.findTetIndex(belem.ParentElement); elemIdx >= 0 {
					faceIdx := elemIdx*4 + belem.ParentFace
					if faceIdx < len(el.BCType) {
						el.BCType[faceIdx] = bcType
					}
				}
			} else {
				// Otherwise, search for matching face
				el.assignBCByFaceMatch(belem.Nodes, bcType)
			}
		}
	}

	// Also mark physical boundary faces (where EToE = -1)
	if el.EToE != nil {
		for k := 0; k < el.K; k++ {
			for f := 0; f < 4; f++ {
				if k < len(el.EToE) && f < len(el.EToE[k]) {
					if el.EToE[k][f] < 0 {
						// This is a boundary face
						bcIdx := k*4 + f
						// Only set if not already set by boundary elements
						if el.BCType[bcIdx] == 0 {
							el.BCType[bcIdx] = 1 // Default to wall
						}
					}
				}
			}
		}
	}

	return nil
}

// findTetIndex finds the local tet index from the original mesh element index
func (el *Element3D) findTetIndex(meshElemIdx int) int {
	// If we have tetIndices mapping, use it
	if len(el.tetIndices) > 0 {
		for localIdx, origIdx := range el.tetIndices {
			if origIdx == meshElemIdx {
				return localIdx
			}
		}
	}

	// Otherwise assume 1-to-1 mapping if within range
	if meshElemIdx < el.K {
		return meshElemIdx
	}

	return -1
}

// assignBCByFaceMatch finds matching faces by comparing nodes
func (el *Element3D) assignBCByFaceMatch(boundaryNodes []int, bcType int) {
	// Get face definitions for tetrahedra
	tetFaces := [][]int{
		{0, 2, 1}, // Face 0
		{0, 1, 3}, // Face 1
		{1, 2, 3}, // Face 2
		{0, 3, 2}, // Face 3
	}

	// For each tet element
	for k := 0; k < el.K; k++ {
		if k >= len(el.EToV) {
			continue
		}

		// For each face of the tet
		for f := 0; f < 4; f++ {
			// Get the nodes of this face
			faceNodes := make([]int, 3)
			for i := 0; i < 3; i++ {
				localNode := tetFaces[f][i]
				if localNode < len(el.EToV[k]) {
					faceNodes[i] = el.EToV[k][localNode]
				}
			}

			// Check if this face matches the boundary element
			if facesMatch(faceNodes, boundaryNodes) {
				bcIdx := k*4 + f
				if bcIdx < len(el.BCType) {
					el.BCType[bcIdx] = bcType
				}
				break // Move to next element
			}
		}
	}
}

// facesMatch checks if two triangular faces have the same nodes (regardless of order)
func facesMatch(face1, face2 []int) bool {
	if len(face1) != 3 || len(face2) != 3 {
		return false
	}

	// Create sorted copies
	f1 := make([]int, 3)
	f2 := make([]int, 3)
	copy(f1, face1)
	copy(f2, face2)
	sort.Ints(f1)
	sort.Ints(f2)

	// Compare sorted nodes
	return f1[0] == f2[0] && f1[1] == f2[1] && f1[2] == f2[2]
}

// extractPeriodicBoundaries processes periodic boundary conditions from the mesh
func (el *Element3D) extractPeriodicBoundaries() error {
	if len(el.Mesh.Periodics) == 0 {
		return nil
	}

	// Process each periodic boundary pair
	for _, periodic := range el.Mesh.Periodics {
		// Skip if not 2D periodic (face periodic)
		if periodic.Dimension != 2 {
			continue
		}

		// Mark all faces involved in this periodic boundary
		// The NodeMap contains slave->master vertex mappings
		for slaveNode, masterNode := range periodic.NodeMap {
			el.markPeriodicFaces(slaveNode, masterNode, 6) // 6 = periodic BC type
		}

		// If affine transformation is provided, we'd need to apply it
		// to properly map face orientations - not implemented here
		if len(periodic.AffineTransform) > 0 {
			// TODO: Apply affine transformation to handle rotated/reflected periodic faces
			// This would modify how we match face nodes in updatePeriodicMapping
		}
	}

	return nil
}

// markPeriodicFaces marks faces containing periodic nodes
func (el *Element3D) markPeriodicFaces(slaveNode, masterNode int, bcType int) {
	// Standard tet face definitions (local vertex indices)
	tetFaces := [][]int{
		{0, 2, 1}, // Face 0
		{0, 1, 3}, // Face 1
		{1, 2, 3}, // Face 2
		{0, 3, 2}, // Face 3
	}

	// Find all faces containing the slave node
	for k := 0; k < el.K; k++ {
		if k >= len(el.EToV) {
			continue
		}

		for f := 0; f < 4; f++ {
			// Get the actual vertex indices for this face
			faceContainsSlaveNode := false
			faceVertices := make([]int, 3)

			for i := 0; i < 3; i++ {
				localVertexIdx := tetFaces[f][i]
				if localVertexIdx < len(el.EToV[k]) {
					globalVertexIdx := el.EToV[k][localVertexIdx]
					faceVertices[i] = globalVertexIdx
					if globalVertexIdx == slaveNode {
						faceContainsSlaveNode = true
					}
				}
			}

			// If this face contains the slave node, mark it as periodic
			if faceContainsSlaveNode {
				bcIdx := k*4 + f
				if bcIdx < len(el.BCType) {
					el.BCType[bcIdx] = bcType
				}

				// Update periodic mapping in VmapP
				el.updatePeriodicMapping(k, f, slaveNode, masterNode)
			}
		}
	}
}

// updatePeriodicMapping updates VmapP for periodic faces
func (el *Element3D) updatePeriodicMapping(slaveElem, slaveFace int, slaveNode, masterNode int) {
	// Find the master element and face that contains the master node
	tetFaces := [][]int{
		{0, 2, 1}, // Face 0
		{0, 1, 3}, // Face 1
		{1, 2, 3}, // Face 2
		{0, 3, 2}, // Face 3
	}

	for k2 := 0; k2 < el.K; k2++ {
		if k2 >= len(el.EToV) {
			continue
		}

		for f2 := 0; f2 < 4; f2++ {
			// Skip if same element/face
			if k2 == slaveElem && f2 == slaveFace {
				continue
			}

			// Check if this face contains the master node
			faceContainsMasterNode := false
			for i := 0; i < 3; i++ {
				localVertexIdx := tetFaces[f2][i]
				if localVertexIdx < len(el.EToV[k2]) {
					if el.EToV[k2][localVertexIdx] == masterNode {
						faceContainsMasterNode = true
						break
					}
				}
			}

			if faceContainsMasterNode {
				// Update VmapP for all nodes on the slave face
				for i := 0; i < el.Nfp; i++ {
					// Index into VmapP for slave face
					slaveIdx := slaveFace*el.Nfp + i + slaveElem*el.Nfp*4
					// Index into volume nodes for master face
					masterVolumeNode := el.Fmask[f2][i] + k2*el.Np

					if slaveIdx < len(el.VmapP) {
						el.VmapP[slaveIdx] = masterVolumeNode
					}
				}

				// Also update MapP to connect face nodes
				for i := 0; i < el.Nfp; i++ {
					slaveMapIdx := slaveFace*el.Nfp + i + slaveElem*el.Nfp*4
					masterMapIdx := f2*el.Nfp + i + k2*el.Nfp*4

					if slaveMapIdx < len(el.MapP) && masterMapIdx < len(el.MapM) {
						el.MapP[slaveMapIdx] = el.MapM[masterMapIdx]
					}
				}

				return
			}
		}
	}
}

// GetBCType returns the boundary condition type for a specific element face
func (el *Element3D) GetBCType(elemIdx, faceIdx int) int {
	idx := elemIdx*4 + faceIdx
	if idx < len(el.BCType) {
		return el.BCType[idx]
	}
	return 0 // Default to interior
}

// SetBCType sets the boundary condition type for a specific element face
func (el *Element3D) SetBCType(elemIdx, faceIdx int, bcType int) error {
	idx := elemIdx*4 + faceIdx
	if idx >= len(el.BCType) {
		return fmt.Errorf("face index out of bounds: elem %d, face %d", elemIdx, faceIdx)
	}
	el.BCType[idx] = bcType
	return nil
}

// GetBoundaryFaces returns all faces with non-zero (boundary) BC types
func (el *Element3D) GetBoundaryFaces() []BoundaryFaceInfo {
	var boundaryFaces []BoundaryFaceInfo

	for k := 0; k < el.K; k++ {
		for f := 0; f < 4; f++ {
			bcType := el.GetBCType(k, f)
			if bcType != 0 { // Non-interior face
				boundaryFaces = append(boundaryFaces, BoundaryFaceInfo{
					ElementIndex: k,
					FaceIndex:    f,
					BCType:       bcType,
				})
			}
		}
	}

	return boundaryFaces
}

// PrintBoundaryStatistics prints summary of boundary conditions
func (el *Element3D) PrintBoundaryStatistics() {
	bcCounts := make(map[int]int)
	totalFaces := el.K * 4

	for i := 0; i < totalFaces; i++ {
		bcCounts[el.BCType[i]]++
	}

	fmt.Println("Boundary Condition Statistics:")
	fmt.Printf("  Total faces: %d\n", totalFaces)

	bcNames := map[int]string{
		0:   "Interior",
		1:   "Wall",
		2:   "Inlet",
		3:   "Outlet",
		4:   "Farfield",
		5:   "Symmetry",
		6:   "Periodic",
		255: "PartitionBoundary",
	}

	for bcType, count := range bcCounts {
		name, ok := bcNames[bcType]
		if !ok {
			name = fmt.Sprintf("Type %d", bcType)
		}
		fmt.Printf("  %s: %d faces (%.1f%%)\n", name, count,
			float64(count)/float64(totalFaces)*100)
	}
}

// GetFaceNodes returns the node indices for a specific face of an element
func (el *Element3D) GetFaceNodes(elemIdx, faceIdx int) []int {
	if elemIdx >= el.K || faceIdx >= 4 {
		return nil
	}

	// Standard tet face definitions
	tetFaces := [][]int{
		{0, 2, 1}, // Face 0
		{0, 1, 3}, // Face 1
		{1, 2, 3}, // Face 2
		{0, 3, 2}, // Face 3
	}

	faceNodes := make([]int, 3)
	for i := 0; i < 3; i++ {
		localNode := tetFaces[faceIdx][i]
		if localNode < len(el.EToV[elemIdx]) {
			faceNodes[i] = el.EToV[elemIdx][localNode]
		}
	}

	return faceNodes
}

// ValidateBoundaryConditions checks consistency of boundary conditions
func (el *Element3D) ValidateBoundaryConditions() error {
	// Check that BCType array is properly sized
	expectedSize := el.K * 4
	if len(el.BCType) != expectedSize {
		return fmt.Errorf("BCType array size mismatch: expected %d, got %d",
			expectedSize, len(el.BCType))
	}

	// Check that boundary faces (EToE = -1) have non-zero BC types
	inconsistencies := 0
	for k := 0; k < el.K; k++ {
		for f := 0; f < 4; f++ {
			if k < len(el.EToE) && f < len(el.EToE[k]) {
				neighborElem := el.EToE[k][f]
				bcType := el.GetBCType(k, f)

				if neighborElem < 0 && bcType == 0 {
					fmt.Printf("Warning: Face %d of element %d is boundary (EToE=-1) but BCType=0\n",
						f, k)
					inconsistencies++
				} else if neighborElem >= 0 && bcType != 0 {
					fmt.Printf("Warning: Face %d of element %d is interior (EToE=%d) but BCType=%d\n",
						f, k, neighborElem, bcType)
					inconsistencies++
				}
			}
		}
	}

	if inconsistencies > 0 {
		return fmt.Errorf("found %d boundary condition inconsistencies", inconsistencies)
	}

	return nil
}

// GetBCTypeMap returns the default mapping of BC names to type integers
func GetBCTypeMap() map[string]int {
	return map[string]int{
		"interior":           0,
		"wall":               1,
		"inlet":              2,
		"outlet":             3,
		"farfield":           4,
		"symmetry":           5,
		"periodic":           6,
		"partition_boundary": 255,
	}
}

// GetBCTypeName returns the name for a given BC type integer
func GetBCTypeName(bcType int) string {
	names := map[int]string{
		0:   "Interior",
		1:   "Wall",
		2:   "Inlet",
		3:   "Outlet",
		4:   "Farfield",
		5:   "Symmetry",
		6:   "Periodic",
		255: "PartitionBoundary",
	}

	if name, ok := names[bcType]; ok {
		return name
	}
	return fmt.Sprintf("Type%d", bcType)
}
