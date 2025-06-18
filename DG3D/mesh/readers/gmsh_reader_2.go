package readers

import (
	"bufio"
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/utils"
	"os"
	"strconv"
	"strings"
)

// ReadGmsh22 reads a Gmsh MSH file format version 2.2
func ReadGmsh22(filename string) (*mesh.Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	msh := mesh.NewMesh()

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		switch line {
		case "$MeshFormat":
			if err := readMeshFormat22(scanner, msh); err != nil {
				return nil, err
			}

		case "$PhysicalNames":
			if err := readPhysicalNames(scanner, msh); err != nil {
				return nil, err
			}

		case "$Nodes":
			if err := readNodes22(scanner, msh); err != nil {
				return nil, err
			}

		case "$Elements":
			if err := readElements22(scanner, msh); err != nil {
				return nil, err
			}

		case "$Periodic":
			if err := readPeriodic22(scanner, msh); err != nil {
				return nil, err
			}

		case "$NodeData", "$ElementData", "$ElementNodeData":
			// Skip data sections
			endMarker := "$End" + line[1:]
			for scanner.Scan() {
				if strings.TrimSpace(scanner.Text()) == endMarker {
					break
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner error: %v", err)
	}

	msh.BuildConnectivity()
	return msh, nil
}

// readMeshFormat22 reads the MeshFormat section
func readMeshFormat22(scanner *bufio.Scanner, mesh *mesh.Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in MeshFormat")
	}

	parts := strings.Fields(scanner.Text())
	if len(parts) < 3 {
		return fmt.Errorf("invalid MeshFormat line")
	}

	mesh.FormatVersion = parts[0]
	fileType, _ := strconv.Atoi(parts[1])
	mesh.IsBinary = fileType == 1
	mesh.DataSize, _ = strconv.Atoi(parts[2])

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndMeshFormat" {
			break
		}
	}

	return nil
}

// readPhysicalNames reads physical group names (common to v2.2 and v4)
func readPhysicalNames(scanner *bufio.Scanner, msh *mesh.Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in PhysicalNames")
	}

	numNames, _ := strconv.Atoi(strings.TrimSpace(scanner.Text()))

	for i := 0; i < numNames; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading physical names")
		}

		parts := strings.Fields(scanner.Text())
		if len(parts) >= 3 {
			dimension, _ := strconv.Atoi(parts[0])
			tag, _ := strconv.Atoi(parts[1])
			name := strings.Trim(parts[2], "\"")

			// Join remaining parts if name contains spaces
			for j := 3; j < len(parts); j++ {
				name += " " + strings.Trim(parts[j], "\"")
			}

			group := &mesh.ElementGroup{
				Dimension: dimension,
				Tag:       tag,
				Name:      name,
				Elements:  []int{},
			}
			msh.ElementGroups[tag] = group
		}
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndPhysicalNames" {
			break
		}
	}

	return nil
}

// readNodes22 reads nodes in v2.2 format
func readNodes22(scanner *bufio.Scanner, mesh *mesh.Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Nodes")
	}

	numNodes, _ := strconv.Atoi(strings.TrimSpace(scanner.Text()))
	mesh.Vertices = make([][]float64, 0, numNodes)

	for i := 0; i < numNodes; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading nodes")
		}

		parts := strings.Fields(scanner.Text())
		if len(parts) < 4 {
			return fmt.Errorf("invalid node line: %s", scanner.Text())
		}

		nodeID, _ := strconv.Atoi(parts[0])
		x, _ := strconv.ParseFloat(parts[1], 64)
		y, _ := strconv.ParseFloat(parts[2], 64)
		z, _ := strconv.ParseFloat(parts[3], 64)

		mesh.AddNode(nodeID, []float64{x, y, z})
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndNodes" {
			break
		}
	}

	return nil
}

// readElements22 reads elements in v2.2 format
func readElements22(scanner *bufio.Scanner, msh *mesh.Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Elements")
	}

	numElements, _ := strconv.Atoi(strings.TrimSpace(scanner.Text()))

	for i := 0; i < numElements; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading elements")
		}

		parts := strings.Fields(scanner.Text())
		if len(parts) < 5 {
			return fmt.Errorf("invalid element line")
		}

		elemID, _ := strconv.Atoi(parts[0])
		elemType, _ := strconv.Atoi(parts[1])
		numTags, _ := strconv.Atoi(parts[2])

		if len(parts) < 3+numTags {
			return fmt.Errorf("invalid element tags")
		}

		// Read tags
		tags := make([]int, numTags)
		for j := 0; j < numTags; j++ {
			tags[j], _ = strconv.Atoi(parts[3+j])
		}

		// Map element type
		etype, ok := gmshElementType22[elemType]
		if !ok {
			// Skip unknown element types
			continue
		}

		// Check if it's a boundary element
		if etype.GetDimension() < msh.GetMeshDimension() {
			// This is a boundary element
			handleBoundaryElement22(msh, etype, tags, parts[3+numTags:])
			continue
		}

		// Read nodes
		expectedNodes := etype.GetNumNodes()
		nodeStart := 3 + numTags
		if len(parts) < nodeStart+expectedNodes {
			return fmt.Errorf("element %d: expected %d nodes, got %d",
				elemID, expectedNodes, len(parts)-nodeStart)
		}

		nodeIDs := make([]int, expectedNodes)
		for j := 0; j < expectedNodes; j++ {
			nodeIDs[j], _ = strconv.Atoi(parts[nodeStart+j])
		}

		if err := msh.AddElement(elemID, etype, tags, nodeIDs); err != nil {
			return err
		}
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndElements" {
			break
		}
	}

	return nil
}

// handleBoundaryElement22 processes boundary elements
func handleBoundaryElement22(msh *mesh.Mesh, etype utils.ElementType, tags []int, nodeStrs []string) {
	// Read nodes
	numNodes := etype.GetNumNodes()
	if len(nodeStrs) < numNodes {
		return // Skip if not enough nodes
	}

	nodes := make([]int, numNodes)
	for i := 0; i < numNodes; i++ {
		nodeID, _ := strconv.Atoi(nodeStrs[i])
		idx, ok := msh.GetNodeIndex(nodeID)
		if !ok {
			return // Skip if node not found
		}
		nodes[i] = idx
	}

	// Get physical tag if present
	var physicalTag int
	if len(tags) > 0 {
		physicalTag = tags[0]
	}

	// Find or create boundary tag name
	var tagName string
	if group, ok := msh.ElementGroups[physicalTag]; ok {
		tagName = group.Name
	} else {
		tagName = fmt.Sprintf("boundary_%d", physicalTag)
	}

	// Create boundary element
	belem := mesh.BoundaryElement{
		ElementType:   etype,
		Nodes:         nodes,
		ParentElement: -1, // Not tracked in v2.2
		ParentFace:    -1, // Not tracked in v2.2
	}

	msh.AddBoundaryElement(tagName, belem)
}

// readPeriodic22 reads periodic boundary conditions in v2.2 format
func readPeriodic22(scanner *bufio.Scanner, msh *mesh.Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Periodic")
	}

	numPeriodic, _ := strconv.Atoi(strings.TrimSpace(scanner.Text()))

	for i := 0; i < numPeriodic; i++ {
		// Read periodic header
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading periodic %d", i)
		}

		parts := strings.Fields(scanner.Text())
		if len(parts) < 3 {
			return fmt.Errorf("invalid periodic header")
		}

		dimension, _ := strconv.Atoi(parts[0])
		slaveTag, _ := strconv.Atoi(parts[1])
		masterTag, _ := strconv.Atoi(parts[2])

		// Check if next line is affine transformation
		var affineTransform []float64
		if scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())

			// Check if line starts with "Affine"
			if strings.HasPrefix(line, "Affine") {
				// Parse affine values
				fields := strings.Fields(line)
				numValues := (dimension + 1) * (dimension + 1)
				affineTransform = make([]float64, numValues)

				// Read values from the Affine line (skip "Affine" keyword)
				affineIdx := 0
				for j := 1; j < len(fields) && affineIdx < numValues; j++ {
					affineTransform[affineIdx], _ = strconv.ParseFloat(fields[j], 64)
					affineIdx++
				}

				// Read remaining affine values if needed
				for affineIdx < numValues {
					if !scanner.Scan() {
						return fmt.Errorf("unexpected EOF reading affine transform")
					}
					affineFields := strings.Fields(scanner.Text())
					for _, field := range affineFields {
						if affineIdx < numValues {
							affineTransform[affineIdx], _ = strconv.ParseFloat(field, 64)
							affineIdx++
						}
					}
				}

				// Now read the number of node correspondences
				if !scanner.Scan() {
					return fmt.Errorf("unexpected EOF reading node count")
				}
				line = strings.TrimSpace(scanner.Text())
			}

			// Parse number of node correspondences
			numNodes, _ := strconv.Atoi(line)

			// Read node correspondences
			nodeMap := make(map[int]int)
			for j := 0; j < numNodes; j++ {
				if !scanner.Scan() {
					return fmt.Errorf("unexpected EOF reading node correspondence")
				}

				fields := strings.Fields(scanner.Text())
				if len(fields) >= 2 {
					slaveNode, _ := strconv.Atoi(fields[0])
					masterNode, _ := strconv.Atoi(fields[1])
					nodeMap[slaveNode] = masterNode
				}
			}

			// Create periodic structure
			periodic := mesh.Periodic{
				Dimension:       dimension,
				SlaveTag:        slaveTag,
				MasterTag:       masterTag,
				NodeMap:         nodeMap,
				AffineTransform: affineTransform,
			}

			msh.Periodics = append(msh.Periodics, periodic)
		}
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndPeriodic" {
			break
		}
	}

	return nil
}

// gmshElementType22 maps Gmsh v2.2 element type numbers to our ElementType
var gmshElementType22 = map[int]utils.ElementType{
	1:  utils.Line,       // 2-node line
	2:  utils.Triangle,   // 3-node triangle
	3:  utils.Quad,       // 4-node quadrangle
	4:  utils.Tet,        // 4-node tetrahedron
	5:  utils.Hex,        // 8-node hexahedron
	6:  utils.Prism,      // 6-node prism
	7:  utils.Pyramid,    // 5-node pyramid
	8:  utils.Line3,      // 3-node line
	9:  utils.Triangle6,  // 6-node triangle
	10: utils.Quad9,      // 9-node quadrangle
	11: utils.Tet10,      // 10-node tetrahedron
	15: utils.Point,      // 1-node point
	16: utils.Quad8,      // 8-node quadrangle
	17: utils.Hex20,      // 20-node hexahedron
	18: utils.Prism15,    // 15-node prism
	19: utils.Pyramid13,  // 13-node pyramid
	20: utils.Triangle9,  // 9-node triangle
	21: utils.Triangle10, // 10-node triangle
	22: utils.Quad9,      // 9-node quadrangle (duplicate)
	26: utils.Line,       // 2-node line (duplicate)
	27: utils.Line3,      // 3-node line (duplicate)
	28: utils.Line,       // 2-node line (duplicate)
	29: utils.Tet10,      // 10-node tetrahedron (duplicate)
	30: utils.Prism18,    // 18-node prism
	31: utils.Pyramid14,  // 14-node pyramid
	92: utils.Hex27,      // 27-node hexahedron
}
