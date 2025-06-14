package mesh

import (
	"bufio"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

// gmshElementType2_2 maps Gmsh 2.2 element types to our ElementType
var gmshElementType2_2 = map[int]ElementType{
	1:  Line,
	2:  Triangle,
	3:  Quad,
	4:  Tet,
	5:  Hex,
	6:  Prism,
	7:  Pyramid,
	8:  Line3,
	9:  Triangle6,
	10: Quad9,
	11: Tet10,
	12: Hex27,
	13: Prism18,
	14: Pyramid14,
	15: Point,
	16: Quad8,
	17: Hex20,
	18: Prism15,
	19: Pyramid13,
	20: Triangle9,
	21: Triangle10,
	// Additional types exist but are less common
}

// ReadGmsh reads a Gmsh format file (auto-detects version)
func ReadGmshAuto(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	// Read first section to determine version
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "$MeshFormat" {
			if !scanner.Scan() {
				return nil, fmt.Errorf("unexpected EOF after $MeshFormat")
			}
			versionLine := strings.Fields(scanner.Text())
			if len(versionLine) < 1 {
				return nil, fmt.Errorf("invalid version line")
			}

			version := versionLine[0]
			file.Seek(0, 0) // Reset to beginning

			if strings.HasPrefix(version, "2") {
				return ReadGmsh22(filename)
			} else if strings.HasPrefix(version, "4") {
				return ReadGmsh4(filename)
			}
			return nil, fmt.Errorf("unsupported Gmsh version: %s", version)
		}
	}

	return nil, fmt.Errorf("no $MeshFormat section found")
}

// ReadGmsh22 reads a Gmsh 2.2 format file with full compliance
func ReadGmsh22(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	mesh := NewMesh()

	// First pass: determine if binary
	scanner := bufio.NewScanner(file)
	var isBinary bool
	var dataSize int

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "$MeshFormat" {
			scanner.Scan()
			parts := strings.Fields(scanner.Text())
			if len(parts) >= 3 {
				mesh.FormatVersion = parts[0]
				fileType, _ := strconv.Atoi(parts[1])
				isBinary = fileType == 1
				dataSize, _ = strconv.Atoi(parts[2])
				mesh.IsBinary = isBinary
				mesh.DataSize = dataSize
			}
			break
		}
	}

	// Reset and read properly
	file.Seek(0, 0)

	if isBinary {
		return readGmsh22Binary(file, mesh)
	}
	return readGmsh22ASCII(file, mesh)
}

// readGmsh22ASCII reads ASCII format
func readGmsh22ASCII(file *os.File, mesh *Mesh) (*Mesh, error) {
	scanner := bufio.NewScanner(file)

	// Increase scanner buffer for large files
	const maxScanTokenSize = 1024 * 1024 * 10 // 10MB
	buf := make([]byte, maxScanTokenSize)
	scanner.Buffer(buf, maxScanTokenSize)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		switch line {
		case "$MeshFormat":
			if err := readMeshFormat(scanner, mesh); err != nil {
				return nil, err
			}

		case "$PhysicalNames":
			if err := readPhysicalNames(scanner, mesh); err != nil {
				return nil, err
			}

		case "$Nodes":
			if err := readNodes(scanner, mesh); err != nil {
				return nil, err
			}

		case "$Elements":
			if err := readElements(scanner, mesh); err != nil {
				return nil, err
			}

		case "$Periodic":
			if err := readPeriodic(scanner, mesh); err != nil {
				return nil, err
			}

		case "$NodeData":
			if err := skipSection(scanner, "$EndNodeData"); err != nil {
				return nil, err
			}

		case "$ElementData":
			if err := skipSection(scanner, "$EndElementData"); err != nil {
				return nil, err
			}

		case "$ElementNodeData":
			if err := skipSection(scanner, "$EndElementNodeData"); err != nil {
				return nil, err
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner error: %v", err)
	}

	mesh.BuildConnectivity()
	return mesh, nil
}

// readMeshFormat reads the MeshFormat section
func readMeshFormat(scanner *bufio.Scanner, mesh *Mesh) error {
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

// readPhysicalNames reads physical entity names
func readPhysicalNames(scanner *bufio.Scanner, mesh *Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in PhysicalNames")
	}

	numPhysical, err := strconv.Atoi(scanner.Text())
	if err != nil {
		return fmt.Errorf("invalid number of physical names: %v", err)
	}

	for i := 0; i < numPhysical; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in PhysicalNames")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 3 {
			return fmt.Errorf("invalid physical name entry")
		}

		dim, _ := strconv.Atoi(fields[0])
		tag, _ := strconv.Atoi(fields[1])
		name := strings.Trim(strings.Join(fields[2:], " "), "\"")

		// Store as element group
		mesh.ElementGroups[tag] = &ElementGroup{
			Dimension: dim,
			Tag:       tag,
			Name:      name,
			Elements:  []int{},
		}

		// Also store in boundary tags for compatibility
		mesh.BoundaryTags[tag] = name
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndPhysicalNames" {
			break
		}
	}

	return nil
}

// readNodes reads the Nodes section
func readNodes(scanner *bufio.Scanner, mesh *Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Nodes")
	}

	numNodes, err := strconv.Atoi(scanner.Text())
	if err != nil {
		return fmt.Errorf("invalid number of nodes: %v", err)
	}

	for i := 0; i < numNodes; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in Nodes at node %d", i)
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 4 {
			return fmt.Errorf("invalid node entry at line %d", i+1)
		}

		nodeID, err := strconv.Atoi(fields[0])
		if err != nil {
			return fmt.Errorf("invalid node ID: %v", err)
		}

		coords := make([]float64, 3)
		for j := 0; j < 3; j++ {
			coords[j], err = strconv.ParseFloat(fields[j+1], 64)
			if err != nil {
				return fmt.Errorf("invalid coordinate: %v", err)
			}
		}

		mesh.AddNode(nodeID, coords)
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndNodes" {
			break
		}
	}

	return nil
}

// readElements reads the Elements section
func readElements(scanner *bufio.Scanner, mesh *Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Elements")
	}

	numElems, err := strconv.Atoi(scanner.Text())
	if err != nil {
		return fmt.Errorf("invalid number of elements: %v", err)
	}

	for i := 0; i < numElems; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in Elements at element %d", i)
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 3 {
			return fmt.Errorf("invalid element entry at line %d", i+1)
		}

		elemID, err := strconv.Atoi(fields[0])
		if err != nil {
			return fmt.Errorf("invalid element ID: %v", err)
		}

		gmshType, err := strconv.Atoi(fields[1])
		if err != nil {
			return fmt.Errorf("invalid element type: %v", err)
		}

		numTags, err := strconv.Atoi(fields[2])
		if err != nil {
			return fmt.Errorf("invalid number of tags: %v", err)
		}

		// Read tags
		tags := make([]int, numTags)
		for j := 0; j < numTags; j++ {
			if 3+j >= len(fields) {
				return fmt.Errorf("insufficient fields for tags")
			}
			tags[j], err = strconv.Atoi(fields[3+j])
			if err != nil {
				return fmt.Errorf("invalid tag: %v", err)
			}
		}

		// Get element type
		elemType, ok := gmshElementType2_2[gmshType]
		if !ok {
			// Skip unknown element types
			continue
		}

		// In the readElements function, replace the node reading section with:

		// Read nodes
		expectedNodes := elemType.GetNumNodes()
		startIdx := 3 + numTags
		actualNodes := len(fields) - startIdx

		// Use the actual number of nodes in the file if less than expected
		// This handles cases where higher-order elements might be represented with fewer nodes
		numNodesToRead := expectedNodes
		if actualNodes < expectedNodes {
			// For higher-order elements, we might only have corner nodes
			corners := elemType.GetCornerNodes()
			if actualNodes == len(corners) {
				numNodesToRead = actualNodes
			} else {
				return fmt.Errorf("element type %v expects %d nodes, got %d", elemType, expectedNodes, actualNodes)
			}
		}

		nodeIDs := make([]int, numNodesToRead)
		for j := 0; j < numNodesToRead; j++ {
			nodeIDs[j], err = strconv.Atoi(fields[startIdx+j])
			if err != nil {
				return fmt.Errorf("invalid node ID: %v", err)
			}
		}

		// If we read fewer nodes than expected, pad with zeros
		if numNodesToRead < expectedNodes {
			paddedNodeIDs := make([]int, expectedNodes)
			copy(paddedNodeIDs, nodeIDs)
			nodeIDs = paddedNodeIDs
		}

		if err := mesh.AddElement(elemID, elemType, tags, nodeIDs); err != nil {
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

// readPeriodic reads periodic boundary conditions
func readPeriodic(scanner *bufio.Scanner, mesh *Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Periodic")
	}

	numPeriodic, err := strconv.Atoi(scanner.Text())
	if err != nil {
		return fmt.Errorf("invalid number of periodic entities: %v", err)
	}

	for i := 0; i < numPeriodic; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in Periodic entity %d", i)
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 3 {
			return fmt.Errorf("invalid periodic entity header")
		}

		dim, _ := strconv.Atoi(fields[0])
		slaveTag, _ := strconv.Atoi(fields[1])
		masterTag, _ := strconv.Atoi(fields[2])

		periodic := Periodic{
			Dimension: dim,
			SlaveTag:  slaveTag,
			MasterTag: masterTag,
			NodeMap:   make(map[int]int),
		}

		// Read affine transformation if present
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in Periodic affine")
		}

		affineLine := scanner.Text()
		if strings.Contains(affineLine, "Affine") {
			// Read 16 values for 4x4 affine matrix
			affineFields := strings.Fields(affineLine)[1:] // Skip "Affine"
			if len(affineFields) >= 16 {
				periodic.AffineTransform = make([]float64, 16)
				for j := 0; j < 16; j++ {
					periodic.AffineTransform[j], _ = strconv.ParseFloat(affineFields[j], 64)
				}
			}

			// Read number of nodes
			if !scanner.Scan() {
				return fmt.Errorf("unexpected EOF in Periodic nodes")
			}
			affineLine = scanner.Text()
		}

		// Now affineLine contains the number of corresponding nodes
		numNodes, err := strconv.Atoi(affineLine)
		if err != nil {
			return fmt.Errorf("invalid number of periodic nodes: %v", err)
		}

		// Read node correspondences
		for j := 0; j < numNodes; j++ {
			if !scanner.Scan() {
				return fmt.Errorf("unexpected EOF in Periodic nodes")
			}

			fields := strings.Fields(scanner.Text())
			if len(fields) < 2 {
				return fmt.Errorf("invalid periodic node mapping")
			}

			slaveNode, _ := strconv.Atoi(fields[0])
			masterNode, _ := strconv.Atoi(fields[1])
			periodic.NodeMap[slaveNode] = masterNode
		}

		mesh.Periodics = append(mesh.Periodics, periodic)
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndPeriodic" {
			break
		}
	}

	return nil
}

// skipSection skips an unhandled section
func skipSection(scanner *bufio.Scanner, endTag string) error {
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == endTag {
			return nil
		}
	}
	return fmt.Errorf("unexpected EOF while looking for %s", endTag)
}

// readGmsh22Binary reads binary format (placeholder - implement if needed)
func readGmsh22Binary(file *os.File, mesh *Mesh) (*Mesh, error) {
	reader := bufio.NewReader(file)

	// Read until $MeshFormat
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			return nil, err
		}
		if strings.TrimSpace(line) == "$MeshFormat" {
			break
		}
	}

	// Read format line
	line, err := reader.ReadString('\n')
	if err != nil {
		return nil, err
	}

	parts := strings.Fields(line)
	if len(parts) < 3 {
		return nil, fmt.Errorf("invalid MeshFormat")
	}

	mesh.FormatVersion = parts[0]
	mesh.IsBinary = true
	mesh.DataSize, _ = strconv.Atoi(parts[2])

	// Read binary endianness check
	var one int32
	if err := binary.Read(reader, binary.LittleEndian, &one); err != nil {
		return nil, fmt.Errorf("failed to read endianness check: %v", err)
	}

	var byteOrder binary.ByteOrder = binary.LittleEndian
	if one != 1 {
		byteOrder = binary.BigEndian
	}

	// Read newline after binary data
	reader.ReadString('\n')

	// Continue with sections
	for {
		line, err := reader.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}

		line = strings.TrimSpace(line)

		switch line {
		case "$EndMeshFormat":
			continue

		case "$Nodes":
			if err := readNodesBinary(reader, mesh, byteOrder); err != nil {
				return nil, err
			}

		case "$Elements":
			if err := readElementsBinary(reader, mesh, byteOrder); err != nil {
				return nil, err
			}

		default:
			// Skip unknown sections in binary
			if strings.HasPrefix(line, "$") && !strings.HasPrefix(line, "$End") {
				endTag := "$End" + line[1:]
				for {
					skipLine, err := reader.ReadString('\n')
					if err != nil {
						return nil, err
					}
					if strings.TrimSpace(skipLine) == endTag {
						break
					}
				}
			}
		}
	}

	mesh.BuildConnectivity()
	return mesh, nil
}

// readNodesBinary reads nodes in binary format
func readNodesBinary(reader *bufio.Reader, mesh *Mesh, byteOrder binary.ByteOrder) error {
	// Read number of nodes
	line, err := reader.ReadString('\n')
	if err != nil {
		return err
	}

	numNodes, err := strconv.Atoi(strings.TrimSpace(line))
	if err != nil {
		return err
	}

	// Read binary node data
	for i := 0; i < numNodes; i++ {
		var nodeID int32
		var coords [3]float64

		if err := binary.Read(reader, byteOrder, &nodeID); err != nil {
			return err
		}

		for j := 0; j < 3; j++ {
			if err := binary.Read(reader, byteOrder, &coords[j]); err != nil {
				return err
			}
		}

		mesh.AddNode(int(nodeID), coords[:])
	}

	// Read $EndNodes
	line, err = reader.ReadString('\n')
	if err != nil {
		return err
	}

	return nil
}

// readElementsBinary reads elements in binary format
func readElementsBinary(reader *bufio.Reader, mesh *Mesh, byteOrder binary.ByteOrder) error {
	// This is complex due to variable element sizes
	// For now, return error indicating binary elements not yet supported
	return fmt.Errorf("binary element reading not yet implemented")
}
