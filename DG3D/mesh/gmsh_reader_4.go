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

// gmshElementType4 maps Gmsh 4.x element types to our ElementType
// Note: Gmsh 4.x uses the same element type numbers as 2.2
var gmshElementType4 = gmshElementType2_2

// EntityInfo stores information about geometric entities
type EntityInfo struct {
	Tag              int
	Dimension        int
	MinX, MinY, MinZ float64
	MaxX, MaxY, MaxZ float64
	PhysicalTags     []int
	NumNodes         int
	NodeTags         []int
}

// ReadGmsh4 reads a Gmsh 4.x format file
func ReadGmsh4(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	mesh := NewMesh()

	// Parse format to determine if binary
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
		return readGmsh4Binary(file, mesh)
	}
	return readGmsh4ASCII(file, mesh)
}

// readGmsh4ASCII reads ASCII format for version 4.x
func readGmsh4ASCII(file *os.File, mesh *Mesh) (*Mesh, error) {
	scanner := bufio.NewScanner(file)

	// Increase scanner buffer for large files
	const maxScanTokenSize = 1024 * 1024 * 10 // 10MB
	buf := make([]byte, maxScanTokenSize)
	scanner.Buffer(buf, maxScanTokenSize)

	// Store entities for reference
	entities := make(map[string]map[int]*EntityInfo) // dim -> tag -> info
	entities["0"] = make(map[int]*EntityInfo)        // Points
	entities["1"] = make(map[int]*EntityInfo)        // Curves
	entities["2"] = make(map[int]*EntityInfo)        // Surfaces
	entities["3"] = make(map[int]*EntityInfo)        // Volumes

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		switch line {
		case "$MeshFormat":
			if err := readMeshFormat4(scanner, mesh); err != nil {
				return nil, err
			}

		case "$PhysicalNames":
			if err := readPhysicalNames(scanner, mesh); err != nil {
				return nil, err
			}

		case "$Entities":
			if err := readEntities4(scanner, entities); err != nil {
				return nil, err
			}

		case "$PartitionedEntities":
			if err := readPartitionedEntities4(scanner, entities); err != nil {
				return nil, err
			}

		case "$Nodes":
			if err := readNodes4(scanner, mesh, entities); err != nil {
				return nil, err
			}

		case "$Elements":
			if err := readElements4(scanner, mesh, entities); err != nil {
				return nil, err
			}

		case "$Periodic":
			if err := readPeriodic4(scanner, mesh); err != nil {
				return nil, err
			}

		case "$GhostElements":
			if err := skipSection(scanner, "$EndGhostElements"); err != nil {
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

// readMeshFormat4 reads the MeshFormat section for v4
func readMeshFormat4(scanner *bufio.Scanner, mesh *Mesh) error {
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

// readEntities4 reads the Entities section (new in v4)
func readEntities4(scanner *bufio.Scanner, entities map[string]map[int]*EntityInfo) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Entities")
	}

	// Read counts: numPoints numCurves numSurfaces numVolumes
	counts := strings.Fields(scanner.Text())
	if len(counts) < 4 {
		return fmt.Errorf("invalid entity counts")
	}

	numPoints, _ := strconv.Atoi(counts[0])
	numCurves, _ := strconv.Atoi(counts[1])
	numSurfaces, _ := strconv.Atoi(counts[2])
	numVolumes, _ := strconv.Atoi(counts[3])

	// Read points (0D)
	for i := 0; i < numPoints; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in point entities")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 5 {
			return fmt.Errorf("invalid point entity")
		}

		tag, _ := strconv.Atoi(fields[0])
		x, _ := strconv.ParseFloat(fields[1], 64)
		y, _ := strconv.ParseFloat(fields[2], 64)
		z, _ := strconv.ParseFloat(fields[3], 64)
		numPhysTags, _ := strconv.Atoi(fields[4])

		entity := &EntityInfo{
			Tag:       tag,
			Dimension: 0,
			MinX:      x, MinY: y, MinZ: z,
			MaxX: x, MaxY: y, MaxZ: z,
		}

		if numPhysTags > 0 && len(fields) >= 5+numPhysTags {
			entity.PhysicalTags = make([]int, numPhysTags)
			for j := 0; j < numPhysTags; j++ {
				entity.PhysicalTags[j], _ = strconv.Atoi(fields[5+j])
			}
		}

		entities["0"][tag] = entity
	}

	// Read curves (1D)
	for i := 0; i < numCurves; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in curve entities")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 8 {
			return fmt.Errorf("invalid curve entity")
		}

		tag, _ := strconv.Atoi(fields[0])
		minX, _ := strconv.ParseFloat(fields[1], 64)
		minY, _ := strconv.ParseFloat(fields[2], 64)
		minZ, _ := strconv.ParseFloat(fields[3], 64)
		maxX, _ := strconv.ParseFloat(fields[4], 64)
		maxY, _ := strconv.ParseFloat(fields[5], 64)
		maxZ, _ := strconv.ParseFloat(fields[6], 64)
		numPhysTags, _ := strconv.Atoi(fields[7])

		entity := &EntityInfo{
			Tag:       tag,
			Dimension: 1,
			MinX:      minX, MinY: minY, MinZ: minZ,
			MaxX: maxX, MaxY: maxY, MaxZ: maxZ,
		}

		offset := 8
		if numPhysTags > 0 && len(fields) >= offset+numPhysTags {
			entity.PhysicalTags = make([]int, numPhysTags)
			for j := 0; j < numPhysTags; j++ {
				entity.PhysicalTags[j], _ = strconv.Atoi(fields[offset+j])
			}
			offset += numPhysTags
		}

		// Skip bounding points for now
		entities["1"][tag] = entity
	}

	// Read surfaces (2D)
	for i := 0; i < numSurfaces; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in surface entities")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 8 {
			return fmt.Errorf("invalid surface entity")
		}

		tag, _ := strconv.Atoi(fields[0])
		minX, _ := strconv.ParseFloat(fields[1], 64)
		minY, _ := strconv.ParseFloat(fields[2], 64)
		minZ, _ := strconv.ParseFloat(fields[3], 64)
		maxX, _ := strconv.ParseFloat(fields[4], 64)
		maxY, _ := strconv.ParseFloat(fields[5], 64)
		maxZ, _ := strconv.ParseFloat(fields[6], 64)
		numPhysTags, _ := strconv.Atoi(fields[7])

		entity := &EntityInfo{
			Tag:       tag,
			Dimension: 2,
			MinX:      minX, MinY: minY, MinZ: minZ,
			MaxX: maxX, MaxY: maxY, MaxZ: maxZ,
		}

		offset := 8
		if numPhysTags > 0 && len(fields) >= offset+numPhysTags {
			entity.PhysicalTags = make([]int, numPhysTags)
			for j := 0; j < numPhysTags; j++ {
				entity.PhysicalTags[j], _ = strconv.Atoi(fields[offset+j])
			}
		}

		// Skip bounding curves for now
		entities["2"][tag] = entity
	}

	// Read volumes (3D)
	for i := 0; i < numVolumes; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in volume entities")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 8 {
			return fmt.Errorf("invalid volume entity")
		}

		tag, _ := strconv.Atoi(fields[0])
		minX, _ := strconv.ParseFloat(fields[1], 64)
		minY, _ := strconv.ParseFloat(fields[2], 64)
		minZ, _ := strconv.ParseFloat(fields[3], 64)
		maxX, _ := strconv.ParseFloat(fields[4], 64)
		maxY, _ := strconv.ParseFloat(fields[5], 64)
		maxZ, _ := strconv.ParseFloat(fields[6], 64)
		numPhysTags, _ := strconv.Atoi(fields[7])

		entity := &EntityInfo{
			Tag:       tag,
			Dimension: 3,
			MinX:      minX, MinY: minY, MinZ: minZ,
			MaxX: maxX, MaxY: maxY, MaxZ: maxZ,
		}

		offset := 8
		if numPhysTags > 0 && len(fields) >= offset+numPhysTags {
			entity.PhysicalTags = make([]int, numPhysTags)
			for j := 0; j < numPhysTags; j++ {
				entity.PhysicalTags[j], _ = strconv.Atoi(fields[offset+j])
			}
		}

		// Skip bounding surfaces for now
		entities["3"][tag] = entity
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndEntities" {
			break
		}
	}

	return nil
}

// readPartitionedEntities4 reads partitioned entity information (for parallel meshes)
func readPartitionedEntities4(scanner *bufio.Scanner, entities map[string]map[int]*EntityInfo) error {
	// For now, skip this section as we handle partitioning separately
	return skipSection(scanner, "$EndPartitionedEntities")
}

// readNodes4 reads nodes in v4 format
func readNodes4(scanner *bufio.Scanner, mesh *Mesh, entities map[string]map[int]*EntityInfo) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Nodes")
	}

	// Format: numEntityBlocks numNodes minNodeTag maxNodeTag
	header := strings.Fields(scanner.Text())
	if len(header) < 4 {
		return fmt.Errorf("invalid Nodes header")
	}

	numEntityBlocks, _ := strconv.Atoi(header[0])
	totalNodes, _ := strconv.Atoi(header[1])
	// minNodeTag, _ := strconv.Atoi(header[2])
	// maxNodeTag, _ := strconv.Atoi(header[3])

	// Pre-allocate
	mesh.Vertices = make([][]float64, 0, totalNodes)

	// Read entity blocks
	for i := 0; i < numEntityBlocks; i++ {
		// Read entity info: entityDim entityTag parametric numNodes
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in node entity block %d", i)
		}

		blockHeader := strings.Fields(scanner.Text())
		if len(blockHeader) < 4 {
			return fmt.Errorf("invalid node block header")
		}

		// entityDim, _ := strconv.Atoi(blockHeader[0])
		// entityTag, _ := strconv.Atoi(blockHeader[1])
		parametric, _ := strconv.Atoi(blockHeader[2])
		numNodesInBlock, _ := strconv.Atoi(blockHeader[3])

		// Read node tags
		nodeTags := make([]int, numNodesInBlock)
		for j := 0; j < numNodesInBlock; j++ {
			if !scanner.Scan() {
				return fmt.Errorf("unexpected EOF reading node tags")
			}
			nodeTags[j], _ = strconv.Atoi(scanner.Text())
		}

		// Read coordinates
		for j := 0; j < numNodesInBlock; j++ {
			if !scanner.Scan() {
				return fmt.Errorf("unexpected EOF reading node coordinates")
			}

			fields := strings.Fields(scanner.Text())
			expectedFields := 3
			if parametric > 0 {
				expectedFields += parametric
			}

			if len(fields) < expectedFields {
				return fmt.Errorf("invalid node coordinate line")
			}

			coords := make([]float64, 3)
			for k := 0; k < 3; k++ {
				coords[k], _ = strconv.ParseFloat(fields[k], 64)
			}

			mesh.AddNode(nodeTags[j], coords)
		}
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndNodes" {
			break
		}
	}

	return nil
}

// readElements4 reads elements in v4 format
func readElements4(scanner *bufio.Scanner, mesh *Mesh, entities map[string]map[int]*EntityInfo) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Elements")
	}

	// Format: numEntityBlocks numElements minElementTag maxElementTag
	header := strings.Fields(scanner.Text())
	if len(header) < 4 {
		return fmt.Errorf("invalid Elements header")
	}

	numEntityBlocks, _ := strconv.Atoi(header[0])
	// totalElements, _ := strconv.Atoi(header[1])
	// minElemTag, _ := strconv.Atoi(header[2])
	// maxElemTag, _ := strconv.Atoi(header[3])

	// Read entity blocks
	for i := 0; i < numEntityBlocks; i++ {
		// Read entity info: entityDim entityTag elementType numElements
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in element entity block %d", i)
		}

		blockHeader := strings.Fields(scanner.Text())
		if len(blockHeader) < 4 {
			return fmt.Errorf("invalid element block header")
		}

		entityDim, _ := strconv.Atoi(blockHeader[0])
		entityTag, _ := strconv.Atoi(blockHeader[1])
		gmshType, _ := strconv.Atoi(blockHeader[2])
		numElemsInBlock, _ := strconv.Atoi(blockHeader[3])

		// Get element type
		elemType, ok := gmshElementType4[gmshType]
		if !ok {
			// Skip unknown element types
			for j := 0; j < numElemsInBlock; j++ {
				scanner.Scan()
			}
			continue
		}

		// Get physical tags from entity
		var physicalTags []int
		if entity, ok := entities[strconv.Itoa(entityDim)][entityTag]; ok {
			physicalTags = entity.PhysicalTags
		}

		// Read elements
		expectedNodes := elemType.GetNumNodes()
		for j := 0; j < numElemsInBlock; j++ {
			if !scanner.Scan() {
				return fmt.Errorf("unexpected EOF reading elements")
			}

			fields := strings.Fields(scanner.Text())
			if len(fields) < 1+expectedNodes {
				return fmt.Errorf("invalid element line: expected %d fields, got %d",
					1+expectedNodes, len(fields))
			}

			elemTag, _ := strconv.Atoi(fields[0])

			// Build tags array
			tags := []int{}
			if len(physicalTags) > 0 {
				tags = append(tags, physicalTags[0]) // Primary physical tag
			}
			tags = append(tags, entityTag) // Elementary entity tag

			// Read nodes
			nodeIDs := make([]int, expectedNodes)
			for k := 0; k < expectedNodes; k++ {
				nodeIDs[k], _ = strconv.Atoi(fields[1+k])
			}

			if err := mesh.AddElement(elemTag, elemType, tags, nodeIDs); err != nil {
				return err
			}
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

// readPeriodic4 reads periodic boundary conditions in v4 format
func readPeriodic4(scanner *bufio.Scanner, mesh *Mesh) error {
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

		// Format: entityDim entityTag entityTagMaster
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

		// Read number of corresponding nodes
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF in Periodic nodes")
		}

		numNodes, err := strconv.Atoi(scanner.Text())
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

		// Check for affine transformation
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF after periodic nodes")
		}

		line := scanner.Text()
		if strings.Contains(line, "Affine") {
			// Parse affine transformation
			affineFields := strings.Fields(line)[1:] // Skip "Affine"
			if len(affineFields) >= 16 {
				periodic.AffineTransform = make([]float64, 16)
				for k := 0; k < 16; k++ {
					periodic.AffineTransform[k], _ = strconv.ParseFloat(affineFields[k], 64)
				}
			}
		} else {
			// Put the line back by checking if it's the next section
			if strings.HasPrefix(strings.TrimSpace(line), "$") {
				// This is the next section, we're done with this periodic entity
				mesh.Periodics = append(mesh.Periodics, periodic)
				if strings.TrimSpace(line) == "$EndPeriodic" {
					return nil
				}
				// Otherwise, we've read into the next section, which is an error
				// in our simple scanner approach. In practice, you'd need to
				// handle this more gracefully.
				continue
			}
		}

		mesh.Periodics = append(mesh.Periodics, periodic)
	}

	// Skip to end if not already there
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndPeriodic" {
			break
		}
	}

	return nil
}

// readGmsh4Binary reads binary format for v4
func readGmsh4Binary(file *os.File, mesh *Mesh) (*Mesh, error) {
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

	// Store entities
	entities := make(map[string]map[int]*EntityInfo)
	entities["0"] = make(map[int]*EntityInfo)
	entities["1"] = make(map[int]*EntityInfo)
	entities["2"] = make(map[int]*EntityInfo)
	entities["3"] = make(map[int]*EntityInfo)

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

		case "$PhysicalNames":
			// Physical names are always ASCII even in binary files
			scanner := bufio.NewScanner(reader)
			if err := readPhysicalNames(scanner, mesh); err != nil {
				return nil, err
			}

		case "$Entities":
			if err := readEntities4Binary(reader, entities, byteOrder); err != nil {
				return nil, err
			}

		case "$Nodes":
			if err := readNodes4Binary(reader, mesh, entities, byteOrder); err != nil {
				return nil, err
			}

		case "$Elements":
			if err := readElements4Binary(reader, mesh, entities, byteOrder); err != nil {
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

// readEntities4Binary reads entities in binary format
func readEntities4Binary(reader *bufio.Reader, entities map[string]map[int]*EntityInfo, byteOrder binary.ByteOrder) error {
	// Read counts line (ASCII)
	line, err := reader.ReadString('\n')
	if err != nil {
		return err
	}

	counts := strings.Fields(line)
	if len(counts) < 4 {
		return fmt.Errorf("invalid entity counts")
	}

	numPoints, _ := strconv.Atoi(counts[0])
	numCurves, _ := strconv.Atoi(counts[1])
	numSurfaces, _ := strconv.Atoi(counts[2])
	numVolumes, _ := strconv.Atoi(counts[3])

	// For binary format, entity data is more complex
	// For now, we'll skip the binary entity reading
	// and just create empty entity maps
	_ = numPoints
	_ = numCurves
	_ = numSurfaces
	_ = numVolumes

	// Skip to $EndEntities
	for {
		line, err := reader.ReadString('\n')
		if err != nil {
			return err
		}
		if strings.TrimSpace(line) == "$EndEntities" {
			break
		}
	}

	return nil
}

// readNodes4Binary reads nodes in binary format for v4
func readNodes4Binary(reader *bufio.Reader, mesh *Mesh, entities map[string]map[int]*EntityInfo, byteOrder binary.ByteOrder) error {
	// Read header line (ASCII)
	line, err := reader.ReadString('\n')
	if err != nil {
		return err
	}

	header := strings.Fields(line)
	if len(header) < 4 {
		return fmt.Errorf("invalid Nodes header")
	}

	numEntityBlocks, _ := strconv.Atoi(header[0])
	totalNodes, _ := strconv.Atoi(header[1])

	// Pre-allocate
	mesh.Vertices = make([][]float64, 0, totalNodes)

	// Read entity blocks
	for i := 0; i < numEntityBlocks; i++ {
		// Read block header (ASCII)
		line, err := reader.ReadString('\n')
		if err != nil {
			return err
		}

		blockHeader := strings.Fields(line)
		if len(blockHeader) < 4 {
			return fmt.Errorf("invalid node block header")
		}

		// entityDim, _ := strconv.Atoi(blockHeader[0])
		// entityTag, _ := strconv.Atoi(blockHeader[1])
		parametric, _ := strconv.Atoi(blockHeader[2])
		numNodesInBlock, _ := strconv.Atoi(blockHeader[3])

		// Read binary node data
		// First read all node tags
		nodeTags := make([]int64, numNodesInBlock)
		if mesh.DataSize == 8 {
			if err := binary.Read(reader, byteOrder, &nodeTags); err != nil {
				return err
			}
		} else {
			// 4-byte tags
			tags32 := make([]int32, numNodesInBlock)
			if err := binary.Read(reader, byteOrder, &tags32); err != nil {
				return err
			}
			for j, tag := range tags32 {
				nodeTags[j] = int64(tag)
			}
		}

		// Then read all coordinates
		for j := 0; j < numNodesInBlock; j++ {
			coords := make([]float64, 3)
			if err := binary.Read(reader, byteOrder, &coords); err != nil {
				return err
			}

			// Skip parametric coordinates if present
			if parametric > 0 {
				paramCoords := make([]float64, parametric)
				if err := binary.Read(reader, byteOrder, &paramCoords); err != nil {
					return err
				}
			}

			mesh.AddNode(int(nodeTags[j]), coords)
		}
	}

	// Read $EndNodes
	line, err = reader.ReadString('\n')
	if err != nil {
		return err
	}

	return nil
}

// readElements4Binary reads elements in binary format for v4
func readElements4Binary(reader *bufio.Reader, mesh *Mesh, entities map[string]map[int]*EntityInfo, byteOrder binary.ByteOrder) error {
	// Read header line (ASCII)
	line, err := reader.ReadString('\n')
	if err != nil {
		return err
	}

	header := strings.Fields(line)
	if len(header) < 4 {
		return fmt.Errorf("invalid Elements header")
	}

	numEntityBlocks, _ := strconv.Atoi(header[0])

	// Read entity blocks
	for i := 0; i < numEntityBlocks; i++ {
		// Read block header (ASCII)
		line, err := reader.ReadString('\n')
		if err != nil {
			return err
		}

		blockHeader := strings.Fields(line)
		if len(blockHeader) < 4 {
			return fmt.Errorf("invalid element block header")
		}

		entityDim, _ := strconv.Atoi(blockHeader[0])
		entityTag, _ := strconv.Atoi(blockHeader[1])
		gmshType, _ := strconv.Atoi(blockHeader[2])
		numElemsInBlock, _ := strconv.Atoi(blockHeader[3])

		// Get element type
		elemType, ok := gmshElementType4[gmshType]
		if !ok {
			// Skip unknown element types
			// Need to read and discard the binary data
			expectedNodes := getNodeCountForGmshType(gmshType)
			if expectedNodes == 0 {
				return fmt.Errorf("unknown element type %d", gmshType)
			}

			for j := 0; j < numElemsInBlock; j++ {
				// Read element tag
				if mesh.DataSize == 8 {
					var tag int64
					binary.Read(reader, byteOrder, &tag)
				} else {
					var tag int32
					binary.Read(reader, byteOrder, &tag)
				}

				// Read nodes
				for k := 0; k < expectedNodes; k++ {
					if mesh.DataSize == 8 {
						var node int64
						binary.Read(reader, byteOrder, &node)
					} else {
						var node int32
						binary.Read(reader, byteOrder, &node)
					}
				}
			}
			continue
		}

		// Get physical tags from entity
		var physicalTags []int
		if entity, ok := entities[strconv.Itoa(entityDim)][entityTag]; ok {
			physicalTags = entity.PhysicalTags
		}

		// Read elements
		expectedNodes := elemType.GetNumNodes()
		for j := 0; j < numElemsInBlock; j++ {
			// Read element tag
			var elemTag int
			if mesh.DataSize == 8 {
				var tag int64
				if err := binary.Read(reader, byteOrder, &tag); err != nil {
					return err
				}
				elemTag = int(tag)
			} else {
				var tag int32
				if err := binary.Read(reader, byteOrder, &tag); err != nil {
					return err
				}
				elemTag = int(tag)
			}

			// Build tags array
			tags := []int{}
			if len(physicalTags) > 0 {
				tags = append(tags, physicalTags[0])
			}
			tags = append(tags, entityTag)

			// Read nodes
			nodeIDs := make([]int, expectedNodes)
			for k := 0; k < expectedNodes; k++ {
				if mesh.DataSize == 8 {
					var node int64
					if err := binary.Read(reader, byteOrder, &node); err != nil {
						return err
					}
					nodeIDs[k] = int(node)
				} else {
					var node int32
					if err := binary.Read(reader, byteOrder, &node); err != nil {
						return err
					}
					nodeIDs[k] = int(node)
				}
			}

			if err := mesh.AddElement(elemTag, elemType, tags, nodeIDs); err != nil {
				return err
			}
		}
	}

	// Read $EndElements
	line, err = reader.ReadString('\n')
	if err != nil {
		return err
	}

	return nil
}

// getNodeCountForGmshType returns the number of nodes for a Gmsh element type
// This is used when skipping unknown element types in binary format
func getNodeCountForGmshType(gmshType int) int {
	// Based on Gmsh documentation
	nodeCounts := map[int]int{
		1:  2,  // 2-node line
		2:  3,  // 3-node triangle
		3:  4,  // 4-node quadrangle
		4:  4,  // 4-node tetrahedron
		5:  8,  // 8-node hexahedron
		6:  6,  // 6-node prism
		7:  5,  // 5-node pyramid
		8:  3,  // 3-node second order line
		9:  6,  // 6-node second order triangle
		10: 9,  // 9-node second order quadrangle
		11: 10, // 10-node second order tetrahedron
		12: 27, // 27-node second order hexahedron
		13: 18, // 18-node second order prism
		14: 14, // 14-node second order pyramid
		15: 1,  // 1-node point
		16: 8,  // 8-node second order quadrangle
		17: 20, // 20-node second order hexahedron
		18: 15, // 15-node second order prism
		19: 13, // 13-node second order pyramid
		20: 9,  // 9-node third order triangle
		21: 10, // 10-node third order triangle
		// Add more as needed
	}

	if count, ok := nodeCounts[gmshType]; ok {
		return count
	}
	return 0
}
