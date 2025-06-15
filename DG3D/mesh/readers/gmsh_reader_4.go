package readers

import (
	"bufio"
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"os"
	"strconv"
	"strings"
)

// EntityInfo stores information about geometric entities
type EntityInfo struct {
	Dimension    int
	Tag          int
	BoundingBox  [2][3]float64
	PhysicalTags []int
}

// ReadGmsh4 reads a Gmsh MSH file format version 4.x
func ReadGmsh4(filename string) (*mesh.Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	mesh := mesh.NewMesh()

	// Store entity information for later use
	entities := make(map[string]map[int]*EntityInfo)
	entities["0"] = make(map[int]*EntityInfo) // Points
	entities["1"] = make(map[int]*EntityInfo) // Curves
	entities["2"] = make(map[int]*EntityInfo) // Surfaces
	entities["3"] = make(map[int]*EntityInfo) // Volumes

	hasEntities := false

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

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
			if err := readEntities4(scanner, mesh, entities); err != nil {
				return nil, err
			}
			hasEntities = true

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
			if err := readGhostElements4(scanner, mesh); err != nil {
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
	_ = hasEntities // $Entities section is optional in MSH 4.1

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner error: %v", err)
	}

	mesh.BuildConnectivity()
	return mesh, nil
}

// readMeshFormat4 reads the MeshFormat section for v4
func readMeshFormat4(scanner *bufio.Scanner, mesh *mesh.Mesh) error {
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
func readEntities4(scanner *bufio.Scanner, msh *mesh.Mesh, entities map[string]map[int]*EntityInfo) error {
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

	// Read points
	for i := 0; i < numPoints; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading point entity")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 4 {
			return fmt.Errorf("invalid point entity")
		}

		tag, _ := strconv.Atoi(fields[0])
		x, _ := strconv.ParseFloat(fields[1], 64)
		y, _ := strconv.ParseFloat(fields[2], 64)
		z, _ := strconv.ParseFloat(fields[3], 64)

		entity := &EntityInfo{
			Dimension: 0,
			Tag:       tag,
			BoundingBox: [2][3]float64{
				{x, y, z},
				{x, y, z},
			},
		}

		// Read physical tags if present
		if len(fields) > 4 {
			numPhysTags, _ := strconv.Atoi(fields[4])
			entity.PhysicalTags = make([]int, numPhysTags)
			for j := 0; j < numPhysTags && 5+j < len(fields); j++ {
				entity.PhysicalTags[j], _ = strconv.Atoi(fields[5+j])
			}
		}

		entities["0"][tag] = entity

		// Store in mesh
		meshEntity := &mesh.Entity{
			Dimension:    entity.Dimension,
			Tag:          entity.Tag,
			BoundingBox:  entity.BoundingBox,
			PhysicalTags: entity.PhysicalTags,
		}
		msh.Entities[tag] = meshEntity
	}

	// Read curves
	for i := 0; i < numCurves; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading curve entity")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 8 {
			return fmt.Errorf("invalid curve entity")
		}

		tag, _ := strconv.Atoi(fields[0])
		entity := &EntityInfo{
			Dimension: 1,
			Tag:       tag,
		}

		// Bounding box
		for j := 0; j < 3; j++ {
			entity.BoundingBox[0][j], _ = strconv.ParseFloat(fields[1+j], 64)
			entity.BoundingBox[1][j], _ = strconv.ParseFloat(fields[4+j], 64)
		}

		// Physical tags
		pos := 7
		if pos < len(fields) {
			numPhysTags, _ := strconv.Atoi(fields[pos])
			pos++
			entity.PhysicalTags = make([]int, numPhysTags)
			for j := 0; j < numPhysTags && pos+j < len(fields); j++ {
				entity.PhysicalTags[j], _ = strconv.Atoi(fields[pos+j])
			}
			pos += numPhysTags
		}

		// Bounding points
		var boundingEntities []int
		if pos < len(fields) {
			numBoundingPoints, _ := strconv.Atoi(fields[pos])
			pos++
			for j := 0; j < numBoundingPoints && pos+j < len(fields); j++ {
				pointTag, _ := strconv.Atoi(fields[pos+j])
				boundingEntities = append(boundingEntities, pointTag)
			}
		}

		entities["1"][tag] = entity

		// Store in mesh
		meshEntity := &mesh.Entity{
			Dimension:        entity.Dimension,
			Tag:              entity.Tag,
			BoundingBox:      entity.BoundingBox,
			PhysicalTags:     entity.PhysicalTags,
			BoundingEntities: boundingEntities,
		}
		msh.Entities[tag] = meshEntity
	}

	// Read surfaces
	for i := 0; i < numSurfaces; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading surface entity")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 8 {
			return fmt.Errorf("invalid surface entity")
		}

		tag, _ := strconv.Atoi(fields[0])
		entity := &EntityInfo{
			Dimension: 2,
			Tag:       tag,
		}

		// Bounding box
		for j := 0; j < 3; j++ {
			entity.BoundingBox[0][j], _ = strconv.ParseFloat(fields[1+j], 64)
			entity.BoundingBox[1][j], _ = strconv.ParseFloat(fields[4+j], 64)
		}

		// Physical tags
		pos := 7
		if pos < len(fields) {
			numPhysTags, _ := strconv.Atoi(fields[pos])
			pos++
			entity.PhysicalTags = make([]int, numPhysTags)
			for j := 0; j < numPhysTags && pos+j < len(fields); j++ {
				entity.PhysicalTags[j], _ = strconv.Atoi(fields[pos+j])
			}
			pos += numPhysTags
		}

		// Bounding curves
		var boundingEntities []int
		if pos < len(fields) {
			numBoundingCurves, _ := strconv.Atoi(fields[pos])
			pos++
			for j := 0; j < numBoundingCurves && pos+j < len(fields); j++ {
				curveTag, _ := strconv.Atoi(fields[pos+j])
				boundingEntities = append(boundingEntities, curveTag)
			}
		}

		entities["2"][tag] = entity

		// Store in mesh
		meshEntity := &mesh.Entity{
			Dimension:        entity.Dimension,
			Tag:              entity.Tag,
			BoundingBox:      entity.BoundingBox,
			PhysicalTags:     entity.PhysicalTags,
			BoundingEntities: boundingEntities,
		}
		msh.Entities[tag] = meshEntity
	}

	// Read volumes
	for i := 0; i < numVolumes; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading volume entity")
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 8 {
			return fmt.Errorf("invalid volume entity")
		}

		tag, _ := strconv.Atoi(fields[0])
		entity := &EntityInfo{
			Dimension: 3,
			Tag:       tag,
		}

		// Bounding box
		for j := 0; j < 3; j++ {
			entity.BoundingBox[0][j], _ = strconv.ParseFloat(fields[1+j], 64)
			entity.BoundingBox[1][j], _ = strconv.ParseFloat(fields[4+j], 64)
		}

		// Physical tags
		pos := 7
		if pos < len(fields) {
			numPhysTags, _ := strconv.Atoi(fields[pos])
			pos++
			entity.PhysicalTags = make([]int, numPhysTags)
			for j := 0; j < numPhysTags && pos+j < len(fields); j++ {
				entity.PhysicalTags[j], _ = strconv.Atoi(fields[pos+j])
			}
			pos += numPhysTags
		}

		// Bounding surfaces
		var boundingEntities []int
		if pos < len(fields) {
			numBoundingSurfaces, _ := strconv.Atoi(fields[pos])
			pos++
			for j := 0; j < numBoundingSurfaces && pos+j < len(fields); j++ {
				surfaceTag, _ := strconv.Atoi(fields[pos+j])
				boundingEntities = append(boundingEntities, surfaceTag)
			}
		}

		entities["3"][tag] = entity

		// Store in mesh
		meshEntity := &mesh.Entity{
			Dimension:        entity.Dimension,
			Tag:              entity.Tag,
			BoundingBox:      entity.BoundingBox,
			PhysicalTags:     entity.PhysicalTags,
			BoundingEntities: boundingEntities,
		}
		msh.Entities[tag] = meshEntity
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndEntities" {
			break
		}
	}

	return nil
}

// readPartitionedEntities4 reads partitioned entities (for parallel meshes)
func readPartitionedEntities4(scanner *bufio.Scanner, entities map[string]map[int]*EntityInfo) error {
	// For now, skip this section
	return skipSection(scanner, "$EndPartitionedEntities")
}

// readNodes4 reads nodes in v4 format
func readNodes4(scanner *bufio.Scanner, msh *mesh.Mesh, entities map[string]map[int]*EntityInfo) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Nodes")
	}

	// Format: numEntityBlocks numNodes minNodeTag maxNodeTag
	header := strings.Fields(scanner.Text())
	if len(header) < 4 {
		return fmt.Errorf("invalid Nodes header")
	}

	numEntityBlocks, _ := strconv.Atoi(header[0])
	// totalNodes, _ := strconv.Atoi(header[1])

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
			nodeTags[j], _ = strconv.Atoi(strings.TrimSpace(scanner.Text()))
		}

		// Read node coordinates (and parametric coords if present)
		for j := 0; j < numNodesInBlock; j++ {
			if !scanner.Scan() {
				return fmt.Errorf("unexpected EOF reading node coordinates")
			}

			fields := strings.Fields(scanner.Text())
			if len(fields) < 3 {
				return fmt.Errorf("invalid node coordinate line")
			}

			coords := make([]float64, 3)
			for k := 0; k < 3; k++ {
				coords[k], _ = strconv.ParseFloat(fields[k], 64)
			}

			// Handle parametric coordinates if present
			if parametric > 0 && len(fields) >= 3+parametric {
				paramCoords := make([]float64, parametric)
				for k := 0; k < parametric; k++ {
					paramCoords[k], _ = strconv.ParseFloat(fields[3+k], 64)
				}
				msh.AddNodeWithParametric(nodeTags[j], coords, paramCoords)
			} else {
				msh.AddNode(nodeTags[j], coords)
			}
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
func readElements4(scanner *bufio.Scanner, msh *mesh.Mesh, entities map[string]map[int]*EntityInfo) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Elements")
	}

	// Format: numEntityBlocks numElements minElementTag maxElementTag
	header := strings.Fields(scanner.Text())
	if len(header) < 4 {
		return fmt.Errorf("invalid Elements header")
	}

	numEntityBlocks, _ := strconv.Atoi(header[0])

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
				return fmt.Errorf("invalid element line: expected at least %d fields, got %d",
					1+expectedNodes, len(fields))
			}

			elemTag, _ := strconv.Atoi(fields[0])

			// Build tags array
			tags := []int{}
			if len(physicalTags) > 0 {
				tags = append(tags, physicalTags...)
			}
			tags = append(tags, entityTag) // Add geometric entity tag

			// Read nodes
			nodeIDs := make([]int, expectedNodes)
			for k := 0; k < expectedNodes; k++ {
				nodeIDs[k], _ = strconv.Atoi(fields[1+k])
			}

			if err := msh.AddElement(elemTag, elemType, tags, nodeIDs); err != nil {
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
func readPeriodic4(scanner *bufio.Scanner, msh *mesh.Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in Periodic")
	}

	numPeriodicLinks, _ := strconv.Atoi(strings.TrimSpace(scanner.Text()))

	for i := 0; i < numPeriodicLinks; i++ {
		// Read entity dimension and tags
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading periodic link %d", i)
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 3 {
			return fmt.Errorf("invalid periodic entity line")
		}

		entityDim, _ := strconv.Atoi(fields[0])
		slaveTag, _ := strconv.Atoi(fields[1])
		masterTag, _ := strconv.Atoi(fields[2])

		// Read affine transformation size
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading affine size")
		}

		numAffine, _ := strconv.Atoi(strings.TrimSpace(scanner.Text()))

		// Read affine transformation if present
		var affineTransform []float64
		if numAffine > 0 {
			if !scanner.Scan() {
				return fmt.Errorf("unexpected EOF reading affine transform")
			}

			affineFields := strings.Fields(scanner.Text())
			affineTransform = make([]float64, numAffine)
			for j := 0; j < numAffine && j < len(affineFields); j++ {
				affineTransform[j], _ = strconv.ParseFloat(affineFields[j], 64)
			}
		}

		// Read number of corresponding nodes
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading node count")
		}

		numNodes, _ := strconv.Atoi(strings.TrimSpace(scanner.Text()))

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
			Dimension:       entityDim,
			SlaveTag:        slaveTag,
			MasterTag:       masterTag,
			NodeMap:         nodeMap,
			AffineTransform: affineTransform,
		}

		msh.Periodics = append(msh.Periodics, periodic)
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndPeriodic" {
			break
		}
	}

	return nil
}

// readGhostElements4 reads ghost elements for parallel meshes
func readGhostElements4(scanner *bufio.Scanner, msh *mesh.Mesh) error {
	if !scanner.Scan() {
		return fmt.Errorf("unexpected EOF in GhostElements")
	}

	// Read header
	header := strings.Fields(scanner.Text())
	if len(header) < 1 {
		return fmt.Errorf("invalid GhostElements header")
	}

	// For v4.0 format: numGhostCells minElementTag maxElementTag numPartitions partitionTag ...
	// For v4.1 format: numGhostElements

	var numGhostElements int
	if len(header) >= 3 {
		// v4.0 format
		numGhostElements, _ = strconv.Atoi(header[0])
	} else {
		// v4.1 format
		numGhostElements, _ = strconv.Atoi(header[0])
	}

	// Read ghost elements
	for i := 0; i < numGhostElements; i++ {
		if !scanner.Scan() {
			return fmt.Errorf("unexpected EOF reading ghost element %d", i)
		}

		fields := strings.Fields(scanner.Text())
		if len(fields) < 3 {
			return fmt.Errorf("invalid ghost element line")
		}

		elementTag, _ := strconv.Atoi(fields[0])
		ownerPartition, _ := strconv.Atoi(fields[1])
		numGhostPartitions, _ := strconv.Atoi(fields[2])

		ghostPartitions := make([]int, numGhostPartitions)
		for j := 0; j < numGhostPartitions && 3+j < len(fields); j++ {
			ghostPartitions[j], _ = strconv.Atoi(fields[3+j])
		}

		ghost := mesh.GhostElement{
			ElementTag:      elementTag,
			OwnerPartition:  ownerPartition,
			GhostPartitions: ghostPartitions,
		}

		msh.GhostElements = append(msh.GhostElements, ghost)
	}

	// Skip to end
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == "$EndGhostElements" {
			break
		}
	}

	return nil
}

// skipSection skips a section until the end marker
func skipSection(scanner *bufio.Scanner, endMarker string) error {
	for scanner.Scan() {
		if strings.TrimSpace(scanner.Text()) == endMarker {
			return nil
		}
	}
	return fmt.Errorf("unexpected EOF while looking for %s", endMarker)
}

// gmshElementType4 maps Gmsh v4 element type numbers to our ElementType
var gmshElementType4 = map[int]mesh.ElementType{
	1:  mesh.Line,      // 2-node line
	2:  mesh.Triangle,  // 3-node triangle
	3:  mesh.Quad,      // 4-node quadrangle
	4:  mesh.Tet,       // 4-node tetrahedron
	5:  mesh.Hex,       // 8-node hexahedron
	6:  mesh.Prism,     // 6-node prism
	7:  mesh.Pyramid,   // 5-node pyramid
	8:  mesh.Line3,     // 3-node line
	9:  mesh.Triangle6, // 6-node triangle
	10: mesh.Quad9,     // 9-node quadrangle
	11: mesh.Tet10,     // 10-node tetrahedron
	12: mesh.Hex27,     // 27-node hexahedron
	13: mesh.Prism18,   // 18-node prism
	14: mesh.Pyramid14, // 14-node pyramid
	// Add more as needed
}

// Helper function to get node count for unknown Gmsh types
func getNodeCountForGmshType(gmshType int) int {
	// This is a fallback for binary reading
	switch gmshType {
	case 1:
		return 2 // Line
	case 2:
		return 3 // Triangle
	case 3:
		return 4 // Quad
	case 4:
		return 4 // Tet
	case 5:
		return 8 // Hex
	case 6:
		return 6 // Prism
	case 7:
		return 5 // Pyramid
	case 8:
		return 3 // Line3
	case 9:
		return 6 // Triangle6
	case 10:
		return 9 // Quad9
	case 11:
		return 10 // Tet10
	case 12:
		return 27 // Hex27
	case 13:
		return 18 // Prism18
	case 14:
		return 14 // Pyramid14
	default:
		return 0
	}
}
