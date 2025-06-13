package tetrahedra

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// GambitMesh represents the mesh data from a Gambit neutral file
type GambitMesh struct {
	// File header info
	Title     string
	DataTitle string
	Program   string
	Version   string
	Date      string

	// Basic mesh data
	Vertices           []Vertex
	Elements           []Element
	Groups             []ElementGroup
	BoundaryConditions []BoundaryCondition

	// Connectivity arrays
	EToV [][]int // Element to vertex connectivity
	EToE [][]int // Element to element connectivity
	EToF [][]int // Element to face connectivity
	EToP []int   // Element to partition mapping (derived from groups)

	// Mesh statistics
	NumVertices   int
	NumElements   int
	NumGroups     int
	NumBCs        int
	NumDimCoord   int // NDFCD
	NumDimVel     int // NDFVL
	HasPartitions bool
	NumPartitions int
}

// Vertex represents a 3D point
type Vertex struct {
	ID      int
	X, Y, Z float64
}

// Element represents a mesh element
type Element struct {
	ID       int
	Type     int // Element type (e.g., 6 for tetrahedron)
	NumNodes int
	Vertices []int
}

// ElementGroup represents a group of elements
type ElementGroup struct {
	Name        string
	NumElements int
	Material    int
	NFlags      int
	Elements    []int // Element IDs in this group
}

// BoundaryCondition represents boundary condition data
type BoundaryCondition struct {
	Name       string
	Type       int // ITYPE
	NumEntries int // NENTRY
	NumValues  int // NVALUES
	Elements   []BCElement
}

// BCElement represents an element/face pair in a boundary condition
type BCElement struct {
	ElementID int
	FaceType  int
	FaceID    int
}

// ConnectivityArrays holds the element connectivity information
// ReadGambitNeutralFile reads a Gambit neutral file according to the standard format
func ReadGambitNeutralFile(filename string) (*GambitMesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("failed to open file: %w", err)
	}
	defer file.Close()

	reader := bufio.NewReader(file)
	mesh := &GambitMesh{}

	// Read control info section
	if err := mesh.readControlInfo(reader); err != nil {
		return nil, fmt.Errorf("failed to read control info: %w", err)
	}

	// Read nodal coordinates section
	if err := mesh.readNodalCoordinates(reader); err != nil {
		return nil, fmt.Errorf("failed to read nodal coordinates: %w", err)
	}

	// Read elements/cells section
	if err := mesh.readElementsCells(reader); err != nil {
		return nil, fmt.Errorf("failed to read elements/cells: %w", err)
	}

	// Read element groups section (may contain partition info)
	if err := mesh.readElementGroups(reader); err != nil {
		// Groups are optional, don't fail
		mesh.NumGroups = 0
	}

	// Read boundary conditions section
	if err := mesh.readBoundaryConditions(reader); err != nil {
		// BCs are optional, don't fail
		mesh.NumBCs = 0
	}

	// Derive partition information from groups if present
	mesh.derivePartitionInfo()

	// Build connectivity arrays
	mesh.buildConnectivityArrays()

	return mesh, nil
}

// readControlInfo reads the CONTROL INFO section
func (m *GambitMesh) readControlInfo(reader *bufio.Reader) error {
	// Read "CONTROL INFO X.X.X"
	line, err := reader.ReadString('\n')
	if err != nil {
		return err
	}
	if !strings.Contains(line, "CONTROL INFO") {
		return fmt.Errorf("expected CONTROL INFO, got: %s", line)
	}

	// Read title
	m.Title, _ = reader.ReadString('\n')
	m.Title = strings.TrimSpace(m.Title)

	// Read data title
	m.DataTitle, _ = reader.ReadString('\n')
	m.DataTitle = strings.TrimSpace(m.DataTitle)

	// Read PROGRAM line
	progLine, _ := reader.ReadString('\n')
	if fields := strings.Fields(progLine); len(fields) >= 4 {
		m.Program = fields[1]
		m.Version = fields[3]
	}

	// Read date
	m.Date, _ = reader.ReadString('\n')
	m.Date = strings.TrimSpace(m.Date)

	// Read NUMNP NELEM NGRPS NBSETS NDFCD NDFVL line
	_, err = reader.ReadString('\n')
	if err != nil {
		return err
	}

	// Read the actual values line
	valuesLine, err := reader.ReadString('\n')
	if err != nil {
		return err
	}

	values := strings.Fields(valuesLine)
	if len(values) >= 6 {
		m.NumVertices, _ = strconv.Atoi(values[0])
		m.NumElements, _ = strconv.Atoi(values[1])
		m.NumGroups, _ = strconv.Atoi(values[2])
		m.NumBCs, _ = strconv.Atoi(values[3])
		m.NumDimCoord, _ = strconv.Atoi(values[4])
		m.NumDimVel, _ = strconv.Atoi(values[5])
	}

	// Read ENDOFSECTION
	endLine, _ := reader.ReadString('\n')
	if !strings.Contains(endLine, "ENDOFSECTION") {
		return fmt.Errorf("expected ENDOFSECTION, got: %s", endLine)
	}

	return nil
}

// readNodalCoordinates reads the NODAL COORDINATES section
func (m *GambitMesh) readNodalCoordinates(reader *bufio.Reader) error {
	// Read "NODAL COORDINATES X.X.X"
	line, err := reader.ReadString('\n')
	if err != nil {
		return err
	}
	if !strings.Contains(line, "NODAL COORDINATES") {
		return fmt.Errorf("expected NODAL COORDINATES, got: %s", line)
	}

	m.Vertices = make([]Vertex, m.NumVertices)

	for i := 0; i < m.NumVertices; i++ {
		line, err := reader.ReadString('\n')
		if err != nil {
			return err
		}

		fields := strings.Fields(line)
		if len(fields) < 4 {
			return fmt.Errorf("invalid vertex line: %s", line)
		}

		id, _ := strconv.Atoi(fields[0])
		x, _ := strconv.ParseFloat(fields[1], 64)
		y, _ := strconv.ParseFloat(fields[2], 64)
		z, _ := strconv.ParseFloat(fields[3], 64)

		m.Vertices[i] = Vertex{
			ID: id,
			X:  x,
			Y:  y,
			Z:  z,
		}
	}

	// Read ENDOFSECTION
	endLine, _ := reader.ReadString('\n')
	if !strings.Contains(endLine, "ENDOFSECTION") {
		return fmt.Errorf("expected ENDOFSECTION after vertices, got: %s", endLine)
	}

	return nil
}

// readElementsCells reads the ELEMENTS/CELLS section
func (m *GambitMesh) readElementsCells(reader *bufio.Reader) error {
	// Read "ELEMENTS/CELLS X.X.X"
	line, err := reader.ReadString('\n')
	if err != nil {
		return err
	}
	if !strings.Contains(line, "ELEMENTS/CELLS") {
		return fmt.Errorf("expected ELEMENTS/CELLS, got: %s", line)
	}

	m.Elements = make([]Element, m.NumElements)
	m.EToV = make([][]int, m.NumElements)

	// Create ID to index mapping for elements
	elemIDToIndex := make(map[int]int)

	for i := 0; i < m.NumElements; i++ {
		line, err := reader.ReadString('\n')
		if err != nil {
			return err
		}

		fields := strings.Fields(line)
		if len(fields) < 3 {
			return fmt.Errorf("invalid element line: %s", line)
		}

		id, _ := strconv.Atoi(fields[0])
		elemType, _ := strconv.Atoi(fields[1])
		numNodes, _ := strconv.Atoi(fields[2])

		if len(fields) < 3+numNodes {
			return fmt.Errorf("insufficient nodes in element line: %s", line)
		}

		vertices := make([]int, numNodes)
		for j := 0; j < numNodes; j++ {
			v, _ := strconv.Atoi(fields[3+j])
			vertices[j] = v - 1 // Convert to 0-based indexing
		}

		m.Elements[i] = Element{
			ID:       id,
			Type:     elemType,
			NumNodes: numNodes,
			Vertices: vertices,
		}

		m.EToV[i] = vertices
		elemIDToIndex[id] = i
	}

	// Read ENDOFSECTION
	endLine, _ := reader.ReadString('\n')
	if !strings.Contains(endLine, "ENDOFSECTION") {
		return fmt.Errorf("expected ENDOFSECTION after elements, got: %s", endLine)
	}

	return nil
}

// readElementGroups reads ELEMENT GROUP sections
func (m *GambitMesh) readElementGroups(reader *bufio.Reader) error {
	m.Groups = make([]ElementGroup, 0)

	for {
		// Look for "ELEMENT GROUP X.X.X"
		line, err := reader.ReadString('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		if strings.Contains(line, "ELEMENT GROUP") {
			// Read GROUP line
			groupLine, err := reader.ReadString('\n')
			if err != nil {
				return err
			}

			// Parse: GROUP: <NGP> ELEMENTS: <NELGP> MATERIAL: <MTYP> NFLAGS: <NFLAGS>
			var group ElementGroup
			fmt.Sscanf(groupLine, "GROUP: %d ELEMENTS: %d MATERIAL: %d NFLAGS: %d",
				&group.NumElements, &group.NumElements, &group.Material, &group.NFlags)

			// Read group name
			nameLine, _ := reader.ReadString('\n')
			group.Name = strings.TrimSpace(nameLine)

			// Read flags line (usually just contains 0)
			_, _ = reader.ReadString('\n')

			// Read element IDs
			group.Elements = make([]int, 0, group.NumElements)
			elementStr := ""

			for {
				elemLine, err := reader.ReadString('\n')
				if err != nil || strings.Contains(elemLine, "ENDOFSECTION") {
					break
				}
				elementStr += " " + strings.TrimSpace(elemLine)
			}

			// Parse all element IDs
			fields := strings.Fields(elementStr)
			for _, field := range fields {
				if id, err := strconv.Atoi(field); err == nil {
					group.Elements = append(group.Elements, id)
				}
			}

			m.Groups = append(m.Groups, group)
		} else if strings.Contains(line, "BOUNDARY CONDITIONS") {
			// Put back the line for BC reading
			// This is a simplification - in production you'd use a more sophisticated approach
			break
		}
	}

	return nil
}

// readBoundaryConditions reads the BOUNDARY CONDITIONS section
func (m *GambitMesh) readBoundaryConditions(reader *bufio.Reader) error {
	// Read "BOUNDARY CONDITIONS X.X.X"
	line, err := reader.ReadString('\n')
	if err == io.EOF {
		return nil // No boundary conditions
	}
	if err != nil {
		return err
	}

	if !strings.Contains(line, "BOUNDARY CONDITIONS") {
		return nil // No boundary conditions
	}

	m.BoundaryConditions = make([]BoundaryCondition, 0)

	for i := 0; i < m.NumBCs; i++ {
		// Read BC header line
		bcLine, err := reader.ReadString('\n')
		if err != nil {
			return err
		}

		fields := strings.Fields(bcLine)
		if len(fields) < 4 {
			continue
		}

		bc := BoundaryCondition{
			Name: fields[0],
		}
		bc.Type, _ = strconv.Atoi(fields[1])
		bc.NumEntries, _ = strconv.Atoi(fields[2])
		bc.NumValues, _ = strconv.Atoi(fields[3])

		// Read element/face pairs
		bc.Elements = make([]BCElement, bc.NumEntries)
		for j := 0; j < bc.NumEntries; j++ {
			elemLine, _ := reader.ReadString('\n')
			elemFields := strings.Fields(elemLine)

			if len(elemFields) >= 3 {
				elem := BCElement{}
				elem.ElementID, _ = strconv.Atoi(elemFields[0])
				elem.ElementID-- // Convert to 0-based
				elem.FaceType, _ = strconv.Atoi(elemFields[1])
				elem.FaceID, _ = strconv.Atoi(elemFields[2])
				elem.FaceID-- // Convert to 0-based

				bc.Elements[j] = elem
			}
		}

		m.BoundaryConditions = append(m.BoundaryConditions, bc)
	}

	return nil
}

// derivePartitionInfo derives partition information from element groups
func (m *GambitMesh) derivePartitionInfo() {
	// Look for groups named p1, p2, p3, etc.
	partitionGroups := make(map[int][]int)

	for _, group := range m.Groups {
		if strings.HasPrefix(group.Name, "p") {
			// Try to extract partition number
			partNumStr := strings.TrimPrefix(group.Name, "p")
			if partNum, err := strconv.Atoi(partNumStr); err == nil {
				partitionGroups[partNum-1] = group.Elements // Convert to 0-based
			}
		}
	}

	if len(partitionGroups) > 0 {
		m.HasPartitions = true
		m.NumPartitions = len(partitionGroups)
		m.EToP = make([]int, m.NumElements)

		// Initialize all elements to partition -1 (unassigned)
		for i := range m.EToP {
			m.EToP[i] = -1
		}

		// Create element ID to index mapping
		elemIDToIndex := make(map[int]int)
		for i, elem := range m.Elements {
			elemIDToIndex[elem.ID] = i
		}

		// Assign partitions based on groups
		for partNum, elemIDs := range partitionGroups {
			for _, elemID := range elemIDs {
				if idx, exists := elemIDToIndex[elemID]; exists {
					m.EToP[idx] = partNum
				}
			}
		}
	}
}

// buildConnectivityArrays builds the EToE and EToF arrays
func (m *GambitMesh) buildConnectivityArrays() {
	K := m.NumElements

	// Face to vertex mapping for different element types
	faceVertices := map[int][][]int{
		6: { // Tetrahedron (4 nodes)
			{0, 1, 2}, // Face 0
			{0, 1, 3}, // Face 1
			{1, 2, 3}, // Face 2
			{0, 2, 3}, // Face 3
		},
		// Add more element types as needed
	}

	// Create face to element+face mapping
	faces := make(map[[3]int]struct{ elem, face int })

	for k := 0; k < K; k++ {
		elemType := m.Elements[k].Type
		faceVerts, exists := faceVertices[elemType]
		if !exists {
			continue // Skip unsupported element types
		}

		Nfaces := len(faceVerts)

		for f := 0; f < Nfaces; f++ {
			// Get vertices of this face
			v := make([]int, len(faceVerts[f]))
			for i := 0; i < len(faceVerts[f]); i++ {
				v[i] = m.EToV[k][faceVerts[f][i]]
			}

			// Sort vertices to create unique key (only for triangular faces)
			if len(v) == 3 {
				sort.Ints(v)
				key := [3]int{v[0], v[1], v[2]}

				if match, exists := faces[key]; exists {
					// Found matching face
					if k != match.elem {
						// Different elements share this face
					}
				} else {
					faces[key] = struct{ elem, face int }{k, f}
				}
			}
		}
	}

	// Build connectivity arrays
	m.EToE = make([][]int, K)
	m.EToF = make([][]int, K)

	for k := 0; k < K; k++ {
		elemType := m.Elements[k].Type
		faceVerts, exists := faceVertices[elemType]
		if !exists {
			continue
		}

		Nfaces := len(faceVerts)
		m.EToE[k] = make([]int, Nfaces)
		m.EToF[k] = make([]int, Nfaces)

		// Initialize with self-connectivity
		for f := 0; f < Nfaces; f++ {
			m.EToE[k][f] = k
			m.EToF[k][f] = f
		}

		// Find actual connections
		for f := 0; f < Nfaces; f++ {
			v := make([]int, len(faceVerts[f]))
			for i := 0; i < len(faceVerts[f]); i++ {
				v[i] = m.EToV[k][faceVerts[f][i]]
			}

			if len(v) == 3 {
				sort.Ints(v)
				key := [3]int{v[0], v[1], v[2]}

				// Search all faces for match
				for fkey, fdata := range faces {
					if fkey == key && fdata.elem != k {
						m.EToE[k][f] = fdata.elem
						m.EToF[k][f] = fdata.face
						break
					}
				}
			}
		}
	}
}

// GetConnectivityArrays returns the connectivity arrays
func (m *GambitMesh) GetConnectivityArrays() *ConnectivityArrays {
	return &ConnectivityArrays{
		EToE: m.EToE,
		EToF: m.EToF,
	}
}

// PrintMeshInfo prints summary information about the mesh
func (m *GambitMesh) PrintMeshInfo() {
	fmt.Printf("Gambit Mesh Information:\n")
	fmt.Printf("  Title: %s\n", m.Title)
	fmt.Printf("  Program: %s %s\n", m.Program, m.Version)
	fmt.Printf("  Vertices: %d\n", m.NumVertices)
	fmt.Printf("  Elements: %d\n", m.NumElements)
	fmt.Printf("  Groups: %d\n", m.NumGroups)
	fmt.Printf("  Boundary Conditions: %d\n", m.NumBCs)

	if m.HasPartitions {
		fmt.Printf("  Partitions: %d\n", m.NumPartitions)

		// Count elements per partition
		partitionCounts := make(map[int]int)
		for _, p := range m.EToP {
			if p >= 0 {
				partitionCounts[p]++
			}
		}

		for p := 0; p < m.NumPartitions; p++ {
			fmt.Printf("    Partition %d: %d elements\n", p, partitionCounts[p])
		}
	} else {
		fmt.Printf("  Partitions: None\n")
	}

	// Print element groups
	fmt.Printf("\n  Element Groups:\n")
	for _, group := range m.Groups {
		fmt.Printf("    %s: %d elements\n", group.Name, len(group.Elements))
	}

	// Print boundary conditions
	if len(m.BoundaryConditions) > 0 {
		fmt.Printf("\n  Boundary Conditions:\n")
		for _, bc := range m.BoundaryConditions {
			fmt.Printf("    %s: %d faces\n", bc.Name, bc.NumEntries)
		}
	}
}
