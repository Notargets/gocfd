package mesh

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

// ElementType represents different element types
type ElementType int

const (
	Line ElementType = iota
	Triangle
	Quad
	Tet
	Hex
	Prism
	Pyramid
)

func (e ElementType) String() string {
	return [...]string{"Line", "Triangle", "Quad", "Tet", "Hex", "Prism", "Pyramid"}[e]
}

// Face represents a face of an element
type Face struct {
	Vertices []int // Sorted vertex indices
	Element  int   // Parent element
	LocalID  int   // Local face ID within element
}

// Mesh represents a complete unstructured mesh with all connectivity
type Mesh struct {
	// Geometry
	Vertices [][]float64 // Vertex coordinates [nvertices][3]

	// Element data
	Elements     [][]int       // Element to vertex connectivity [nelems][nverts_per_elem]
	ElementTypes []ElementType // Element type for each element
	ElementTags  []int         // Physical group/tag for each element

	// Connectivity (built during initialization)
	EToE [][]int // Element to element connectivity [nelems][nfaces_per_elem]
	EToF [][]int // Element to face connectivity [nelems][nfaces_per_elem]
	EToP []int   // Element to partition mapping (set after partitioning)

	// Face data
	Faces        []Face         // All unique faces in mesh
	FaceMap      map[string]int // Map from sorted vertex string to face ID
	BoundaryTags map[int]string // Boundary condition tags

	// Mesh statistics
	NumElements int
	NumVertices int
	NumFaces    int
}

// NewMesh creates a new mesh and builds connectivity
func NewMesh() *Mesh {
	return &Mesh{
		FaceMap:      make(map[string]int),
		BoundaryTags: make(map[int]string),
	}
}

// ReadMeshFile reads a mesh file based on extension
func ReadMeshFile(filename string) (*Mesh, error) {
	ext := strings.ToLower(filepath.Ext(filename))

	switch ext {
	case ".neu":
		return ReadGambitNeutral(filename)
	case ".msh":
		return ReadGmsh(filename)
	case ".su2":
		return ReadSU2(filename)
	default:
		return nil, fmt.Errorf("unsupported mesh format: %s", ext)
	}
}

// ReadGambitNeutral reads a Gambit neutral file
func ReadGambitNeutral(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	mesh := NewMesh()
	scanner := bufio.NewScanner(file)

	var numnp, nelem int

	// Read header
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if strings.Contains(line, "NUMNP") {
			fields := strings.Fields(line)
			numnp, _ = strconv.Atoi(fields[1])
			nelem, _ = strconv.Atoi(fields[2])
			break
		}
	}

	// Read sections
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if strings.HasPrefix(line, "NODAL COORDINATES") {
			scanner.Scan() // Skip version line
			mesh.Vertices = make([][]float64, numnp)

			for i := 0; i < numnp; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())
				id, _ := strconv.Atoi(fields[0])
				x, _ := strconv.ParseFloat(fields[1], 64)
				y, _ := strconv.ParseFloat(fields[2], 64)
				z, _ := strconv.ParseFloat(fields[3], 64)
				mesh.Vertices[id-1] = []float64{x, y, z}
			}

		} else if strings.HasPrefix(line, "ELEMENTS/CELLS") {
			scanner.Scan() // Skip version line
			mesh.Elements = make([][]int, nelem)
			mesh.ElementTypes = make([]ElementType, nelem)
			mesh.ElementTags = make([]int, nelem)

			for i := 0; i < nelem; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())

				id, _ := strconv.Atoi(fields[0])
				elemType, _ := strconv.Atoi(fields[1])
				numNodes, _ := strconv.Atoi(fields[2])

				// Map Gambit element types
				switch elemType {
				case 1: // Bar
					mesh.ElementTypes[id-1] = Line
				case 2: // Tri
					mesh.ElementTypes[id-1] = Triangle
				case 3: // Quad
					mesh.ElementTypes[id-1] = Quad
				case 4: // Hex
					mesh.ElementTypes[id-1] = Hex
				case 5: // Wedge/Prism
					mesh.ElementTypes[id-1] = Prism
				case 6: // Tet
					mesh.ElementTypes[id-1] = Tet
				case 7: // Pyramid
					mesh.ElementTypes[id-1] = Pyramid
				}

				// Read vertices (convert to 0-indexed)
				verts := make([]int, numNodes)
				for j := 0; j < numNodes; j++ {
					v, _ := strconv.Atoi(fields[3+j])
					verts[j] = v - 1
				}
				mesh.Elements[id-1] = verts
			}
		}
	}

	mesh.NumElements = len(mesh.Elements)
	mesh.NumVertices = len(mesh.Vertices)
	mesh.BuildConnectivity()

	return mesh, nil
}

// ReadGmsh reads a Gmsh format file (version 2.2)
func ReadGmsh(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	mesh := NewMesh()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		switch line {
		case "$Nodes":
			scanner.Scan()
			numNodes, _ := strconv.Atoi(scanner.Text())
			mesh.Vertices = make([][]float64, numNodes)

			for i := 0; i < numNodes; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())
				id, _ := strconv.Atoi(fields[0])
				x, _ := strconv.ParseFloat(fields[1], 64)
				y, _ := strconv.ParseFloat(fields[2], 64)
				z, _ := strconv.ParseFloat(fields[3], 64)
				mesh.Vertices[id-1] = []float64{x, y, z}
			}

		case "$Elements":
			scanner.Scan()
			numElems, _ := strconv.Atoi(scanner.Text())

			// First pass: count volume elements
			elemData := [][]string{}
			for i := 0; i < numElems; i++ {
				scanner.Scan()
				elemData = append(elemData, strings.Fields(scanner.Text()))
			}

			// Second pass: extract volume elements only
			for _, fields := range elemData {
				elemType, _ := strconv.Atoi(fields[1])

				// Only process 3D elements
				var etype ElementType
				var numNodes int
				switch elemType {
				case 4: // 4-node tet
					etype = Tet
					numNodes = 4
				case 5: // 8-node hex
					etype = Hex
					numNodes = 8
				case 6: // 6-node prism
					etype = Prism
					numNodes = 6
				case 7: // 5-node pyramid
					etype = Pyramid
					numNodes = 5
				case 11: // 10-node tet (2nd order)
					etype = Tet
					numNodes = 4 // Only use corner nodes
				default:
					continue // Skip non-3D elements
				}

				numTags, _ := strconv.Atoi(fields[2])
				physTag := 0
				if numTags > 0 {
					physTag, _ = strconv.Atoi(fields[3])
				}

				// Read vertices (convert to 0-indexed)
				verts := make([]int, numNodes)
				offset := 3 + numTags
				for j := 0; j < numNodes; j++ {
					v, _ := strconv.Atoi(fields[offset+j])
					verts[j] = v - 1
				}

				mesh.Elements = append(mesh.Elements, verts)
				mesh.ElementTypes = append(mesh.ElementTypes, etype)
				mesh.ElementTags = append(mesh.ElementTags, physTag)
			}
		}
	}

	mesh.NumElements = len(mesh.Elements)
	mesh.NumVertices = len(mesh.Vertices)
	mesh.BuildConnectivity()

	return mesh, nil
}

// ReadSU2 reads an SU2 native format file
func ReadSU2(filename string) (*Mesh, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	mesh := NewMesh()
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		if strings.HasPrefix(line, "NELEM=") {
			var nelem int
			fmt.Sscanf(line, "NELEM=%d", &nelem)

			mesh.Elements = make([][]int, nelem)
			mesh.ElementTypes = make([]ElementType, nelem)
			mesh.ElementTags = make([]int, nelem)

			for i := 0; i < nelem; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())

				su2Type, _ := strconv.Atoi(fields[0])

				// Map SU2 element types
				var etype ElementType
				var numNodes int
				switch su2Type {
				case 3: // Line
					etype = Line
					numNodes = 2
				case 5: // Triangle
					etype = Triangle
					numNodes = 3
				case 9: // Quad
					etype = Quad
					numNodes = 4
				case 10: // Tet
					etype = Tet
					numNodes = 4
				case 12: // Hex
					etype = Hex
					numNodes = 8
				case 13: // Prism
					etype = Prism
					numNodes = 6
				case 14: // Pyramid
					etype = Pyramid
					numNodes = 5
				}

				// Read vertices
				verts := make([]int, numNodes)
				for j := 0; j < numNodes; j++ {
					verts[j], _ = strconv.Atoi(fields[1+j])
				}

				// Element ID is last field
				elemID, _ := strconv.Atoi(fields[len(fields)-1])

				mesh.Elements[elemID] = verts
				mesh.ElementTypes[elemID] = etype
			}

		} else if strings.HasPrefix(line, "NPOIN=") {
			var npoin int
			fmt.Sscanf(line, "NPOIN=%d", &npoin)

			mesh.Vertices = make([][]float64, npoin)

			for i := 0; i < npoin; i++ {
				scanner.Scan()
				fields := strings.Fields(scanner.Text())

				x, _ := strconv.ParseFloat(fields[0], 64)
				y, _ := strconv.ParseFloat(fields[1], 64)
				z := 0.0
				if len(fields) > 3 {
					z, _ = strconv.ParseFloat(fields[2], 64)
				}

				// Point ID is last field
				ptID, _ := strconv.Atoi(fields[len(fields)-1])
				mesh.Vertices[ptID] = []float64{x, y, z}
			}
		}
	}

	mesh.NumElements = len(mesh.Elements)
	mesh.NumVertices = len(mesh.Vertices)
	mesh.BuildConnectivity()

	return mesh, nil
}

// BuildConnectivity builds element-to-element and face connectivity
func (m *Mesh) BuildConnectivity() {
	m.EToE = make([][]int, m.NumElements)
	m.EToF = make([][]int, m.NumElements)

	// Build face connectivity
	for elemID := 0; elemID < m.NumElements; elemID++ {
		elemType := m.ElementTypes[elemID]
		vertices := m.Elements[elemID]

		// Get faces for this element type
		faceVertices := GetElementFaces(elemType, vertices)

		m.EToE[elemID] = make([]int, len(faceVertices))
		m.EToF[elemID] = make([]int, len(faceVertices))

		// Initialize to -1 (boundary)
		for i := range m.EToE[elemID] {
			m.EToE[elemID][i] = -1
			m.EToF[elemID][i] = -1
		}

		// Process each face
		for localFaceID, faceVerts := range faceVertices {
			// Create sorted vertex key for face
			sorted := make([]int, len(faceVerts))
			copy(sorted, faceVerts)
			sort.Ints(sorted)

			key := fmt.Sprintf("%v", sorted)

			if faceID, exists := m.FaceMap[key]; exists {
				// Face already exists - this is an interior face
				face := &m.Faces[faceID]
				neighborElem := face.Element
				neighborLocalID := face.LocalID

				// Set connectivity
				m.EToE[elemID][localFaceID] = neighborElem
				m.EToE[neighborElem][neighborLocalID] = elemID

				m.EToF[elemID][localFaceID] = faceID
				m.EToF[neighborElem][neighborLocalID] = faceID
			} else {
				// New face
				face := Face{
					Vertices: sorted,
					Element:  elemID,
					LocalID:  localFaceID,
				}

				faceID := len(m.Faces)
				m.Faces = append(m.Faces, face)
				m.FaceMap[key] = faceID
				m.EToF[elemID][localFaceID] = faceID
			}
		}
	}

	m.NumFaces = len(m.Faces)
}

// GetElementFaces returns the face vertices for each element type
func GetElementFaces(elemType ElementType, vertices []int) [][]int {
	switch elemType {
	case Tet:
		return [][]int{
			{vertices[0], vertices[2], vertices[1]}, // Face 0
			{vertices[0], vertices[1], vertices[3]}, // Face 1
			{vertices[1], vertices[2], vertices[3]}, // Face 2
			{vertices[0], vertices[3], vertices[2]}, // Face 3
		}
	case Hex:
		return [][]int{
			{vertices[0], vertices[3], vertices[2], vertices[1]}, // Face 0 (bottom)
			{vertices[4], vertices[5], vertices[6], vertices[7]}, // Face 1 (top)
			{vertices[0], vertices[1], vertices[5], vertices[4]}, // Face 2
			{vertices[1], vertices[2], vertices[6], vertices[5]}, // Face 3
			{vertices[2], vertices[3], vertices[7], vertices[6]}, // Face 4
			{vertices[3], vertices[0], vertices[4], vertices[7]}, // Face 5
		}
	case Prism:
		return [][]int{
			{vertices[0], vertices[2], vertices[1]},              // Face 0 (bottom tri)
			{vertices[3], vertices[4], vertices[5]},              // Face 1 (top tri)
			{vertices[0], vertices[1], vertices[4], vertices[3]}, // Face 2 (quad)
			{vertices[1], vertices[2], vertices[5], vertices[4]}, // Face 3 (quad)
			{vertices[2], vertices[0], vertices[3], vertices[5]}, // Face 4 (quad)
		}
	case Pyramid:
		return [][]int{
			{vertices[0], vertices[3], vertices[2], vertices[1]}, // Face 0 (base quad)
			{vertices[0], vertices[1], vertices[4]},              // Face 1 (tri)
			{vertices[1], vertices[2], vertices[4]},              // Face 2 (tri)
			{vertices[2], vertices[3], vertices[4]},              // Face 3 (tri)
			{vertices[3], vertices[0], vertices[4]},              // Face 4 (tri)
		}
	default:
		return [][]int{}
	}
}

// GetFaceArea computes the area of a face
func (m *Mesh) GetFaceArea(faceID int) float64 {
	face := m.Faces[faceID]
	if len(face.Vertices) == 3 {
		// Triangle
		v0 := m.Vertices[face.Vertices[0]]
		v1 := m.Vertices[face.Vertices[1]]
		v2 := m.Vertices[face.Vertices[2]]

		// Cross product
		dx1 := v1[0] - v0[0]
		dy1 := v1[1] - v0[1]
		dz1 := v1[2] - v0[2]

		dx2 := v2[0] - v0[0]
		dy2 := v2[1] - v0[1]
		dz2 := v2[2] - v0[2]

		cx := dy1*dz2 - dz1*dy2
		cy := dz1*dx2 - dx1*dz2
		cz := dx1*dy2 - dy1*dx2

		return 0.5 * math.Sqrt(cx*cx+cy*cy+cz*cz)
	} else if len(face.Vertices) == 4 {
		// Quad - split into two triangles
		// This is approximate for non-planar quads
		area1 := m.GetFaceArea(-1) // Would need to create temp face
		area2 := m.GetFaceArea(-1) // Would need to create temp face
		return area1 + area2
	}
	return 0.0
}

// PrintStatistics prints mesh statistics
func (m *Mesh) PrintStatistics() {
	fmt.Printf("Mesh Statistics:\n")
	fmt.Printf("  Vertices: %d\n", m.NumVertices)
	fmt.Printf("  Elements: %d\n", m.NumElements)
	fmt.Printf("  Faces: %d\n", m.NumFaces)

	// Count element types
	typeCounts := make(map[ElementType]int)
	for _, t := range m.ElementTypes {
		typeCounts[t]++
	}

	fmt.Printf("  Element types:\n")
	for t, count := range typeCounts {
		fmt.Printf("    %s: %d\n", t, count)
	}

	// Count boundary faces
	boundaryFaces := 0
	for i := 0; i < m.NumElements; i++ {
		for _, neighbor := range m.EToE[i] {
			if neighbor < 0 {
				boundaryFaces++
			}
		}
	}
	fmt.Printf("  Boundary faces: %d\n", boundaryFaces)
}
