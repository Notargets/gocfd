package mesh

import (
	"fmt"
	"path/filepath"
	"sort"
	"strings"
)

// ElementType represents different element types
type ElementType int

const (
	// 0D elements
	Point ElementType = iota

	// 1D elements
	Line
	Line3 // 3-node line

	// 2D elements
	Triangle
	Triangle6  // 6-node triangle
	Triangle9  // 9-node triangle
	Triangle10 // 10-node triangle
	Quad
	Quad8 // 8-node quad
	Quad9 // 9-node quad

	// 3D elements
	Tet
	Tet10 // 10-node tetrahedron
	Hex
	Hex20 // 20-node hexahedron
	Hex27 // 27-node hexahedron
	Prism
	Prism15 // 15-node prism
	Prism18 // 18-node prism
	Pyramid
	Pyramid13 // 13-node pyramid
	Pyramid14 // 14-node pyramid
)

func (e ElementType) String() string {
	names := []string{
		"Point", "Line", "Line3",
		"Triangle", "Triangle6", "Triangle9", "Triangle10",
		"Quad", "Quad8", "Quad9",
		"Tet", "Tet10", "Hex", "Hex20", "Hex27",
		"Prism", "Prism15", "Prism18",
		"Pyramid", "Pyramid13", "Pyramid14",
	}
	if int(e) < len(names) {
		return names[e]
	}
	return fmt.Sprintf("Unknown(%d)", e)
}

// GetDimension returns the spatial dimension of the element
func (e ElementType) GetDimension() int {
	switch e {
	case Point:
		return 0
	case Line, Line3:
		return 1
	case Triangle, Triangle6, Triangle9, Triangle10, Quad, Quad8, Quad9:
		return 2
	default:
		return 3
	}
}

// GetNumNodes returns the number of nodes for the element type
func (e ElementType) GetNumNodes() int {
	switch e {
	case Point:
		return 1
	case Line:
		return 2
	case Line3:
		return 3
	case Triangle:
		return 3
	case Triangle6:
		return 6
	case Triangle9:
		return 9
	case Triangle10:
		return 10
	case Quad:
		return 4
	case Quad8:
		return 8
	case Quad9:
		return 9
	case Tet:
		return 4
	case Tet10:
		return 10
	case Hex:
		return 8
	case Hex20:
		return 20
	case Hex27:
		return 27
	case Prism:
		return 6
	case Prism15:
		return 15
	case Prism18:
		return 18
	case Pyramid:
		return 5
	case Pyramid13:
		return 13
	case Pyramid14:
		return 14
	default:
		return 0
	}
}

// IsHigherOrder returns true if this is a higher-order element
func (e ElementType) IsHigherOrder() bool {
	switch e {
	case Line3, Triangle6, Triangle9, Triangle10, Quad8, Quad9,
		Tet10, Hex20, Hex27, Prism15, Prism18, Pyramid13, Pyramid14:
		return true
	default:
		return false
	}
}

// GetCornerNodes returns the indices of corner nodes for higher-order elements
func (e ElementType) GetCornerNodes() []int {
	switch e {
	case Line3:
		return []int{0, 1}
	case Triangle6, Triangle9, Triangle10:
		return []int{0, 1, 2}
	case Quad8, Quad9:
		return []int{0, 1, 2, 3}
	case Tet10:
		return []int{0, 1, 2, 3}
	case Hex20, Hex27:
		return []int{0, 1, 2, 3, 4, 5, 6, 7}
	case Prism15, Prism18:
		return []int{0, 1, 2, 3, 4, 5}
	case Pyramid13, Pyramid14:
		return []int{0, 1, 2, 3, 4}
	default:
		// For linear elements, all nodes are corner nodes
		n := e.GetNumNodes()
		nodes := make([]int, n)
		for i := 0; i < n; i++ {
			nodes[i] = i
		}
		return nodes
	}
}

// Face represents a face of an element
type Face struct {
	Vertices []int // Sorted vertex indices
	Element  int   // Parent element
	LocalID  int   // Local face ID within element
}

// ElementGroup represents a physical or elementary entity group
type ElementGroup struct {
	Dimension int
	Tag       int
	Name      string
	Elements  []int // Element indices in this group
}

// NodeGroup represents a group of nodes (for boundary conditions, etc.)
type NodeGroup struct {
	Tag   int
	Name  string
	Nodes []int // Node indices in this group
}

// Periodic represents periodic boundary condition data
type Periodic struct {
	Dimension       int
	SlaveTag        int
	MasterTag       int
	NodeMap         map[int]int // slave node -> master node
	AffineTransform []float64   // Optional affine transformation matrix
}

// Mesh represents a complete unstructured mesh with all connectivity
type Mesh struct {
	// Geometry
	Vertices     [][]float64 // Vertex coordinates [nvertices][3]
	NodeIDMap    map[int]int // Maps original node IDs to array indices
	NodeArrayMap map[int]int // Maps array indices to original node IDs

	// Element data
	EtoV         [][]int       // Element to vertex connectivity [nelems][nverts_per_elem]
	ElementTypes []ElementType // Element type for each element
	ElementTags  [][]int       // All tags for each element (physical, elementary, etc.)
	ElementIDMap map[int]int   // Maps original element IDs to array indices

	// Groups (physical entities)
	ElementGroups map[int]*ElementGroup // Physical/elementary groups by tag
	NodeGroups    map[int]*NodeGroup    // Node groups by tag

	// Connectivity (built during initialization)
	EToE [][]int // Element to element connectivity [nelems][nfaces_per_elem]
	EToF [][]int // Element to face connectivity [nelems][nfaces_per_elem]
	EToP []int   // Element to partition mapping (set after partitioning)

	// Face data
	Faces   []Face         // All unique faces in mesh
	FaceMap map[string]int // Map from sorted vertex string to face ID

	// Boundary conditions
	BoundaryTags map[int]string // Boundary condition tags
	Periodics    []Periodic     // Periodic boundary conditions

	// Mesh statistics
	NumElements int
	NumVertices int
	NumFaces    int

	// Format-specific data
	FormatVersion string // e.g., "2.2", "4.1"
	IsBinary      bool
	DataSize      int // Size of data fields (for binary)
}

// NewMesh creates a new mesh and initializes maps
func NewMesh() *Mesh {
	return &Mesh{
		NodeIDMap:     make(map[int]int),
		NodeArrayMap:  make(map[int]int),
		ElementIDMap:  make(map[int]int),
		ElementGroups: make(map[int]*ElementGroup),
		NodeGroups:    make(map[int]*NodeGroup),
		FaceMap:       make(map[string]int),
		BoundaryTags:  make(map[int]string),
	}
}

// AddNode adds a node with the given ID and coordinates
func (m *Mesh) AddNode(nodeID int, coords []float64) {
	idx := len(m.Vertices)
	m.Vertices = append(m.Vertices, coords)
	m.NodeIDMap[nodeID] = idx
	m.NodeArrayMap[idx] = nodeID
	m.NumVertices = len(m.Vertices)
}

// AddElement adds an element with the given ID, type, and connectivity
func (m *Mesh) AddElement(elemID int, elemType ElementType, tags []int, nodeIDs []int) error {
	// Convert node IDs to array indices
	nodes := make([]int, len(nodeIDs))
	for i, nid := range nodeIDs {
		// Skip zero node IDs (padding)
		if nid == 0 {
			nodes[i] = -1
			continue
		}
		idx, ok := m.NodeIDMap[nid]
		if !ok {
			return fmt.Errorf("element %d references unknown node %d", elemID, nid)
		}
		nodes[i] = idx
	}

	idx := len(m.EtoV)
	m.EtoV = append(m.EtoV, nodes)
	m.ElementTypes = append(m.ElementTypes, elemType)
	m.ElementTags = append(m.ElementTags, tags)
	m.ElementIDMap[elemID] = idx
	m.NumElements = len(m.EtoV)

	// Add to physical groups if applicable
	if len(tags) > 0 && tags[0] > 0 {
		physTag := tags[0]
		if group, ok := m.ElementGroups[physTag]; ok {
			group.Elements = append(group.Elements, idx)
		}
	}

	return nil
}

// GetNodeIndex returns the array index for a given node ID
func (m *Mesh) GetNodeIndex(nodeID int) (int, bool) {
	idx, ok := m.NodeIDMap[nodeID]
	return idx, ok
}

// GetNodeID returns the original node ID for a given array index
func (m *Mesh) GetNodeID(index int) (int, bool) {
	id, ok := m.NodeArrayMap[index]
	return id, ok
}

// ReadMeshFile reads a mesh file based on extension
func ReadMeshFile(filename string) (*Mesh, error) {
	ext := strings.ToLower(filepath.Ext(filename))

	switch ext {
	case ".neu":
		return ReadGambitNeutral(filename)
	case ".msh":
		// Try to detect version by reading first few lines
		return ReadGmshAuto(filename)
	case ".su2":
		return ReadSU2(filename)
	default:
		return nil, fmt.Errorf("unsupported mesh format: %s", ext)
	}
}

// FilterByDimension returns elements of specified dimension
func (m *Mesh) FilterByDimension(dim int) ([]int, [][]int, []ElementType) {
	var indices []int
	var elements [][]int
	var types []ElementType

	for i, etype := range m.ElementTypes {
		if etype.GetDimension() == dim {
			indices = append(indices, i)
			elements = append(elements, m.EtoV[i])
			types = append(types, etype)
		}
	}

	return indices, elements, types
}

// BuildConnectivity builds element-to-element and face connectivity
func (m *Mesh) BuildConnectivity() {
	m.EToE = make([][]int, m.NumElements)
	m.EToF = make([][]int, m.NumElements)

	// Build face connectivity for 3D elements only
	for elemID := 0; elemID < m.NumElements; elemID++ {
		elemType := m.ElementTypes[elemID]

		// Skip non-3D elements
		if elemType.GetDimension() != 3 {
			continue
		}

		vertices := m.EtoV[elemID]

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
// For higher-order elements, returns only corner nodes of faces
func GetElementFaces(elemType ElementType, vertices []int) [][]int {
	switch elemType {
	case Tet, Tet10:
		// Use corner nodes for all tet variants
		v := vertices
		if elemType == Tet10 {
			v = vertices[:4]
		}
		return [][]int{
			{v[0], v[2], v[1]}, // Face 0
			{v[0], v[1], v[3]}, // Face 1
			{v[1], v[2], v[3]}, // Face 2
			{v[0], v[3], v[2]}, // Face 3
		}

	case Hex, Hex20, Hex27:
		// Use corner nodes for all hex variants
		v := vertices
		if elemType == Hex20 || elemType == Hex27 {
			v = vertices[:8]
		}
		return [][]int{
			{v[0], v[3], v[2], v[1]}, // Face 0 (bottom)
			{v[4], v[5], v[6], v[7]}, // Face 1 (top)
			{v[0], v[1], v[5], v[4]}, // Face 2
			{v[1], v[2], v[6], v[5]}, // Face 3
			{v[2], v[3], v[7], v[6]}, // Face 4
			{v[3], v[0], v[4], v[7]}, // Face 5
		}

	case Prism, Prism15, Prism18:
		// Use corner nodes for all prism variants
		v := vertices
		if elemType == Prism15 || elemType == Prism18 {
			v = vertices[:6]
		}
		return [][]int{
			{v[0], v[2], v[1]},       // Face 0 (bottom tri)
			{v[3], v[4], v[5]},       // Face 1 (top tri)
			{v[0], v[1], v[4], v[3]}, // Face 2 (quad)
			{v[1], v[2], v[5], v[4]}, // Face 3 (quad)
			{v[2], v[0], v[3], v[5]}, // Face 4 (quad)
		}

	case Pyramid, Pyramid13, Pyramid14:
		// Use corner nodes for all pyramid variants
		v := vertices
		if elemType == Pyramid13 || elemType == Pyramid14 {
			v = vertices[:5]
		}
		return [][]int{
			{v[0], v[3], v[2], v[1]}, // Face 0 (base quad)
			{v[0], v[1], v[4]},       // Face 1 (tri)
			{v[1], v[2], v[4]},       // Face 2 (tri)
			{v[2], v[3], v[4]},       // Face 3 (tri)
			{v[3], v[0], v[4]},       // Face 4 (tri)
		}

	default:
		return [][]int{}
	}
}

// PrintStatistics prints mesh statistics
func (m *Mesh) PrintStatistics() {
	fmt.Printf("Mesh Statistics:\n")
	fmt.Printf("  Format version: %s\n", m.FormatVersion)
	fmt.Printf("  Binary: %v\n", m.IsBinary)
	fmt.Printf("  Vertices: %d\n", m.NumVertices)
	fmt.Printf("  Elements: %d\n", m.NumElements)
	fmt.Printf("  Faces: %d\n", m.NumFaces)

	// Count element types by dimension
	dimCounts := make(map[int]int)
	typeCounts := make(map[ElementType]int)
	for _, t := range m.ElementTypes {
		typeCounts[t]++
		dimCounts[t.GetDimension()]++
	}

	fmt.Printf("  Element dimensions:\n")
	for dim := 0; dim <= 3; dim++ {
		if count, ok := dimCounts[dim]; ok && count > 0 {
			fmt.Printf("    %dD: %d\n", dim, count)
		}
	}

	fmt.Printf("  Element types:\n")
	for t, count := range typeCounts {
		if count > 0 {
			fmt.Printf("    %s: %d\n", t, count)
		}
	}

	// Physical groups
	if len(m.ElementGroups) > 0 {
		fmt.Printf("  Physical groups: %d\n", len(m.ElementGroups))
		for tag, group := range m.ElementGroups {
			fmt.Printf("    Tag %d (%s): %d elements\n", tag, group.Name, len(group.Elements))
		}
	}

	// Count boundary faces
	boundaryFaces := 0
	for i := 0; i < m.NumElements; i++ {
		if m.EToE != nil && i < len(m.EToE) {
			for _, neighbor := range m.EToE[i] {
				if neighbor < 0 {
					boundaryFaces++
				}
			}
		}
	}
	fmt.Printf("  Boundary faces: %d\n", boundaryFaces)

	// Periodic boundaries
	if len(m.Periodics) > 0 {
		fmt.Printf("  Periodic boundaries: %d\n", len(m.Periodics))
	}
}
