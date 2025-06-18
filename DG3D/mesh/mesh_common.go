package mesh

import (
	"fmt"
	"sort"
	"strconv"
	"strings"
)

// ElementType represents different finite element types
type ElementType int

const (
	Unknown ElementType = iota
	// 0D elements
	Point
	// 1D elements
	Line
	Line3 // 3-node line (quadratic)
	// 2D elements
	Triangle
	Quad
	Triangle6  // 6-node triangle (quadratic)
	Triangle9  // 9-node triangle
	Triangle10 // 10-node triangle
	Quad8      // 8-node quad (quadratic)
	Quad9      // 9-node quad
	// 3D elements
	Tet
	Hex
	Prism
	Pyramid
	Tet10     // 10-node tetrahedron (quadratic)
	Hex20     // 20-node hexahedron (quadratic)
	Hex27     // 27-node hexahedron
	Prism15   // 15-node prism (quadratic)
	Prism18   // 18-node prism
	Pyramid13 // 13-node pyramid
	Pyramid14 // 14-node pyramid
)

// String representation of element types
func (e ElementType) String() string {
	names := []string{
		"Unknown",
		"Point",
		"Line", "Line3",
		"Triangle", "Quad", "Triangle6", "Triangle9", "Triangle10", "Quad8", "Quad9",
		"Tet", "Hex", "Prism", "Pyramid",
		"Tet10", "Hex20", "Hex27", "Prism15", "Prism18", "Pyramid13", "Pyramid14",
	}
	if int(e) < len(names) {
		return names[e]
	}
	return "Invalid"
}

// GetDimension returns the spatial dimension of the element
func (e ElementType) GetDimension() int {
	switch e {
	case Point:
		return 0
	case Line, Line3:
		return 1
	case Triangle, Quad, Triangle6, Triangle9, Triangle10, Quad8, Quad9:
		return 2
	case Tet, Hex, Prism, Pyramid, Tet10, Hex20, Hex27, Prism15, Prism18, Pyramid13, Pyramid14:
		return 3
	default:
		return -1
	}
}

// GetNumNodes returns the number of nodes for each element type
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
	case Quad:
		return 4
	case Triangle6:
		return 6
	case Triangle9:
		return 9
	case Triangle10:
		return 10
	case Quad8:
		return 8
	case Quad9:
		return 9
	case Tet:
		return 4
	case Hex:
		return 8
	case Prism:
		return 6
	case Pyramid:
		return 5
	case Tet10:
		return 10
	case Hex20:
		return 20
	case Hex27:
		return 27
	case Prism15:
		return 15
	case Prism18:
		return 18
	case Pyramid13:
		return 13
	case Pyramid14:
		return 14
	default:
		return 0
	}
}

// GetNumFaces returns the number of faces for 3D elements
func (e ElementType) GetNumFaces() int {
	switch e {
	case Tet, Tet10:
		return 4
	case Hex, Hex20, Hex27:
		return 6
	case Prism, Prism15, Prism18:
		return 5
	case Pyramid, Pyramid13, Pyramid14:
		return 5
	default:
		return 0
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
	Dimension  int
	Tag        int
	Name       string
	Elements   []int // Element indices in this group
	MaterialID int   // Material identifier (for Gambit)
	Flags      []int // Additional flags (for Gambit NFLAGS)
}

// NodeGroup represents a group of nodes (for boundary conditions, etc.)
type NodeGroup struct {
	Tag   int
	Name  string
	Nodes []int // Node indices in this group
}

// BoundaryElement represents a boundary element/face
type BoundaryElement struct {
	ElementType   ElementType // Type of boundary element (Line, Triangle, Quad)
	Nodes         []int       // Node indices forming the boundary element
	ParentElement int         // Parent volume element (-1 if none)
	ParentFace    int         // Local face ID within parent element (-1 if none)
}

// Entity represents a geometric entity (point, curve, surface, volume)
type Entity struct {
	Dimension        int
	Tag              int
	BoundingBox      [2][3]float64
	PhysicalTags     []int
	BoundingEntities []int // Lower-dimension entities bounding this one
}

// GhostElement represents ghost element data for parallel meshes
type GhostElement struct {
	ElementTag      int
	OwnerPartition  int
	GhostPartitions []int
}

// NonConformalFace represents non-conformal mesh connections
type NonConformalFace struct {
	MasterElement int
	MasterFace    int
	SlaveElements []int
	SlaveFaces    []int
}

// Periodic represents periodic boundary condition data
type Periodic struct {
	Dimension       int
	SlaveTag        int
	MasterTag       int
	NodeMap         map[int]int // slave node -> master node
	AffineTransform []float64   // 4x4 matrix (3D) or 3x3 matrix (2D) in row-major order
}

// Mesh represents a complete unstructured mesh with all connectivity
type Mesh struct {
	// ===== Node Data =====
	Vertices         [][]float64 // Vertex coordinates [nvertices][3]
	ParametricCoords [][]float64 // Parametric coordinates [nvertices][up to 3] (u,v,w)
	HasParametric    []bool      // Whether node has parametric coordinates
	NodeIDMap        map[int]int // Maps original node IDs to array indices
	NodeArrayMap     map[int]int // Maps array indices to original node IDs

	// ===== Element Connectivity =====
	EtoV         [][]int       // Element to vertex connectivity [nelems][nverts_per_elem]
	ElementTypes []ElementType // Element type for each element
	ElementTags  [][]int       // All tags for each element (physical, elementary, etc.)
	ElementIDMap map[int]int   // Maps original element IDs to array indices

	// ===== Groups and Sets =====
	ElementGroups map[int]*ElementGroup // Physical/elementary groups by tag
	NodeGroups    map[int]*NodeGroup    // Node groups by tag

	// ===== Face Connectivity =====
	EToE    [][]int        // Element to element connectivity [nelems][nfaces_per_elem]
	EToF    [][]int        // Element to face connectivity [nelems][nfaces_per_elem]
	Faces   []Face         // All unique faces in mesh
	FaceMap map[string]int // Map from sorted vertex string to face ID

	// ===== Boundary Conditions =====
	BoundaryTags     map[int]string               // Boundary condition tag names
	BoundaryElements map[string][]BoundaryElement // tag name -> boundary elements
	Periodics        []Periodic                   // Periodic boundary conditions

	// ===== Geometric Model =====
	Entities map[int]*Entity // Geometric entities by tag

	// ===== Parallel/Partitioning Data =====
	EToP              []int              // Element to partition mapping
	GhostElements     []GhostElement     // Ghost element data
	NonConformalFaces []NonConformalFace // Non-conformal interface data

	// ===== Mesh Statistics =====
	NumElements int
	NumVertices int
	NumFaces    int

	// ===== Format-Specific Data =====
	FormatVersion string // e.g., "2.2", "4.1"
	IsBinary      bool
	DataSize      int                    // Size of data fields (for binary)
	Extensions    map[string]interface{} // Format-specific extensions
}

// NewMesh creates a new mesh and initializes maps
func NewMesh() *Mesh {
	return &Mesh{
		NodeIDMap:        make(map[int]int),
		NodeArrayMap:     make(map[int]int),
		ElementIDMap:     make(map[int]int),
		ElementGroups:    make(map[int]*ElementGroup),
		NodeGroups:       make(map[int]*NodeGroup),
		FaceMap:          make(map[string]int),
		BoundaryTags:     make(map[int]string),
		BoundaryElements: make(map[string][]BoundaryElement),
		Entities:         make(map[int]*Entity),
		Extensions:       make(map[string]interface{}),
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

// AddNodeWithParametric adds a node with parametric coordinates
func (m *Mesh) AddNodeWithParametric(nodeID int, coords []float64, paramCoords []float64) {
	m.AddNode(nodeID, coords)

	// Ensure parametric arrays are sized correctly
	for len(m.ParametricCoords) < len(m.Vertices) {
		m.ParametricCoords = append(m.ParametricCoords, nil)
		m.HasParametric = append(m.HasParametric, false)
	}

	idx := len(m.Vertices) - 1
	m.ParametricCoords[idx] = paramCoords
	m.HasParametric[idx] = true
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

// AddBoundaryElement adds a boundary element to the mesh
func (m *Mesh) AddBoundaryElement(tagName string, elem BoundaryElement) {
	m.BoundaryElements[tagName] = append(m.BoundaryElements[tagName], elem)
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
//
// For each element and face, this function populates:
// - EToE[elem][face] = neighbor element ID (or -1 for boundary)
// - EToF[elem][face] = neighbor's LOCAL face index (or -1 for boundary)
//
// This creates reciprocal connectivity where:
// - If element A's face i connects to element B's face j
// - Then element B's face j connects to element A's face i
//
// CRITICAL: EToF stores the neighbor's LOCAL face index, not global face IDs
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
			// Create sorted vertex key
			sortedVerts := make([]int, len(faceVerts))
			copy(sortedVerts, faceVerts)
			sort.Ints(sortedVerts)

			// Create string key
			parts := make([]string, len(sortedVerts))
			for i, v := range sortedVerts {
				parts[i] = strconv.Itoa(v)
			}
			key := strings.Join(parts, ",")

			// Check if face already exists
			if faceID, exists := m.FaceMap[key]; exists {
				// Face already exists, find the other element
				face := &m.Faces[faceID]
				otherElem := face.Element
				otherLocalID := face.LocalID

				// Set connectivity - CRITICAL FIX HERE
				// EToE stores the neighbor element ID
				// EToF stores the neighbor's LOCAL face ID (not global face ID)
				m.EToE[elemID][localFaceID] = otherElem
				m.EToF[elemID][localFaceID] = otherLocalID // Store neighbor's local face ID

				// Set reverse connectivity
				m.EToE[otherElem][otherLocalID] = elemID
				m.EToF[otherElem][otherLocalID] = localFaceID // Store this element's local face ID
			} else {
				// New face
				faceID := len(m.Faces)
				m.Faces = append(m.Faces, Face{
					Vertices: sortedVerts,
					Element:  elemID,
					LocalID:  localFaceID,
				})
				m.FaceMap[key] = faceID
				// Note: EToE and EToF remain -1 for boundary faces
			}
		}
	}

	m.NumFaces = len(m.Faces)
}

// Required imports for this function:
// import (
//     "sort"
//     "strconv"
//     "strings"
// )
// GetElementFaces returns the faces of an element as vertex lists
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
			{v[0], v[3], v[2]}, // Face 2
			{v[1], v[2], v[3]}, // Face 3
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

// GetMeshDimension returns the highest dimension of elements in the mesh
func (m *Mesh) GetMeshDimension() int {
	maxDim := 0
	for _, etype := range m.ElementTypes {
		dim := etype.GetDimension()
		if dim > maxDim {
			maxDim = dim
		}
	}
	return maxDim
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

	// Boundary elements
	if len(m.BoundaryElements) > 0 {
		fmt.Printf("  Boundary conditions:\n")
		for tag, elems := range m.BoundaryElements {
			fmt.Printf("    %s: %d boundary elements\n", tag, len(elems))
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

	// Entities
	if len(m.Entities) > 0 {
		entCounts := make(map[int]int)
		for _, ent := range m.Entities {
			entCounts[ent.Dimension]++
		}
		fmt.Printf("  Geometric entities:\n")
		for dim := 0; dim <= 3; dim++ {
			if count, ok := entCounts[dim]; ok && count > 0 {
				fmt.Printf("    %dD: %d\n", dim, count)
			}
		}
	}
}
