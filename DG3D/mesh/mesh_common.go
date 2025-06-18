package mesh

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
	"sort"
	"strconv"
	"strings"
)

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
	ElementType   utils.ElementType // Type of boundary element (Line, Triangle, Quad)
	Nodes         []int             // Node indices forming the boundary element
	ParentElement int               // Parent volume element (-1 if none)
	ParentFace    int               // Local face ID within parent element (-1 if none)
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
	EtoV         [][]int             // Element to vertex connectivity [nelems][nverts_per_elem]
	ElementTypes []utils.ElementType // Element type for each element
	ElementTags  [][]int             // All tags for each element (physical, elementary, etc.)
	ElementIDMap map[int]int         // Maps original element IDs to array indices

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
func (m *Mesh) AddElement(elemID int, elemType utils.ElementType, tags []int, nodeIDs []int) error {
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
func (m *Mesh) FilterByDimension(dim int) ([]int, [][]int, []utils.ElementType) {
	var indices []int
	var elements [][]int
	var types []utils.ElementType

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

func GetElementFaces(elemType utils.ElementType, vertices []int) [][]int {
	return utils.GetElementFaces(elemType, vertices)
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
	typeCounts := make(map[utils.ElementType]int)
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

// ConvertToMesh converts a CompleteMesh to an actual Mesh structure
func ConvertToMesh(cm utils.CompleteMesh) *Mesh {
	mesh := NewMesh()

	// Add nodes
	for name, idx := range cm.Nodes.NodeMap {
		nodeID := cm.Nodes.NodeIDMap[name]
		coords := cm.Nodes.Nodes[idx]
		mesh.AddNode(nodeID, coords)
	}

	// Add elements
	elemID := 1
	for _, elemSet := range cm.Elements {
		for i, elemNodes := range elemSet.Elements {
			// Convert logical names to node IDs
			nodeIDs := make([]int, len(elemNodes))
			for j, nodeName := range elemNodes {
				nodeIDs[j] = cm.Nodes.NodeIDMap[nodeName]
			}

			// Get properties
			props := utils.ElementProps{}
			if i < len(elemSet.Properties) {
				props = elemSet.Properties[i]
			}

			tags := []int{props.PhysicalTag, props.GeometricTag}
			if props.PartitionTag > 0 {
				tags = append(tags, 1, props.PartitionTag)
			}

			mesh.AddElement(elemID, elemSet.Type, tags, nodeIDs)
			elemID++
		}
	}

	mesh.BuildConnectivity()
	return mesh
}
