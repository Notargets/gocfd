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
