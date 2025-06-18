package utils

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
