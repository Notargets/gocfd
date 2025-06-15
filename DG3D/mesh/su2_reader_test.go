package mesh

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// SU2TestBuilder helps build SU2 format test files
type SU2TestBuilder struct {
	tm *TestMeshes
}

// NewSU2TestBuilder creates a new builder with standard test meshes
func NewSU2TestBuilder() *SU2TestBuilder {
	return &SU2TestBuilder{
		tm: GetStandardTestMeshes(),
	}
}

// Helper function to create temporary test files
func createTempSU2File(t *testing.T, content string) string {
	t.Helper()
	tmpFile := filepath.Join(t.TempDir(), "test.su2")
	if err := os.WriteFile(tmpFile, []byte(content), 0644); err != nil {
		t.Fatalf("Failed to create temp file: %v", err)
	}
	return tmpFile
}

// TestReadSU2Dimension tests reading the dimension section
func TestReadSU2Dimension(t *testing.T) {
	testCases := []struct {
		name     string
		content  string
		expected int
		errMsg   string
	}{
		{
			name: "2D mesh",
			content: `NDIME= 2
NPOIN= 3
0.0 0.0
1.0 0.0
0.0 1.0
NELEM= 0`,
			expected: 2,
		},
		{
			name: "3D mesh",
			content: `NDIME= 3
NPOIN= 4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
NELEM= 0`,
			expected: 3,
		},
		{
			name: "Invalid dimension",
			content: `NDIME= 4
NPOIN= 0`,
			errMsg: "unsupported dimension",
		},
		{
			name:    "Missing dimension",
			content: `NPOIN= 0`,
			errMsg:  "NDIME",
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			tmpFile := createTempSU2File(t, tc.content)
			defer os.Remove(tmpFile)

			mesh, err := ReadSU2(tmpFile)
			if tc.errMsg != "" {
				if err == nil {
					t.Error("Expected error, got nil")
				} else if !strings.Contains(err.Error(), tc.errMsg) {
					t.Errorf("Expected error containing '%s', got '%v'", tc.errMsg, err)
				}
			} else {
				if err != nil {
					t.Fatalf("Unexpected error: %v", err)
				}
				// Verify mesh was read successfully
				if mesh == nil {
					t.Fatal("Expected valid mesh, got nil")
				}
			}
		})
	}
}

// TestReadSU2StandardMeshes tests reading standard test meshes
func TestReadSU2StandardMeshes(t *testing.T) {
	builder := NewSU2TestBuilder()

	t.Run("SingleTet", func(t *testing.T) {
		content := builder.BuildSingleTetTest()
		tmpFile := createTempSU2File(t, content)
		defer os.Remove(tmpFile)

		mesh, err := ReadSU2(tmpFile)
		if err != nil {
			t.Fatalf("Failed to read SU2 file: %v", err)
		}

		// Verify vertices
		if mesh.NumVertices != 4 {
			t.Errorf("Expected 4 vertices, got %d", mesh.NumVertices)
		}

		// Verify element
		if mesh.NumElements != 1 {
			t.Errorf("Expected 1 element, got %d", mesh.NumElements)
		}
		if mesh.ElementTypes[0] != Tet {
			t.Errorf("Expected Tet element type, got %v", mesh.ElementTypes[0])
		}
		if len(mesh.EtoV[0]) != 4 {
			t.Errorf("Expected 4 nodes for tet, got %d", len(mesh.EtoV[0]))
		}
	})

	t.Run("MixedMesh", func(t *testing.T) {
		content := builder.BuildMixedElementTest()
		tmpFile := createTempSU2File(t, content)
		defer os.Remove(tmpFile)

		mesh, err := ReadSU2(tmpFile)
		if err != nil {
			t.Fatalf("Failed to read SU2 file: %v", err)
		}

		// Get expected element counts from test mesh
		tm := GetStandardTestMeshes()
		var expectedTypes []ElementType
		for _, elemSet := range tm.MixedMesh.Elements {
			for range elemSet.Elements {
				expectedTypes = append(expectedTypes, elemSet.Type)
			}
		}

		// Verify element count
		if mesh.NumElements != len(expectedTypes) {
			t.Errorf("Expected %d elements, got %d", len(expectedTypes), mesh.NumElements)
		}

		// Verify element types
		for i, expectedType := range expectedTypes {
			if i >= len(mesh.ElementTypes) {
				t.Errorf("Missing element %d", i)
				continue
			}
			if mesh.ElementTypes[i] != expectedType {
				t.Errorf("Element %d: expected type %v, got %v",
					i, expectedType, mesh.ElementTypes[i])
			}
		}
	})

	t.Run("2DMesh", func(t *testing.T) {
		content := builder.Build2DTriangleTest()
		tmpFile := createTempSU2File(t, content)
		defer os.Remove(tmpFile)

		mesh, err := ReadSU2(tmpFile)
		if err != nil {
			t.Fatalf("Failed to read SU2 file: %v", err)
		}

		// Verify 2D elements
		for i, elemType := range mesh.ElementTypes {
			if elemType.GetDimension() != 2 {
				t.Errorf("Element %d: expected 2D element, got %dD (%v)",
					i, elemType.GetDimension(), elemType)
			}
		}
	})
}

// TestReadSU2ElementTypes tests all supported element types
func TestReadSU2ElementTypes(t *testing.T) {
	testCases := []struct {
		name         string
		su2Type      int
		expectedType ElementType
		numNodes     int
		dimension    int
	}{
		// 2D elements
		{"Line", 3, Line, 2, 2},
		{"Triangle", 5, Triangle, 3, 2},
		{"Quadrilateral", 9, Quad, 4, 2},
		// 3D elements
		{"Tetrahedron", 10, Tet, 4, 3},
		{"Hexahedron", 12, Hex, 8, 3},
		{"Prism", 13, Prism, 6, 3},
		{"Pyramid", 14, Pyramid, 5, 3},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Create minimal file with one element of the specified type
			content := fmt.Sprintf(`NDIME= %d
NPOIN= %d`, tc.dimension, tc.numNodes)

			// Add nodes (generic positions)
			for i := 0; i < tc.numNodes; i++ {
				content += fmt.Sprintf("\n%e %e", float64(i)*0.1, float64(i)*0.2)
				if tc.dimension == 3 {
					content += fmt.Sprintf(" %e", float64(i)*0.3)
				}
			}

			// Add element
			content += fmt.Sprintf("\nNELEM= 1\n%d", tc.su2Type)
			for i := 0; i < tc.numNodes; i++ {
				content += fmt.Sprintf(" %d", i)
			}

			tmpFile := createTempSU2File(t, content)
			defer os.Remove(tmpFile)

			mesh, err := ReadSU2(tmpFile)
			if err != nil {
				t.Fatalf("Failed to read SU2 file: %v", err)
			}

			if mesh.NumElements != 1 {
				t.Errorf("Expected 1 element, got %d", mesh.NumElements)
			}
			if mesh.ElementTypes[0] != tc.expectedType {
				t.Errorf("Expected element type %v, got %v",
					tc.expectedType, mesh.ElementTypes[0])
			}
			if len(mesh.EtoV[0]) != tc.numNodes {
				t.Errorf("Expected %d nodes, got %d",
					tc.numNodes, len(mesh.EtoV[0]))
			}
		})
	}
}

// TestReadSU2BoundaryMarkers tests reading boundary markers
func TestReadSU2BoundaryMarkers(t *testing.T) {
	content := `NDIME= 3
NPOIN= 4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
NELEM= 1
10 0 1 2 3
NMARK= 4
MARKER_TAG= lower
MARKER_ELEMS= 1
5 0 1 2
MARKER_TAG= right
MARKER_ELEMS= 1
5 1 2 3
MARKER_TAG= upper
MARKER_ELEMS= 1
5 2 3 0
MARKER_TAG= left
MARKER_ELEMS= 1
5 3 0 1`

	tmpFile := createTempSU2File(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadSU2(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read SU2 file: %v", err)
	}

	// Check boundary tags
	expectedTags := []string{"lower", "right", "upper", "left"}
	if len(mesh.BoundaryTags) != len(expectedTags) {
		t.Errorf("Expected %d boundary tags, got %d",
			len(expectedTags), len(mesh.BoundaryTags))
	}

	// Check that all expected tags are present
	foundTags := make(map[string]bool)
	for _, tag := range mesh.BoundaryTags {
		foundTags[tag] = true
	}

	for _, expected := range expectedTags {
		if !foundTags[expected] {
			t.Errorf("Missing expected boundary tag: %s", expected)
		}
	}
}

// TestReadSU2NodeOrdering tests that nodes maintain correct ordering
func TestReadSU2NodeOrdering(t *testing.T) {
	// Test that nodes are stored with 0-based indices matching their order
	content := `NDIME= 3
NPOIN= 4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
NELEM= 1
10 0 1 2 3`

	tmpFile := createTempSU2File(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadSU2(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read SU2 file: %v", err)
	}

	// Check node coordinates
	expectedCoords := [][]float64{
		{0.0, 0.0, 0.0},
		{1.0, 0.0, 0.0},
		{0.0, 1.0, 0.0},
		{0.0, 0.0, 1.0},
	}

	for i, expected := range expectedCoords {
		if i >= len(mesh.Vertices) {
			t.Errorf("Missing vertex %d", i)
			continue
		}
		for j, coord := range expected {
			if mesh.Vertices[i][j] != coord {
				t.Errorf("Vertex %d coord %d: expected %f, got %f",
					i, j, coord, mesh.Vertices[i][j])
			}
		}
	}

	// Check element connectivity
	if len(mesh.EtoV[0]) != 4 {
		t.Fatalf("Expected 4 nodes in element, got %d", len(mesh.EtoV[0]))
	}
	for i := 0; i < 4; i++ {
		if mesh.EtoV[0][i] != i {
			t.Errorf("Element node %d: expected %d, got %d",
				i, i, mesh.EtoV[0][i])
		}
	}
}

// TestReadSU2ErrorHandling tests various error conditions
func TestReadSU2ErrorHandling(t *testing.T) {
	testCases := []struct {
		name    string
		content string
		errMsg  string
	}{
		{
			name:    "NonexistentFile",
			content: "",
			errMsg:  "no such file",
		},
		{
			name: "InvalidElementType",
			content: `NDIME= 3
NPOIN= 4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
NELEM= 1
999 0 1 2 3`,
			errMsg: "unknown element type",
		},
		{
			name: "TooFewNodes",
			content: `NDIME= 3
NPOIN= 3
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
NELEM= 1
10 0 1 2 3`,
			errMsg: "", // Should handle gracefully or error
		},
		{
			name: "InvalidNodeIndex",
			content: `NDIME= 3
NPOIN= 4
0.0 0.0 0.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
NELEM= 1
10 0 1 2 10`,
			errMsg: "", // Should handle invalid indices
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			var tmpFile string
			if tc.content != "" {
				tmpFile = createTempSU2File(t, tc.content)
				defer os.Remove(tmpFile)
			} else {
				tmpFile = "nonexistent.su2"
			}

			_, err := ReadSU2(tmpFile)
			if tc.errMsg != "" {
				if err == nil {
					t.Error("Expected error, got nil")
				} else if !strings.Contains(err.Error(), tc.errMsg) {
					t.Errorf("Expected error containing '%s', got '%v'", tc.errMsg, err)
				}
			}
		})
	}
}

// TestReadSU2SquareExample tests the documented square mesh example from SU2 docs
func TestReadSU2SquareExample(t *testing.T) {
	// This is the exact example from the SU2 documentation
	content := `NDIME= 2
NPOIN= 9
0.00000000000000 0.00000000000000
0.50000000000000 0.00000000000000
1.00000000000000 0.00000000000000
0.00000000000000 0.50000000000000
0.50000000000000 0.50000000000000
1.00000000000000 0.50000000000000
0.00000000000000 1.00000000000000
0.50000000000000 1.00000000000000
1.00000000000000 1.00000000000000
NELEM= 8
5 0 1 3
5 1 4 3
5 1 2 4
5 2 5 4
5 3 4 6
5 4 7 6
5 4 5 7
5 5 8 7
NMARK= 4
MARKER_TAG= lower
MARKER_ELEMS= 2
3 0 1
3 1 2
MARKER_TAG= right
MARKER_ELEMS= 2
3 2 5
3 5 8
MARKER_TAG= upper
MARKER_ELEMS= 2
3 8 7
3 7 6
MARKER_TAG= left
MARKER_ELEMS= 2
3 6 3
3 3 0`

	tmpFile := createTempSU2File(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadSU2(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read SU2 file: %v", err)
	}

	// Verify the documented example
	if mesh.NumVertices != 9 {
		t.Errorf("Expected 9 vertices, got %d", mesh.NumVertices)
	}
	if mesh.NumElements != 8 {
		t.Errorf("Expected 8 elements, got %d", mesh.NumElements)
	}

	// All elements should be triangles
	for i, elemType := range mesh.ElementTypes {
		if elemType != Triangle {
			t.Errorf("Element %d: expected Triangle, got %v", i, elemType)
		}
	}

	// Check boundary markers
	if len(mesh.BoundaryTags) != 4 {
		t.Errorf("Expected 4 boundary markers, got %d", len(mesh.BoundaryTags))
	}
}

// TestReadSU2LegacyFormat tests handling of legacy format with explicit indices
func TestReadSU2LegacyFormat(t *testing.T) {
	// Legacy format includes explicit indices at end of lines (should be ignored)
	content := `NDIME= 2
NPOIN= 4
0.0 0.0 0
1.0 0.0 1
1.0 1.0 2
0.0 1.0 3
NELEM= 2
5 0 1 2 0
5 0 2 3 1`

	tmpFile := createTempSU2File(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadSU2(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read SU2 file: %v", err)
	}

	// Should still read correctly, ignoring the extra indices
	if mesh.NumVertices != 4 {
		t.Errorf("Expected 4 vertices, got %d", mesh.NumVertices)
	}
	if mesh.NumElements != 2 {
		t.Errorf("Expected 2 elements, got %d", mesh.NumElements)
	}
}

// Builder methods

func (b *SU2TestBuilder) BuildSingleTetTest() string {
	tm := b.tm
	mesh := CompleteMesh{
		Nodes:       tm.TetraNodes,
		Elements:    []ElementSet{tm.SingleTet},
		Dimension:   3,
		BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 1}},
	}
	filtered := filterToUsedNodesSU2(&mesh)
	return b.BuildFromCompleteMesh(&filtered)
}

func (b *SU2TestBuilder) BuildMixedElementTest() string {
	tm := b.tm
	filtered := filterToUsedNodesSU2(&tm.MixedMesh)
	return b.BuildFromCompleteMesh(&filtered)
}

func (b *SU2TestBuilder) Build2DTriangleTest() string {
	// Create a simple 2D triangle mesh
	nodes := NodeSet{
		Nodes: [][]float64{
			{0, 0},
			{1, 0},
			{0, 1},
		},
		NodeMap: map[string]int{
			"v0": 0, "v1": 1, "v2": 2,
		},
	}

	nodes.NodeIDMap = make(map[string]int)
	for name, idx := range nodes.NodeMap {
		nodes.NodeIDMap[name] = idx + 1
	}

	elements := []ElementSet{
		{
			Type: Triangle,
			Elements: [][]string{
				{"v0", "v1", "v2"},
			},
		},
	}

	mesh := CompleteMesh{
		Nodes:       nodes,
		Elements:    elements,
		Dimension:   2,
		BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 0}},
	}

	return b.BuildFromCompleteMesh(&mesh)
}

func (b *SU2TestBuilder) BuildFromCompleteMesh(mesh *CompleteMesh) string {
	var lines []string

	// Dimension
	lines = append(lines, fmt.Sprintf("NDIME= %d", mesh.Dimension))

	// Nodes
	lines = append(lines, fmt.Sprintf("NPOIN= %d", len(mesh.Nodes.Nodes)))
	for _, node := range mesh.Nodes.Nodes {
		if mesh.Dimension == 2 {
			lines = append(lines, fmt.Sprintf("%e %e", node[0], node[1]))
		} else {
			lines = append(lines, fmt.Sprintf("%e %e %e", node[0], node[1], node[2]))
		}
	}

	// Elements
	totalElements := 0
	for _, elemSet := range mesh.Elements {
		totalElements += len(elemSet.Elements)
	}

	lines = append(lines, fmt.Sprintf("NELEM= %d", totalElements))

	for _, elemSet := range mesh.Elements {
		su2Type := elementTypeToSU2(elemSet.Type)
		for _, elem := range elemSet.Elements {
			line := fmt.Sprintf("%d", su2Type)
			// Convert logical node names to indices
			for _, nodeName := range elem {
				nodeIdx := mesh.Nodes.NodeMap[nodeName]
				line += fmt.Sprintf(" %d", nodeIdx)
			}
			lines = append(lines, line)
		}
	}

	// Add some boundary markers for testing
	if mesh.Dimension == 2 {
		lines = append(lines, b.add2DBoundaryMarkers(mesh)...)
	} else {
		lines = append(lines, b.add3DBoundaryMarkers(mesh)...)
	}

	return strings.Join(lines, "\n")
}

func (b *SU2TestBuilder) add2DBoundaryMarkers(mesh *CompleteMesh) []string {
	var lines []string

	// For 2D meshes, add line boundary elements
	lines = append(lines, "NMARK= 2")
	lines = append(lines, "MARKER_TAG= boundary1")
	lines = append(lines, "MARKER_ELEMS= 1")
	lines = append(lines, "3 0 1") // Line element
	lines = append(lines, "MARKER_TAG= boundary2")
	lines = append(lines, "MARKER_ELEMS= 1")
	lines = append(lines, "3 1 2") // Line element

	return lines
}

func (b *SU2TestBuilder) add3DBoundaryMarkers(mesh *CompleteMesh) []string {
	var lines []string

	// For 3D meshes, add triangular boundary elements
	lines = append(lines, "NMARK= 2")
	lines = append(lines, "MARKER_TAG= wall")
	lines = append(lines, "MARKER_ELEMS= 1")
	lines = append(lines, "5 0 1 2") // Triangle element
	lines = append(lines, "MARKER_TAG= inlet")
	lines = append(lines, "MARKER_ELEMS= 1")
	lines = append(lines, "5 0 2 3") // Triangle element

	return lines
}

// elementTypeToSU2 converts internal element types to SU2/VTK type codes
func elementTypeToSU2(et ElementType) int {
	switch et {
	case Line:
		return 3
	case Triangle:
		return 5
	case Quad:
		return 9
	case Tet:
		return 10
	case Hex:
		return 12
	case Prism:
		return 13
	case Pyramid:
		return 14
	default:
		return 0
	}
}

// Helper function to filter a CompleteMesh to only include used nodes
func filterToUsedNodesSU2(mesh *CompleteMesh) CompleteMesh {
	// Collect all nodes used by elements
	usedNodes := make(map[string]bool)
	for _, elemSet := range mesh.Elements {
		for _, elem := range elemSet.Elements {
			for _, nodeName := range elem {
				usedNodes[nodeName] = true
			}
		}
	}

	// Build filtered node set
	var nodes [][]float64
	nodeMap := make(map[string]int)
	nodeIDMap := make(map[string]int)
	idx := 0

	// Create ordered list of used nodes
	var orderedNames []string
	for name := range mesh.Nodes.NodeMap {
		if usedNodes[name] {
			orderedNames = append(orderedNames, name)
		}
	}

	// Sort for consistent ordering
	sort.Strings(orderedNames)

	// Build new node arrays
	for _, name := range orderedNames {
		origIdx := mesh.Nodes.NodeMap[name]
		nodes = append(nodes, mesh.Nodes.Nodes[origIdx])
		nodeMap[name] = idx
		nodeIDMap[name] = idx + 1
		idx++
	}

	filteredNodes := NodeSet{
		Nodes:     nodes,
		NodeMap:   nodeMap,
		NodeIDMap: nodeIDMap,
	}

	return CompleteMesh{
		Nodes:       filteredNodes,
		Elements:    mesh.Elements,
		Dimension:   mesh.Dimension,
		BoundingBox: mesh.BoundingBox,
	}
}
