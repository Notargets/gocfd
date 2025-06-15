package readers

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// GambitTestBuilder helps build Gambit neutral format test files
type GambitTestBuilder struct {
	tm *mesh.TestMeshes
}

// NewGambitTestBuilder creates a new builder with standard test meshes
func NewGambitTestBuilder() *GambitTestBuilder {
	return &GambitTestBuilder{
		tm: mesh.GetStandardTestMeshes(),
	}
}

// Helper function to create temporary test files
func createTempNeuFile(t *testing.T, content string) string {
	t.Helper()
	tmpFile := filepath.Join(t.TempDir(), "test.neu")
	if err := os.WriteFile(tmpFile, []byte(content), 0644); err != nil {
		t.Fatalf("Failed to create temp file: %v", err)
	}
	return tmpFile
}

// TestReadGambitNeutralHeader tests reading the header section
func TestReadGambitNeutralHeader(t *testing.T) {
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test mesh for unit testing
PROGRAM:                  Gmsh     VERSION:  4.13.1
Sat Jun  7 21:41:35 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         8        1         1         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   1.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         5   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         6   1.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         7   1.00000000000e+00   1.00000000000e+00   1.00000000000e+00
         8   0.00000000000e+00   1.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         4         8         1         2         3         4         5         6         7         8
ENDOFSECTION
       ELEMENT GROUP 2.0.0
GROUP:           1 ELEMENTS:           1 MATERIAL:           2 NFLAGS:           0
fluid
         1
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGambitNeutral(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gambit neutral file: %v", err)
	}

	// Verify problem size parameters
	if mesh.NumVertices != 8 {
		t.Errorf("Expected 8 vertices, got %d", mesh.NumVertices)
	}
	if mesh.NumElements != 1 {
		t.Errorf("Expected 1 element, got %d", mesh.NumElements)
	}
}

// TestReadGambitNeutralStandardMeshes tests reading standard test meshes
func TestReadGambitNeutralStandardMeshes(t *testing.T) {
	builder := NewGambitTestBuilder()

	t.Run("SingleTet", func(t *testing.T) {
		content := builder.BuildSingleTetTest()
		tmpFile := createTempNeuFile(t, content)
		defer os.Remove(tmpFile)

		msh, err := ReadGambitNeutral(tmpFile)
		if err != nil {
			t.Fatalf("Failed to read Gambit neutral file: %v", err)
		}

		// Verify vertices
		if msh.NumVertices != 4 {
			t.Errorf("Expected 4 vertices, got %d", msh.NumVertices)
		}

		// Verify element
		if msh.NumElements != 1 {
			t.Errorf("Expected 1 element, got %d", msh.NumElements)
		}
		if msh.ElementTypes[0] != mesh.Tet {
			t.Errorf("Expected Tet element type, got %v", msh.ElementTypes[0])
		}
		if len(msh.EtoV[0]) != 4 {
			t.Errorf("Expected 4 nodes for tet, got %d", len(msh.EtoV[0]))
		}
	})

	t.Run("MixedMesh", func(t *testing.T) {
		content := builder.BuildMixedElementTest()
		tmpFile := createTempNeuFile(t, content)
		defer os.Remove(tmpFile)

		msh, err := ReadGambitNeutral(tmpFile)
		if err != nil {
			t.Fatalf("Failed to read Gambit neutral file: %v", err)
		}

		// Get expected element counts
		tm := mesh.GetStandardTestMeshes()
		var expectedTypes []mesh.ElementType
		for _, elemSet := range tm.MixedMesh.Elements {
			for range elemSet.Elements {
				expectedTypes = append(expectedTypes, elemSet.Type)
			}
		}

		// Verify element count
		if msh.NumElements != len(expectedTypes) {
			t.Errorf("Expected %d elements, got %d", len(expectedTypes), msh.NumElements)
		}

		// Verify element types
		for i, expectedType := range expectedTypes {
			if i >= len(msh.ElementTypes) {
				t.Errorf("Missing element %d", i)
				continue
			}
			if msh.ElementTypes[i] != expectedType {
				t.Errorf("Element %d: expected type %v, got %v",
					i, expectedType, msh.ElementTypes[i])
			}
		}
	})
}

// TestReadGambitNeutralElementTypes tests all supported 3D element types
func TestReadGambitNeutralElementTypes(t *testing.T) {
	testCases := []struct {
		name         string
		gambitType   int
		expectedType mesh.ElementType
		numNodes     int
	}{
		{"Tetrahedron", 6, mesh.Tet, 4},
		{"Hexahedron", 4, mesh.Hex, 8},
		{"Wedge/Prism", 5, mesh.Prism, 6},
		{"Pyramid", 7, mesh.Pyramid, 5},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Create minimal file with one element of the specified type
			content := fmt.Sprintf(`        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test %s element
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
        %d         1         1         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0`, tc.name, tc.numNodes)

			// Add nodes (generic positions)
			for i := 1; i <= tc.numNodes; i++ {
				content += fmt.Sprintf("\n         %d   %e   %e   %e",
					i, float64(i)*0.1, float64(i)*0.2, float64(i)*0.3)
			}
			content += "\nENDOFSECTION\n   ELEMENTS/CELLS 2.0.0\n"

			// Add element
			content += fmt.Sprintf("         1         %d         %d", tc.gambitType, tc.numNodes)
			for i := 1; i <= tc.numNodes; i++ {
				content += fmt.Sprintf("         %d", i)
			}
			content += "\nENDOFSECTION"

			tmpFile := createTempNeuFile(t, content)
			defer os.Remove(tmpFile)

			mesh, err := ReadGambitNeutral(tmpFile)
			if err != nil {
				t.Fatalf("Failed to read Gambit neutral file: %v", err)
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

// TestReadGambitNeutral2DElementFiltering tests that 2D elements are filtered out
func TestReadGambitNeutral2DElementFiltering(t *testing.T) {
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test 2D element filtering
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         8         5         1         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   1.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         5   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         6   1.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         7   1.00000000000e+00   1.00000000000e+00   1.00000000000e+00
         8   0.00000000000e+00   1.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         1         2         1         2
         2         2         4         1         2         3         4
         3         3         3         1         2         3
         4         4         8         1         2         3         4         5         6         7         8
         5         6         4         1         2         3         5
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	msh, err := ReadGambitNeutral(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gambit neutral file: %v", err)
	}
	// Debug output
	t.Logf("Total elements read: %d\n", msh.NumElements)
	for i, elemType := range msh.ElementTypes {
		t.Logf("Element %d: type=%v, dimension=%d\n", i, elemType,
			elemType.GetDimension())
	}

	// Count 3D elements
	count3D := 0
	for _, elemType := range msh.ElementTypes {
		if elemType.GetDimension() == 3 {
			count3D++
		}
	}
	fmt.Printf("3D elements found: %d\n", count3D)
	// Should only have 2 3D elements (hex and tet)
	if count3D != 2 {
		t.Errorf("Expected 2 3D elements (filtering out 2D), got %d", count3D)
	}

	// Check that we have the right element types
	hasHex := false
	hasTet := false
	for _, elemType := range msh.ElementTypes {
		if elemType == mesh.Hex {
			hasHex = true
		} else if elemType == mesh.Tet {
			hasTet = true
		}
	}

	if !hasHex {
		t.Error("Expected to find Hex element")
	}
	if !hasTet {
		t.Error("Expected to find Tet element")
	}
}

// TestReadGambitNeutralElementGroups tests reading element group information
func TestReadGambitNeutralElementGroups(t *testing.T) {
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test element groups
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         8         2         2         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   1.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         5   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         6   1.00000000000e+00   0.00000000000e+00   1.00000000000e+00
         7   1.00000000000e+00   1.00000000000e+00   1.00000000000e+00
         8   0.00000000000e+00   1.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         6         4         1         2         3         5
         2         6         4         2         3         4         6
ENDOFSECTION
       ELEMENT GROUP 2.0.0
GROUP:           1 ELEMENTS:           1 MATERIAL:           2 NFLAGS:           0
fluid
         1
GROUP:           2 ELEMENTS:           1 MATERIAL:           3 NFLAGS:           0
solid
         2
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGambitNeutral(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gambit neutral file: %v", err)
	}

	// Check element tags
	if len(mesh.ElementTags) != 2 {
		t.Fatalf("Expected 2 element tag entries, got %d", len(mesh.ElementTags))
	}

	// Element 0 should have group 1
	if len(mesh.ElementTags[0]) == 0 || mesh.ElementTags[0][0] != 1 {
		t.Errorf("Element 0: expected tag [1], got %v", mesh.ElementTags[0])
	}

	// Element 1 should have group 2
	if len(mesh.ElementTags[1]) == 0 || mesh.ElementTags[1][0] != 2 {
		t.Errorf("Element 1: expected tag [2], got %v", mesh.ElementTags[1])
	}
}

// TestReadGambitNeutralBoundaryConditions tests reading boundary conditions
func TestReadGambitNeutralBoundaryConditions(t *testing.T) {
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test boundary conditions
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         1         2         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         6         4         1         2         3         4
ENDOFSECTION
       BOUNDARY CONDITIONS 2.0.0
inlet           1        10         0         0         0         0         0         0
wall            1        20         0         0         0         0         0         0
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGambitNeutral(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gambit neutral file: %v", err)
	}

	// Check boundary tags
	if len(mesh.BoundaryTags) != 2 {
		t.Errorf("Expected 2 boundary conditions, got %d", len(mesh.BoundaryTags))
	}

	// Check that we have the expected BC names
	foundInlet := false
	foundWall := false
	for _, name := range mesh.BoundaryTags {
		if name == "inlet" {
			foundInlet = true
		} else if name == "wall" {
			foundWall = true
		}
	}

	if !foundInlet {
		t.Error("Expected to find 'inlet' boundary condition")
	}
	if !foundWall {
		t.Error("Expected to find 'wall' boundary condition")
	}
}

// TestReadGambitNeutralNodeOrdering tests node ordering conversion
func TestReadGambitNeutralNodeOrdering(t *testing.T) {
	// Test that 1-based node IDs in file are converted to 0-based indices
	content := `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test node ordering
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         4         1         1         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   0.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         2   1.00000000000e+00   0.00000000000e+00   0.00000000000e+00
         3   0.00000000000e+00   1.00000000000e+00   0.00000000000e+00
         4   0.00000000000e+00   0.00000000000e+00   1.00000000000e+00
ENDOFSECTION
   ELEMENTS/CELLS 2.0.0
         1         6         4         1         2         3         4
ENDOFSECTION`

	tmpFile := createTempNeuFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGambitNeutral(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gambit neutral file: %v", err)
	}

	// Check element connectivity uses 0-based indices
	if len(mesh.EtoV) != 1 {
		t.Fatalf("Expected 1 element, got %d", len(mesh.EtoV))
	}

	expectedConnectivity := []int{0, 1, 2, 3} // 0-based
	elem := mesh.EtoV[0]
	if len(elem) != len(expectedConnectivity) {
		t.Fatalf("Expected %d nodes, got %d", len(expectedConnectivity), len(elem))
	}

	for i, expected := range expectedConnectivity {
		if elem[i] != expected {
			t.Errorf("Node %d: expected index %d, got %d", i, expected, elem[i])
		}
	}
}

// TestReadGambitNeutralErrorHandling tests various error conditions
func TestReadGambitNeutralErrorHandling(t *testing.T) {
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
			name: "MissingEndOfSection",
			content: `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         1         0         0         0         3         3`,
			errMsg: "", // Should handle gracefully
		},
		{
			name: "InvalidNodeData",
			content: `        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
         1         0         0         0         3         3
ENDOFSECTION
   NODAL COORDINATES 2.0.0
         1   not_a_number   0.0   0.0
ENDOFSECTION`,
			errMsg: "", // Should skip invalid nodes
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			var tmpFile string
			if tc.content != "" {
				tmpFile = createTempNeuFile(t, tc.content)
				defer os.Remove(tmpFile)
			} else {
				tmpFile = "nonexistent.neu"
			}

			_, err := ReadGambitNeutral(tmpFile)
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

// TestReadGambitNeutralRealFile tests reading the actual cube-partitioned.neu file
func TestReadGambitNeutralRealFile(t *testing.T) {
	// This test requires the actual file to be present
	filename := "cube-partitioned.neu"
	if _, err := os.Stat(filename); os.IsNotExist(err) {
		t.Skip("Skipping test - cube-partitioned.neu not found")
	}

	mesh, err := ReadGambitNeutral(filename)
	if err != nil {
		t.Fatalf("Failed to read %s: %v", filename, err)
	}

	// Verify expected values from the original test
	if mesh.NumElements != 565 {
		t.Errorf("Expected 565 elements, got %d", mesh.NumElements)
	}
	if mesh.NumVertices != 175 {
		t.Errorf("Expected 175 vertices, got %d", mesh.NumVertices)
	}

	// Additional validation
	// Check that all elements are 3D
	for i, elemType := range mesh.ElementTypes {
		if elemType.GetDimension() != 3 {
			t.Errorf("Element %d: expected 3D element, got %dD (%v)",
				i, elemType.GetDimension(), elemType)
		}
	}

	// Check node coordinates are reasonable
	for i, vertex := range mesh.Vertices {
		if len(vertex) != 3 {
			t.Errorf("Vertex %d: expected 3 coordinates, got %d", i, len(vertex))
		}
		// Check for NaN or Inf
		for j, coord := range vertex {
			if fmt.Sprintf("%f", coord) == "NaN" || fmt.Sprintf("%f", coord) == "+Inf" || fmt.Sprintf("%f", coord) == "-Inf" {
				t.Errorf("Vertex %d coord %d: invalid value %f", i, j, coord)
			}
		}
	}
}

// Builder methods

// Helper function to filter a CompleteMesh to only include used nodes
func filterToUsedNodes(msh *mesh.CompleteMesh) mesh.CompleteMesh {
	// Collect all nodes used by elements
	usedNodes := make(map[string]bool)
	for _, elemSet := range msh.Elements {
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
	for name := range msh.Nodes.NodeMap {
		if usedNodes[name] {
			orderedNames = append(orderedNames, name)
		}
	}

	// Sort for consistent ordering
	sort.Strings(orderedNames)

	// Build new node arrays
	for _, name := range orderedNames {
		origIdx := msh.Nodes.NodeMap[name]
		nodes = append(nodes, msh.Nodes.Nodes[origIdx])
		nodeMap[name] = idx
		nodeIDMap[name] = idx + 1
		idx++
	}

	filteredNodes := mesh.NodeSet{
		Nodes:     nodes,
		NodeMap:   nodeMap,
		NodeIDMap: nodeIDMap,
	}

	return mesh.CompleteMesh{
		Nodes:       filteredNodes,
		Elements:    msh.Elements,
		Dimension:   msh.Dimension,
		BoundingBox: msh.BoundingBox,
	}
}

// Builder methods

func (b *GambitTestBuilder) BuildSingleTetTest() string {
	tm := b.tm
	mesh := mesh.CompleteMesh{
		Nodes:       tm.TetraNodes,
		Elements:    []mesh.ElementSet{tm.SingleTet},
		Dimension:   3,
		BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 1}},
	}
	filtered := filterToUsedNodes(&mesh)
	return b.BuildFromCompleteMesh(&filtered)
}

func (b *GambitTestBuilder) BuildMixedElementTest() string {
	tm := b.tm
	filtered := filterToUsedNodes(&tm.MixedMesh)
	return b.BuildFromCompleteMesh(&filtered)
}

func (b *GambitTestBuilder) BuildFromCompleteMesh(mesh *mesh.CompleteMesh) string {
	var sections []string

	// Count total elements
	totalElements := 0
	for _, elemSet := range mesh.Elements {
		totalElements += len(elemSet.Elements)
	}

	// Header section
	header := fmt.Sprintf(`        CONTROL INFO 2.0.0
** GAMBIT NEUTRAL FILE
Test mesh from test helpers
PROGRAM:                  Test     VERSION:  1.0
Mon Jan  1 00:00:00 2025
     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL
%10d%10d%10d%10d%10d%10d
ENDOFSECTION`,
		len(mesh.Nodes.Nodes), totalElements, 1, 0, 3, 3)
	sections = append(sections, header)

	// Nodal coordinates section
	sections = append(sections, "   NODAL COORDINATES 2.0.0")
	for i, node := range mesh.Nodes.Nodes {
		sections = append(sections, fmt.Sprintf("%10d   %e   %e   %e",
			i+1, node[0], node[1], node[2]))
	}
	sections = append(sections, "ENDOFSECTION")

	// Elements section
	sections = append(sections, "   ELEMENTS/CELLS 2.0.0")
	elemID := 1
	for _, elemSet := range mesh.Elements {
		gambitType := elementTypeToGambit(elemSet.Type)
		numNodes := elemSet.Type.GetNumNodes()

		for _, elem := range elemSet.Elements {
			line := fmt.Sprintf("%10d%10d%10d", elemID, gambitType, numNodes)
			// Convert logical node names to indices
			for _, nodeName := range elem {
				nodeIdx := mesh.Nodes.NodeMap[nodeName] + 1 // 1-based
				line += fmt.Sprintf("%10d", nodeIdx)
			}
			sections = append(sections, line)
			elemID++
		}
	}
	sections = append(sections, "ENDOFSECTION")

	// Element group section
	sections = append(sections, "       ELEMENT GROUP 2.0.0")
	sections = append(sections, fmt.Sprintf("GROUP:%11d ELEMENTS:%11d MATERIAL:%11d NFLAGS:%11d",
		1, totalElements, 2, 0))
	sections = append(sections, "fluid")

	// List all element IDs
	line := ""
	for i := 1; i <= totalElements; i++ {
		if len(line) > 0 && len(line)+11 > 80 {
			sections = append(sections, line)
			line = ""
		}
		line += fmt.Sprintf("%11d", i)
	}
	if line != "" {
		sections = append(sections, line)
	}
	sections = append(sections, "ENDOFSECTION")

	return strings.Join(sections, "\n")
}

// elementTypeToGambit converts internal element types to Gambit type codes
func elementTypeToGambit(et mesh.ElementType) int {
	switch et {
	case mesh.Hex:
		return 4 // Brick
	case mesh.Prism:
		return 5 // Wedge
	case mesh.Tet:
		return 6 // Tetrahedron
	case mesh.Pyramid:
		return 7 // Pyramid
	default:
		return 0
	}
}
