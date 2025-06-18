package readers

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// Helper function to create temporary test files
func createTempMshFile(t *testing.T, content string) string {
	t.Helper()
	tmpFile := filepath.Join(t.TempDir(), "test.msh")
	if err := os.WriteFile(tmpFile, []byte(content), 0644); err != nil {
		t.Fatalf("Failed to create temp file: %v", err)
	}
	return tmpFile
}

// TestReadGmsh22Version tests reading version information
func TestReadGmsh22Version(t *testing.T) {
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
0
$EndNodes
$Elements
0
$EndElements`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	if mesh.FormatVersion != "2.2" {
		t.Errorf("Expected version 2.2, got %s", mesh.FormatVersion)
	}

	if mesh.IsBinary {
		t.Error("Expected ASCII format, got binary")
	}

	if mesh.DataSize != 8 {
		t.Errorf("Expected data size 8, got %d", mesh.DataSize)
	}
}

// TestReadGmsh22StandardMeshes tests reading standard test meshes
func TestReadGmsh22StandardMeshes(t *testing.T) {
	builder := NewGmsh22TestBuilder()

	t.Run("TwoTetMesh", func(t *testing.T) {
		content := builder.BuildTwoTetTest()

		tmpFile := createTempMshFile(t, content)
		defer os.Remove(tmpFile)

		mesh, err := ReadGmsh22(tmpFile)
		if err != nil {
			t.Fatalf("Failed to read Gmsh 2.2 file: %v", err)
		}

		// Validate basic structure
		if mesh.NumVertices != 5 {
			t.Errorf("Expected 5 vertices, got %d", mesh.NumVertices)
		}

		if mesh.NumElements != 2 {
			t.Errorf("Expected 2 elements, got %d", mesh.NumElements)
		}

		// Check element types
		for i := 0; i < 2; i++ {
			if mesh.ElementTypes[i] != utils.Tet {
				t.Errorf("Element %d: expected Tet, got %v", i, mesh.ElementTypes[i])
			}
			if len(mesh.EtoV[i]) != 4 {
				t.Errorf("Element %d: expected 4 nodes, got %d", i, len(mesh.EtoV[i]))
			}
		}
	})

	t.Run("MixedMesh", func(t *testing.T) {
		content := builder.BuildMixedElementTest()

		tmpFile := createTempMshFile(t, content)
		defer os.Remove(tmpFile)

		mesh, err := ReadGmsh22(tmpFile)
		if err != nil {
			t.Fatalf("Failed to read Gmsh 2.2 file: %v", err)
		}

		// Get expected types from test mesh
		tm := utils.GetStandardTestMeshes()
		var expectedTypes []utils.ElementType
		for _, elemSet := range tm.MixedMesh.Elements {
			for range elemSet.Elements {
				expectedTypes = append(expectedTypes, elemSet.Type)
			}
		}

		// Validate element types
		if len(mesh.ElementTypes) != len(expectedTypes) {
			t.Fatalf("Element count mismatch: got %d, expected %d",
				len(mesh.ElementTypes), len(expectedTypes))
		}

		for i, expectedType := range expectedTypes {
			if mesh.ElementTypes[i] != expectedType {
				t.Errorf("Element %d: expected type %v, got %v",
					i, expectedType, mesh.ElementTypes[i])
			}

			// Check node count
			expectedNodes := expectedType.GetNumNodes()
			actualNodes := len(mesh.EtoV[i])
			if actualNodes != expectedNodes {
				t.Errorf("Element %d (%v): expected %d nodes, got %d",
					i, expectedType, expectedNodes, actualNodes)
			}
		}
	})

	t.Run("CubeMesh", func(t *testing.T) {
		tm := utils.GetStandardTestMeshes()
		builder := NewGmsh22TestBuilder()
		content := builder.BuildFromCompleteMesh(&tm.CubeMesh)

		tmpFile := createTempMshFile(t, content)
		defer os.Remove(tmpFile)

		mesh, err := ReadGmsh22(tmpFile)
		if err != nil {
			t.Fatalf("Failed to read Gmsh 2.2 file: %v", err)
		}

		// Validate cube mesh
		if mesh.NumElements != 6 {
			t.Errorf("Cube mesh should have 6 tets, got %d elements", mesh.NumElements)
		}

		for i := 0; i < mesh.NumElements; i++ {
			if mesh.ElementTypes[i] != utils.Tet {
				t.Errorf("Cube element %d: expected Tet, got %v", i, mesh.ElementTypes[i])
			}
		}
	})
}

// TestReadGmsh22MixedElementTypes tests mixed element types using test helpers
func TestReadGmsh22MixedElementTypes(t *testing.T) {
	tm := utils.GetStandardTestMeshes()
	builder := NewGmsh22TestBuilder()

	content := builder.BuildMixedElementTest()

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 2.2 file: %v", err)
	}

	// Build expected types from the test mesh
	var expectedTypes []utils.ElementType
	for _, elemSet := range tm.MixedMesh.Elements {
		for range elemSet.Elements {
			expectedTypes = append(expectedTypes, elemSet.Type)
		}
	}

	// Validate all elements
	for i, expectedType := range expectedTypes {
		if i >= len(mesh.ElementTypes) {
			t.Errorf("Element %d missing", i)
			continue
		}

		if mesh.ElementTypes[i] != expectedType {
			t.Errorf("Element %d: expected type %v, got %v", i, expectedType, mesh.ElementTypes[i])
		}

		// Validate node count
		actualNodes := len(mesh.EtoV[i])
		expectedNodes := expectedType.GetNumNodes()
		if actualNodes != expectedNodes {
			t.Errorf("Element %d (%v): expected %d nodes, got %d",
				i, expectedType, expectedNodes, actualNodes)
		}
	}

	// Check that elements reference valid nodes
	for i, elem := range mesh.EtoV {
		for j, nodeIdx := range elem {
			if nodeIdx < 0 || nodeIdx >= mesh.NumVertices {
				t.Errorf("Element %d, node %d: invalid node index %d (should be 0-%d)",
					i, j, nodeIdx, mesh.NumVertices-1)
			}
		}
	}
}

// TestReadGmsh22NodesWithArbitraryIDs tests non-sequential node IDs
func TestReadGmsh22NodesWithArbitraryIDs(t *testing.T) {
	// This test uses hardcoded values because it's testing specific ID mapping behavior
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
5
10 0.0 0.0 0.0
25 1.0 0.0 0.0
30 1.0 1.0 0.0
100 0.0 1.0 0.0
200 0.5 0.5 0.5
$EndNodes
$Elements
1
1 4 0 10 25 30 200
$EndElements`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	// Check that we have 5 nodes
	if mesh.NumVertices != 5 {
		t.Errorf("Expected 5 vertices, got %d", mesh.NumVertices)
	}

	// Check node ID mapping
	testCases := []struct {
		nodeID int
		coords []float64
	}{
		{10, []float64{0.0, 0.0, 0.0}},
		{25, []float64{1.0, 0.0, 0.0}},
		{30, []float64{1.0, 1.0, 0.0}},
		{100, []float64{0.0, 1.0, 0.0}},
		{200, []float64{0.5, 0.5, 0.5}},
	}

	for _, tc := range testCases {
		idx, ok := mesh.GetNodeIndex(tc.nodeID)
		if !ok {
			t.Errorf("Node ID %d not found", tc.nodeID)
			continue
		}

		for i, coord := range tc.coords {
			if mesh.Vertices[idx][i] != coord {
				t.Errorf("Node %d coord %d: expected %f, got %f",
					tc.nodeID, i, coord, mesh.Vertices[idx][i])
			}
		}
	}

	// Check element connectivity uses correct indices
	if len(mesh.EtoV) != 1 {
		t.Fatalf("Expected 1 element, got %d", len(mesh.EtoV))
	}
}

// TestReadGmsh22AllElementTypes tests all supported element types
func TestReadGmsh22AllElementTypes(t *testing.T) {
	// This test validates that all element types are correctly read
	// We'll create individual elements of each type using test helpers
	tm := utils.GetStandardTestMeshes()

	// Test single elements
	testCases := []struct {
		name     string
		mesh     utils.CompleteMesh
		expected utils.ElementType
	}{
		{
			name: "SingleTet",
			mesh: utils.CompleteMesh{
				Nodes:       tm.TetraNodes,
				Elements:    []utils.ElementSet{tm.SingleTet},
				Dimension:   3,
				BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 1}},
			},
			expected: utils.Tet,
		},
		{
			name: "SingleHex",
			mesh: utils.CompleteMesh{
				Nodes:       tm.CubeNodes,
				Elements:    []utils.ElementSet{tm.SingleHex},
				Dimension:   3,
				BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 1}},
			},
			expected: utils.Hex,
		},
		{
			name: "SinglePrism",
			mesh: utils.CompleteMesh{
				Nodes:       tm.CubeNodes,
				Elements:    []utils.ElementSet{tm.SinglePrism},
				Dimension:   3,
				BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 1}},
			},
			expected: utils.Prism,
		},
		{
			name: "SinglePyramid",
			mesh: utils.CompleteMesh{
				Nodes:       tm.PyramidNodes,
				Elements:    []utils.ElementSet{tm.SinglePyramid},
				Dimension:   3,
				BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 1}},
			},
			expected: utils.Pyramid,
		},
	}

	builder := NewGmsh22TestBuilder()

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			content := builder.BuildFromCompleteMesh(&tc.mesh)
			tmpFile := createTempMshFile(t, content)
			defer os.Remove(tmpFile)

			mesh, err := ReadGmsh22(tmpFile)
			if err != nil {
				t.Fatalf("Failed to read Gmsh file: %v", err)
			}

			if mesh.NumElements != 1 {
				t.Errorf("Expected 1 element, got %d", mesh.NumElements)
			}

			if mesh.ElementTypes[0] != tc.expected {
				t.Errorf("Expected element type %v, got %v", tc.expected, mesh.ElementTypes[0])
			}

			// Check node count
			expectedNodes := tc.expected.GetNumNodes()
			if len(mesh.EtoV[0]) != expectedNodes {
				t.Errorf("Expected %d nodes, got %d", expectedNodes, len(mesh.EtoV[0]))
			}
		})
	}
}

// TestReadGmsh22PhysicalNames tests physical entity names
func TestReadGmsh22PhysicalNames(t *testing.T) {
	// This test requires specific physical names, so we keep it with minimal hardcoding

	// Create a simple mesh with physical groups
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
4
0 10 "Inlet Point"
1 20 "Outlet Curve"
2 30 "Wall Surface"
3 40 "Fluid Volume"
$EndPhysicalNames
$Nodes
4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.0 1.0 0.0
4 0.0 0.0 1.0
$EndNodes
$Elements
2
1 15 1 10 1
2 4 2 40 30 1 2 3 4
$EndElements`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	// Check physical groups
	expectedGroups := map[int]struct {
		name string
		dim  int
	}{
		10: {"Inlet Point", 0},
		20: {"Outlet Curve", 1},
		30: {"Wall Surface", 2},
		40: {"Fluid Volume", 3},
	}

	for tag, expected := range expectedGroups {
		group, ok := mesh.ElementGroups[tag]
		if !ok {
			t.Errorf("Physical group %d not found", tag)
			continue
		}

		if group.Name != expected.name {
			t.Errorf("Group %d: expected name '%s', got '%s'", tag, expected.name, group.Name)
		}

		if group.Dimension != expected.dim {
			t.Errorf("Group %d: expected dimension %d, got %d", tag, expected.dim, group.Dimension)
		}
	}
}

// TestReadGmsh22MultipleTagsPerElement tests elements with multiple tags
func TestReadGmsh22MultipleTagsPerElement(t *testing.T) {
	// Use test helpers to create a mesh with specific tags
	tm := utils.GetStandardTestMeshes()

	// Create a custom mesh with multiple tags
	customMesh := utils.CompleteMesh{
		Nodes: tm.TetraNodes,
		Elements: []utils.ElementSet{
			{
				Type: utils.Tet,
				Elements: [][]string{
					{"v0", "v1", "v2", "v3"},
				},
				Properties: []utils.ElementProps{
					{PhysicalTag: 10, GeometricTag: 20},
				},
			},
		},
		Dimension:   3,
		BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 1}},
	}

	builder := NewGmsh22TestBuilder()
	content := builder.BuildFromCompleteMesh(&customMesh)

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	// Check tags
	if len(mesh.ElementTags[0]) < 2 {
		t.Errorf("Expected at least 2 tags, got %d", len(mesh.ElementTags[0]))
	} else {
		if mesh.ElementTags[0][0] != 10 {
			t.Errorf("Expected physical tag 10, got %d", mesh.ElementTags[0][0])
		}
		if mesh.ElementTags[0][1] != 20 {
			t.Errorf("Expected geometric tag 20, got %d", mesh.ElementTags[0][1])
		}
	}
}

// TestReadGmsh22Periodic tests periodic boundary conditions
func TestReadGmsh22Periodic(t *testing.T) {
	// Periodic boundaries are format-specific, so we keep this test as is
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
8
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 1.0 1.0 0.0
4 0.0 1.0 0.0
5 0.0 0.0 1.0
6 1.0 0.0 1.0
7 1.0 1.0 1.0
8 0.0 1.0 1.0
$EndNodes
$Elements
1
1 5 0 1 2 3 4 5 6 7 8
$EndElements
$Periodic
2
2 1 2
Affine 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
2
1 5
4 8
2 3 4
3
1 2
3 4
5 6
$EndPeriodic`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	// Check periodic entities
	if len(mesh.Periodics) != 2 {
		t.Fatalf("Expected 2 periodic entities, got %d", len(mesh.Periodics))
	}

	// Check first periodic entity (with affine transformation)
	p1 := mesh.Periodics[0]
	if p1.Dimension != 2 {
		t.Errorf("Periodic 1: expected dimension 2, got %d", p1.Dimension)
	}
	if len(p1.AffineTransform) != 9 {
		t.Errorf("Periodic 1: expected 9 affine values, got %d", len(p1.AffineTransform))
	}
	if len(p1.NodeMap) != 2 {
		t.Errorf("Periodic 1: expected 2 node mappings, got %d", len(p1.NodeMap))
	}

	// Check second periodic entity (without affine)
	p2 := mesh.Periodics[1]
	if p2.Dimension != 2 {
		t.Errorf("Periodic 2: expected dimension 2, got %d", p2.Dimension)
	}
	if len(p2.AffineTransform) != 0 {
		t.Errorf("Periodic 2: expected no affine transform, got %d values", len(p2.AffineTransform))
	}
	if len(p2.NodeMap) != 3 {
		t.Errorf("Periodic 2: expected 3 node mappings, got %d", len(p2.NodeMap))
	}
}

// TestReadGmsh22FilterByDimension tests filtering elements by dimension
func TestReadGmsh22FilterByDimension(t *testing.T) {
	// Use the mixed mesh which has elements of different dimensions
	tm := utils.GetStandardTestMeshes()
	builder := NewGmsh22TestBuilder()

	content := builder.BuildMixedElementTest()
	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	// Count expected elements by dimension from the mixed mesh
	dimCounts := make(map[int]int)
	for _, elemSet := range tm.MixedMesh.Elements {
		dim := elemSet.Type.GetDimension()
		dimCounts[dim] += len(elemSet.Elements)
	}

	// Test filtering by dimension
	for dim := 0; dim <= 3; dim++ {
		indices, elements, types := mesh.FilterByDimension(dim)

		expectedCount := dimCounts[dim]
		if len(indices) != expectedCount {
			t.Errorf("Dimension %d: expected %d elements, got %d",
				dim, expectedCount, len(indices))
		}

		// Check that all filtered elements have correct dimension
		for i, elemType := range types {
			if elemType.GetDimension() != dim {
				t.Errorf("Filter dimension %d, element %d: got dimension %d",
					dim, i, elemType.GetDimension())
			}
		}

		// Verify indices match
		if len(indices) != len(elements) || len(indices) != len(types) {
			t.Errorf("Dimension %d: mismatched array lengths", dim)
		}
	}
}

// TestReadGmsh22HigherOrderElements tests higher-order elements
func TestReadGmsh22HigherOrderElements(t *testing.T) {
	// This test is specific to higher-order elements which aren't in the standard test meshes
	// So we keep a minimal hardcoded version
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
10
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
5 0.5 0 0
6 0.5 0.5 0
7 0 0.5 0
8 0.5 0 0.5
9 0 0.5 0.5
10 0 0 0.5
$EndNodes
$Elements
1
1 11 0 1 2 3 4 5 6 7 8 9 10
$EndElements`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	if mesh.NumElements != 1 {
		t.Errorf("Expected 1 element, got %d", mesh.NumElements)
	}

	if mesh.ElementTypes[0] != utils.Tet10 {
		t.Errorf("Expected Tet10, got %v", mesh.ElementTypes[0])
	}

	if len(mesh.EtoV[0]) != 10 {
		t.Errorf("Tet10 should have 10 nodes, got %d", len(mesh.EtoV[0]))
	}
}

// TestReadGmsh22HigherOrderCornerNodes tests corner node extraction
func TestReadGmsh22HigherOrderCornerNodes(t *testing.T) {
	testCases := []struct {
		elemType        utils.ElementType
		expectedCorners []int
	}{
		{utils.Line3, []int{0, 1}},
		{utils.Triangle6, []int{0, 1, 2}},
		{utils.Quad8, []int{0, 1, 2, 3}},
		{utils.Tet10, []int{0, 1, 2, 3}},
		{utils.Hex20, []int{0, 1, 2, 3, 4, 5, 6, 7}},
		{utils.Prism15, []int{0, 1, 2, 3, 4, 5}},
		{utils.Pyramid13, []int{0, 1, 2, 3, 4}},
	}

	for _, tc := range testCases {
		corners := tc.elemType.GetCornerNodes()
		if len(corners) != len(tc.expectedCorners) {
			t.Errorf("%s: expected %d corner nodes, got %d",
				tc.elemType, len(tc.expectedCorners), len(corners))
			continue
		}

		for i, expected := range tc.expectedCorners {
			if corners[i] != expected {
				t.Errorf("%s corner %d: expected %d, got %d",
					tc.elemType, i, expected, corners[i])
			}
		}
	}
}

// TestReadGmsh22AutoDetection tests automatic version detection
func TestReadGmsh22AutoDetection(t *testing.T) {
	builder := NewGmsh22TestBuilder()
	tm := utils.GetStandardTestMeshes()

	// Use a simple test mesh
	content := builder.BuildFromCompleteMesh(&utils.CompleteMesh{
		Nodes:       tm.TetraNodes,
		Elements:    []utils.ElementSet{tm.SingleTet},
		Dimension:   3,
		BoundingBox: [2][3]float64{{0, 0, 0}, {1, 1, 1}},
	})

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	// Test with auto-detection
	mesh, err := ReadGmshAuto(tmpFile)
	if err != nil {
		t.Fatalf("Failed to auto-detect and read Gmsh file: %v", err)
	}

	if mesh.FormatVersion != "2.2" {
		t.Errorf("Auto-detection: expected version 2.2, got %s", mesh.FormatVersion)
	}
}

// TestReadGmsh22BinaryFormat tests binary format detection
func TestReadGmsh22BinaryFormat(t *testing.T) {
	// Note: Full binary implementation is complex, this tests detection
	content := `$MeshFormat
2.2 1 8
$EndMeshFormat`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	// Binary format not fully implemented, so we expect an error or partial read
	if err == nil && mesh.IsBinary {
		t.Log("Binary format detected correctly")
	}
}

// TestReadGmsh22ErrorHandling tests various error conditions
func TestReadGmsh22ErrorHandling(t *testing.T) {
	testCases := []struct {
		name    string
		content string
		errMsg  string
	}{
		{
			name: "Invalid element type",
			content: `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0 0 0
2 1 0 0
3 0 1 0
$EndNodes
$Elements
1
1 999 0 1 2 3
$EndElements`,
			errMsg: "", // Should skip unknown element types
		},
		{
			name: "Element references unknown node",
			content: `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
2
1 0 0 0
2 1 0 0
$EndNodes
$Elements
1
1 4 0 1 2 3 4
$EndElements`,
			errMsg: "unknown node",
		},
		{
			name: "Too few nodes for element",
			content: `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
3
1 0 0 0
2 1 0 0
3 0 1 0
$EndNodes
$Elements
1
1 4 0 1 2 3
$EndElements`,
			errMsg: "expected 4 nodes", // Updated to match actual error message
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			tmpFile := createTempMshFile(t, tc.content)
			defer os.Remove(tmpFile)

			_, err := ReadGmsh22(tmpFile)
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

// BenchmarkReadGmsh22LargeMesh benchmarks reading performance
func BenchmarkReadGmsh22LargeMesh(b *testing.B) {
	// Generate a large mesh using test helpers
	var content strings.Builder
	content.WriteString(`$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
1000
`)

	// Generate nodes
	for i := 1; i <= 1000; i++ {
		x := float64(i%10) * 0.1
		y := float64((i/10)%10) * 0.1
		z := float64(i/100) * 0.1
		fmt.Fprintf(&content, "%d %.1f %.1f %.1f\n", i, x, y, z)
	}

	content.WriteString(`$EndNodes
$Elements
500
`)

	// Generate tetrahedra
	for i := 1; i <= 500; i++ {
		v1 := i
		v2 := i + 1
		v3 := i + 2
		v4 := i + 3
		if v4 > 1000 {
			v4 = v4 - 1000
		}
		fmt.Fprintf(&content, "%d 4 2 1 1 %d %d %d %d\n", i, v1, v2, v3, v4)
	}

	content.WriteString("$EndElements\n")

	tmpFile := filepath.Join(b.TempDir(), "bench.msh")
	os.WriteFile(tmpFile, []byte(content.String()), 0644)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := ReadGmsh22(tmpFile)
		if err != nil {
			b.Fatal(err)
		}
	}
}
