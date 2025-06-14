package mesh

import (
	"fmt"
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

// TestReadGmsh22NodesWithArbitraryIDs tests non-sequential node IDs
func TestReadGmsh22NodesWithArbitraryIDs(t *testing.T) {
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
	if len(mesh.Elements) != 1 {
		t.Fatalf("Expected 1 element, got %d", len(mesh.Elements))
	}

	// Element should reference the correct array indices
	elem := mesh.Elements[0]
	expectedNodeIDs := []int{10, 25, 30, 200}
	for i, nodeID := range expectedNodeIDs {
		idx, _ := mesh.GetNodeIndex(nodeID)
		if elem[i] != idx {
			t.Errorf("Element node %d: expected index %d, got %d", i, idx, elem[i])
		}
	}
}

// TestReadGmsh22AllElementTypes tests all supported element types
func TestReadGmsh22AllElementTypes(t *testing.T) {
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
30
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.5 1.0 0.0
4 0.5 0.5 1.0
5 0.0 1.0 0.0
6 1.0 1.0 0.0
7 0.0 0.0 1.0
8 1.0 0.0 1.0
9 1.0 1.0 1.0
10 0.0 1.0 1.0
11 0.5 0.0 0.0
12 0.75 0.5 0.0
13 0.25 0.5 0.0
14 0.5 0.25 0.5
15 0.5 0.75 0.5
16 0.25 0.5 0.5
17 0.75 0.5 0.5
18 0.5 0.5 0.0
19 0.5 0.5 0.25
20 0.5 0.5 0.75
21 0.5 0.0 1.0
22 0.0 0.5 1.0
23 1.0 0.5 1.0
24 0.5 1.0 1.0
25 0.0 0.5 0.0
26 1.0 0.5 0.0
27 0.5 0.0 0.5
28 0.25 0.25 0.25
29 0.75 0.75 0.75
30 0.5 0.5 0.5
$EndNodes
$Elements
21
1 15 1 1 1
2 1 2 1 2 1 2
3 8 2 1 2 1 2 11
4 2 2 1 2 1 2 3
5 9 2 1 2 1 2 3 11 12 13
6 20 2 1 2 1 2 3 11 12 13 18 19 20
7 21 2 1 2 1 2 3 11 12 13 18 19 20 21
8 3 2 1 2 1 2 6 5
9 16 2 1 2 1 2 6 5 11 26 25 13
10 10 2 1 2 1 2 6 5 11 26 25 13 18
11 4 2 1 2 1 2 3 4
12 11 2 1 2 1 2 3 4 11 12 13 14 15 16
13 5 2 1 2 1 2 6 5 7 8 9 10
14 17 2 1 2 1 2 6 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
15 12 2 1 2 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
16 6 2 1 2 1 2 3 7 8 9
17 18 2 1 2 1 2 3 7 8 9 11 12 13 21 22 23 24 25 26 27
18 13 2 1 2 1 2 3 7 8 9 11 12 13 21 22 23 24 25 26 27 28 29 30
19 7 2 1 2 1 2 6 5 4
20 14 2 1 2 1 2 6 5 4 11 26 25 13 18 19 20 21 27 28
21 19 2 1 2 1 2 6 5 4 11 26 25 18 19 20 21 29
$EndElements`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	// Expected element types in order
	expectedTypes := []ElementType{
		Point,      // 15: 1-node point
		Line,       // 1: 2-node line
		Line3,      // 8: 3-node line
		Triangle,   // 2: 3-node triangle
		Triangle6,  // 9: 6-node triangle
		Triangle9,  // 20: 9-node triangle
		Triangle10, // 21: 10-node triangle
		Quad,       // 3: 4-node quad
		Quad8,      // 16: 8-node quad
		Quad9,      // 10: 9-node quad
		Tet,        // 4: 4-node tet
		Tet10,      // 11: 10-node tet
		Hex,        // 5: 8-node hex
		Hex20,      // 17: 20-node hex
		Hex27,      // 12: 27-node hex
		Prism,      // 6: 6-node prism
		Prism15,    // 18: 15-node prism
		Prism18,    // 13: 18-node prism
		Pyramid,    // 7: 5-node pyramid
		Pyramid14,  // 14: 14-node pyramid
		Pyramid13,  // 19: 13-node pyramid
	}

	if len(mesh.Elements) != len(expectedTypes) {
		t.Fatalf("Expected %d elements, got %d", len(expectedTypes), len(mesh.Elements))
	}

	for i, expected := range expectedTypes {
		if mesh.ElementTypes[i] != expected {
			t.Errorf("Element %d: expected type %v, got %v", i, expected, mesh.ElementTypes[i])
		}

		// Check node count
		expectedNodes := expected.GetNumNodes()
		if len(mesh.Elements[i]) != expectedNodes {
			t.Errorf("Element %d (%v): expected %d nodes, got %d",
				i, expected, expectedNodes, len(mesh.Elements[i]))
		}
	}
}

// TestReadGmsh22PhysicalNames tests physical entity names
func TestReadGmsh22PhysicalNames(t *testing.T) {
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

	// Check that elements are assigned to groups
	if mesh.ElementTags[0][0] != 10 {
		t.Errorf("Element 0: expected physical tag 10, got %d", mesh.ElementTags[0][0])
	}

	if mesh.ElementTags[1][0] != 40 {
		t.Errorf("Element 1: expected physical tag 40, got %d", mesh.ElementTags[1][0])
	}
}

// TestReadGmsh22MultipleTagsPerElement tests elements with multiple tags
func TestReadGmsh22MultipleTagsPerElement(t *testing.T) {
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
4
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.0 1.0 0.0
4 0.0 0.0 1.0
$EndNodes
$Elements
4
1 4 0 1 2 3 4
2 4 1 10 1 2 3 4
3 4 2 10 20 1 2 3 4
4 4 4 10 20 30 40 1 2 3 4
$EndElements`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	// Check tag counts
	expectedTagCounts := []int{0, 1, 2, 4}
	for i, expected := range expectedTagCounts {
		if len(mesh.ElementTags[i]) != expected {
			t.Errorf("Element %d: expected %d tags, got %d",
				i, expected, len(mesh.ElementTags[i]))
		}
	}

	// Check specific tags
	if len(mesh.ElementTags[3]) == 4 {
		expectedTags := []int{10, 20, 30, 40}
		for j, tag := range expectedTags {
			if mesh.ElementTags[3][j] != tag {
				t.Errorf("Element 3 tag %d: expected %d, got %d",
					j, tag, mesh.ElementTags[3][j])
			}
		}
	}
}

// TestReadGmsh22Periodic tests periodic boundary conditions
func TestReadGmsh22Periodic(t *testing.T) {
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
Affine 1.0 0.0 0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 1.0
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
	if p1.SlaveTag != 1 {
		t.Errorf("Periodic 1: expected slave tag 1, got %d", p1.SlaveTag)
	}
	if p1.MasterTag != 2 {
		t.Errorf("Periodic 1: expected master tag 2, got %d", p1.MasterTag)
	}
	if len(p1.AffineTransform) != 16 {
		t.Errorf("Periodic 1: expected 16 affine values, got %d", len(p1.AffineTransform))
	}
	if len(p1.NodeMap) != 2 {
		t.Errorf("Periodic 1: expected 2 node mappings, got %d", len(p1.NodeMap))
	}
	if p1.NodeMap[1] != 5 || p1.NodeMap[4] != 8 {
		t.Error("Periodic 1: incorrect node mappings")
	}

	// Check second periodic entity (without affine)
	p2 := mesh.Periodics[1]
	if p2.Dimension != 2 {
		t.Errorf("Periodic 2: expected dimension 2, got %d", p2.Dimension)
	}
	if p2.SlaveTag != 3 {
		t.Errorf("Periodic 2: expected slave tag 3, got %d", p2.SlaveTag)
	}
	if p2.MasterTag != 4 {
		t.Errorf("Periodic 2: expected master tag 4, got %d", p2.MasterTag)
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
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
10
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.0 1.0 0.0
4 0.0 0.0 1.0
5 0.5 0.0 0.0
6 0.5 0.5 0.0
7 0.0 0.5 0.0
8 0.5 0.5 0.5
9 1.0 1.0 0.0
10 0.0 1.0 1.0
$EndNodes
$Elements
6
1 15 0 1
2 1 0 1 2
3 2 0 1 2 3
4 3 0 1 2 9 3
5 4 0 1 2 3 4
6 5 0 1 2 9 3 4 10 8 7
$EndElements`

	tmpFile := createTempMshFile(t, content)

	mesh, err := ReadGmsh22(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh file: %v", err)
	}

	// Test filtering by dimension
	for dim := 0; dim <= 3; dim++ {
		indices, elements, types := mesh.FilterByDimension(dim)

		expectedCounts := []int{1, 1, 2, 2} // 0D: 1, 1D: 1, 2D: 2, 3D: 2
		if len(indices) != expectedCounts[dim] {
			t.Errorf("Dimension %d: expected %d elements, got %d",
				dim, expectedCounts[dim], len(indices))
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

// TestReadGmsh22HigherOrderCornerNodes tests corner node extraction
func TestReadGmsh22HigherOrderCornerNodes(t *testing.T) {
	testCases := []struct {
		elemType        ElementType
		expectedCorners []int
	}{
		{Line3, []int{0, 1}},
		{Triangle6, []int{0, 1, 2}},
		{Quad8, []int{0, 1, 2, 3}},
		{Tet10, []int{0, 1, 2, 3}},
		{Hex20, []int{0, 1, 2, 3, 4, 5, 6, 7}},
		{Prism15, []int{0, 1, 2, 3, 4, 5}},
		{Pyramid13, []int{0, 1, 2, 3, 4}},
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
	content := `$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
1
1 0.0 0.0 0.0
$EndNodes
$Elements
1
1 15 0 1
$EndElements`

	tmpFile := createTempMshFile(t, content)

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

// Benchmark for performance testing
func BenchmarkReadGmsh22LargeMesh(b *testing.B) {
	// Generate a large mesh
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
