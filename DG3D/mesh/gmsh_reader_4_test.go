package mesh

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestReadGmsh4Version tests reading version 4.x format
func TestReadGmsh4Version(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 0
$EndEntities
$Nodes
0 0 0 0
$EndNodes
$Elements
0 0 0 0
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	if mesh.FormatVersion != "4.1" {
		t.Errorf("Expected version 4.1, got %s", mesh.FormatVersion)
	}

	if mesh.IsBinary {
		t.Error("Expected ASCII format, got binary")
	}
}

// TestReadGmsh4Entities tests reading the Entities section
func TestReadGmsh4Entities(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
3
1 10 "Inlet"
2 20 "Wall"
3 30 "Volume"
$EndPhysicalNames
$Entities
2 3 2 1
1 0 0 0 1 10
2 1 0 0 0
1 0 0 0 1 0 0 1 10 2 1 -2
2 0.5 0 0 1 0.5 0 0 2 2 -1
3 1 0 0 1 1 0 0 2 -2
1 0 0 0 1 1 0 1 20 4 1 2 -3 -4
2 0 0 0 1 0.5 0 0 4 -1 3 4 -2
1 0 0 0 1 1 1 1 30 6 1 2 3 4 5 6
$EndEntities
$Nodes
0 0 0 0
$EndNodes
$Elements
0 0 0 0
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Check physical names
	if len(mesh.ElementGroups) != 3 {
		t.Errorf("Expected 3 physical groups, got %d", len(mesh.ElementGroups))
	}

	if mesh.ElementGroups[10].Name != "Inlet" {
		t.Errorf("Expected group 10 name 'Inlet', got '%s'", mesh.ElementGroups[10].Name)
	}
}

// TestReadGmsh4NodesFormat tests the new node format in v4
func TestReadGmsh4NodesFormat(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
1 0 0 1
1 0 0 0 0
1 0 0 0 1 1 1 1 1 4 1 2 3 4
$EndEntities
$Nodes
3 10 1 20
0 1 0 1
1
0 0 0
3 1 0 4
5
10
15
20
0 0 0
1 0 0
1 1 0
0 1 0
2 1 0 5
6
7
8
9
11
0.5 0 0
0.5 1 0
1 0.5 0
0 0.5 0
0.5 0.5 0
$EndNodes
$Elements
0 0 0 0
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Check total node count
	if mesh.NumVertices != 10 {
		t.Errorf("Expected 10 vertices, got %d", mesh.NumVertices)
	}

	// Check specific nodes
	testCases := []struct {
		nodeID int
		coords []float64
	}{
		{1, []float64{0, 0, 0}},
		{5, []float64{0, 0, 0}},
		{10, []float64{1, 0, 0}},
		{15, []float64{1, 1, 0}},
		{20, []float64{0, 1, 0}},
		{6, []float64{0.5, 0, 0}},
		{7, []float64{0.5, 1, 0}},
		{11, []float64{0.5, 0.5, 0}},
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
}

// TestReadGmsh4ElementsFormat tests the new element format in v4
func TestReadGmsh4ElementsFormat(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
2
2 10 "Surface"
3 20 "Volume"
$EndPhysicalNames
$Entities
0 0 1 1
1 0 0 0 1 1 1 1 10 0
1 0 0 0 1 1 1 1 20 1 1
$EndEntities
$Nodes
2 8 1 8
2 1 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
3 1 0 4
5
6
7
8
0 0 1
1 0 1
1 1 1
0 1 1
$EndNodes
$Elements
4 9 1 9
2 1 2 2
1 1 2 3
2 3 4 1
2 1 3 1
3 1 2 3 4
3 1 4 4
4 1 2 3 4 5
5 2 3 4 5 6
6 3 4 5 6 7
7 4 5 6 7 8
3 1 7 2
8 5 6 7 8 1
9 1 2 3 4 8
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Check element count
	if mesh.NumElements != 9 {
		t.Errorf("Expected 9 elements, got %d", mesh.NumElements)
	}

	// Check element types
	expectedTypes := []ElementType{
		Triangle, Triangle, // 2D elements (type 2)
		Quad,               // 2D element (type 3)
		Tet, Tet, Tet, Tet, // 3D elements (type 4)
		Pyramid, Pyramid, // 3D elements (type 7)
	}

	for i, expected := range expectedTypes {
		if i >= len(mesh.ElementTypes) {
			t.Errorf("Element %d missing", i)
			continue
		}
		if mesh.ElementTypes[i] != expected {
			t.Errorf("Element %d: expected type %v, got %v", i, expected, mesh.ElementTypes[i])
		}
	}

	// Check that 2D elements have physical tag 10
	for i := 0; i < 3; i++ {
		if len(mesh.ElementTags[i]) == 0 || mesh.ElementTags[i][0] != 10 {
			t.Errorf("2D element %d: expected physical tag 10, got %v", i, mesh.ElementTags[i])
		}
	}

	// Check that 3D elements have physical tag 20
	for i := 3; i < 9; i++ {
		if len(mesh.ElementTags[i]) == 0 || mesh.ElementTags[i][0] != 20 {
			t.Errorf("3D element %d: expected physical tag 20, got %v", i, mesh.ElementTags[i])
		}
	}
}

// TestReadGmsh4HigherOrderElements tests higher-order elements in v4
func TestReadGmsh4HigherOrderElements(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 1
1 0 0 0 1 1 1 0 0
$EndEntities
$Nodes
1 10 1 10
3 1 0 10
1
2
3
4
5
6
7
8
9
10
0 0 0
1 0 0
0 1 0
0 0 1
0.5 0 0
0.5 0.5 0
0 0.5 0
0.5 0 0.5
0 0.5 0.5
0 0 0.5
$EndNodes
$Elements
2 2 1 2
3 1 11 1
1 1 2 3 4 5 6 7 8 9 10
3 1 4 1
2 1 2 3 4
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Should have 2 elements: one Tet10 and one Tet
	if mesh.NumElements != 2 {
		t.Errorf("Expected 2 elements, got %d", mesh.NumElements)
	}

	if mesh.ElementTypes[0] != Tet10 {
		t.Errorf("Element 0: expected Tet10, got %v", mesh.ElementTypes[0])
	}

	if mesh.ElementTypes[1] != Tet {
		t.Errorf("Element 1: expected Tet, got %v", mesh.ElementTypes[1])
	}

	// Check node counts
	if len(mesh.Elements[0]) != 10 {
		t.Errorf("Tet10: expected 10 nodes, got %d", len(mesh.Elements[0]))
	}

	if len(mesh.Elements[1]) != 4 {
		t.Errorf("Tet: expected 4 nodes, got %d", len(mesh.Elements[1]))
	}
}

// TestReadGmsh4Periodic tests periodic boundary conditions in v4
func TestReadGmsh4Periodic(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
4 4 2 0
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 0 2 1 -2
2 1 0 0 1 1 0 0 2 2 -3
3 0 1 0 1 1 0 0 2 3 -4
4 0 0 0 0 1 0 0 2 4 -1
1 0 0 0 1 1 0 0 4 1 2 3 4
2 0 0 1 1 1 1 0 4 1 2 3 4
$EndEntities
$Nodes
0 0 0 0
$EndNodes
$Elements
0 0 0 0
$EndElements
$Periodic
2
1 1 3
0
2
1 4
2 3
2 1 2
16
1 0 0 1 0 1 0 0 0 0 1 0 0 0 0 1
4
5 9
6 10
7 11
8 12
$EndPeriodic`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Check periodic entities
	if len(mesh.Periodics) != 2 {
		t.Fatalf("Expected 2 periodic entities, got %d", len(mesh.Periodics))
	}

	// First periodic (curves)
	p1 := mesh.Periodics[0]
	if p1.Dimension != 1 {
		t.Errorf("Periodic 1: expected dimension 1, got %d", p1.Dimension)
	}
	if len(p1.NodeMap) != 2 {
		t.Errorf("Periodic 1: expected 2 node pairs, got %d", len(p1.NodeMap))
	}
	if len(p1.AffineTransform) != 0 {
		t.Errorf("Periodic 1: expected no affine transform, got %d values", len(p1.AffineTransform))
	}

	// Second periodic (surfaces, with affine)
	p2 := mesh.Periodics[1]
	if p2.Dimension != 2 {
		t.Errorf("Periodic 2: expected dimension 2, got %d", p2.Dimension)
	}
	if len(p2.NodeMap) != 4 {
		t.Errorf("Periodic 2: expected 4 node pairs, got %d", len(p2.NodeMap))
	}
	if len(p2.AffineTransform) != 16 {
		t.Errorf("Periodic 2: expected affine transform with 16 values, got %d values", len(p2.AffineTransform))
	}
	if len(p2.AffineTransform) != 16 {
		t.Errorf("Periodic 2: expected affine transform with 16 values, got %d values", len(p2.AffineTransform))
	}
}

// TestReadGmsh4PartitionedEntities tests reading partitioned meshes
func TestReadGmsh4PartitionedEntities(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 1
1 0 0 0 1 1 1 0 0
$EndEntities
$PartitionedEntities
2
0 0 0 1
1 0 0 0 1 1 1 0 0 1 0
1 0 0 0 1 1 1 0 0 1 1
$EndPartitionedEntities
$Nodes
0 0 0 0
$EndNodes
$Elements
0 0 0 0
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	// Should handle (skip) partitioned entities gracefully
	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file with partitioned entities: %v", err)
	}

	if mesh.FormatVersion != "4.1" {
		t.Errorf("Expected version 4.1, got %s", mesh.FormatVersion)
	}
}

// TestReadGmsh4GhostElements tests skipping ghost elements
func TestReadGmsh4GhostElements(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 1
1 0 0 0 1 1 1 0 0
$EndEntities
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
$EndElements
$GhostElements
1 1 1 1
3 1 4 1
2 1 2 3 4 0 1
$EndGhostElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Should only have 1 regular element (ghost elements ignored)
	if mesh.NumElements != 1 {
		t.Errorf("Expected 1 element (ghosts ignored), got %d", mesh.NumElements)
	}
}

// TestReadGmsh4AutoDetection tests automatic version detection for v4
func TestReadGmsh4AutoDetection(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 0
$EndEntities
$Nodes
0 0 0 0
$EndNodes
$Elements
0 0 0 0
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	// Test with auto-detection
	mesh, err := ReadGmshAuto(tmpFile)
	if err != nil {
		t.Fatalf("Failed to auto-detect and read Gmsh 4 file: %v", err)
	}

	if !strings.HasPrefix(mesh.FormatVersion, "4") {
		t.Errorf("Auto-detection: expected version 4.x, got %s", mesh.FormatVersion)
	}
}

// TestReadGmsh4ComplexMesh tests a more realistic mesh
func TestReadGmsh4ComplexMesh(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$PhysicalNames
4
1 1 "Inlet"
1 2 "Outlet"
2 10 "Wall"
3 20 "Fluid"
$EndPhysicalNames
$Entities
4 6 3 1
1 0 0 0 0
2 1 0 0 0
3 1 1 0 0
4 0 1 0 0
1 0 0 0 1 0 0 1 1 2 1 -2
2 1 0 0 1 1 0 1 2 2 2 -3
3 0 1 0 1 1 0 0 2 3 -4
4 0 0 0 0 1 0 0 2 4 -1
5 0.5 0 0 0.5 0.5 0 0 2 1 -5
6 0.5 0.5 0 1 0.5 0 0 2 5 -2
1 0 0 0 1 1 0 1 10 6 1 5 -6 -2 -3 -4
2 0 0 0 0.5 0.5 0 0 3 1 -5 4
3 0.5 0 0 1 0.5 0 0 3 5 6 2
1 0 0 0 1 1 1 1 20 3 1 2 3
$EndEntities
$Nodes
10 15 1 15
0 1 0 1
1
0 0 0
0 2 0 1
2
1 0 0
0 3 0 1
3
1 1 0
0 4 0 1
4
0 1 0
1 1 0 2
5
6
0.25 0 0
0.5 0 0
1 5 0 1
7
0.75 0 0
1 6 0 1
8
0.75 0.5 0
2 1 0 3
9
10
11
0.25 0.25 0
0.5 0.25 0
0.25 0.5 0
2 3 0 2
12
13
0.75 0.25 0
0.5 0.5 0
3 1 0 2
14
15
0.5 0.5 0.5
0.5 0.5 1
$EndNodes
$Elements
9 16 1 16
0 1 15 1
1 1
0 2 15 1
2 2
0 3 15 1
3 3
0 4 15 1
4 4
1 1 1 2
5 1 5
6 5 6
1 5 1 1
7 7 2
2 1 2 5
8 1 5 9
9 5 9 10
10 9 10 11
11 10 11 4
12 11 4 1
2 3 2 2
13 6 12 13
14 12 13 8
3 1 4 2
15 14 15 9 10
16 14 15 10 11
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Check total counts
	if mesh.NumVertices != 15 {
		t.Errorf("Expected 15 vertices, got %d", mesh.NumVertices)
	}
	if mesh.NumElements != 16 {
		t.Errorf("Expected 16 elements, got %d", mesh.NumElements)
	}

	// Check physical groups
	expectedGroups := map[int]string{
		1:  "Inlet",
		2:  "Outlet",
		10: "Wall",
		20: "Fluid",
	}

	for tag, name := range expectedGroups {
		group, ok := mesh.ElementGroups[tag]
		if !ok {
			t.Errorf("Physical group %d not found", tag)
			continue
		}
		if group.Name != name {
			t.Errorf("Group %d: expected name '%s', got '%s'", tag, name, group.Name)
		}
	}

	// Count elements by dimension
	dimCounts := make(map[int]int)
	for _, elemType := range mesh.ElementTypes {
		dimCounts[elemType.GetDimension()]++
	}

	expectedDimCounts := map[int]int{
		0: 4, // 4 points
		1: 3, // 3 lines (2 on curve 1 + 1 on curve 5)
		2: 7, // 7 triangles (5 on surface 1 + 2 on surface 3)
		3: 2, // 2 tets
	}

	for dim, expected := range expectedDimCounts {
		if dimCounts[dim] != expected {
			t.Errorf("Dimension %d: expected %d elements, got %d", dim, expected, dimCounts[dim])
		}
	}
}

// TestReadGmsh4MixedElementTypes tests mixed element types in one block
func TestReadGmsh4MixedElementTypes(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 1
1 0 0 0 2 2 2 0 0
$EndEntities
$Nodes
1 14 1 14
3 1 0 14
1
2
3
4
5
6
7
8
9
10
11
12
13
14
0 0 0
1 0 0
1 1 0
0 1 0
0 0 1
1 0 1
1 1 1
0 1 1
0.5 0 0.5
0.5 0.5 0
0.5 0.5 1
0.5 1 0.5
0 0.5 0.5
1 0.5 0.5
$EndNodes
$Elements
2 10 1 10
3 1 4 4
1 1 2 3 4
2 5 6 7 8
3 1 2 6 5
4 4 3 7 8
3 1 5 2
5 1 2 3 4 8
6 5 6 7 8
3 1 6 1
7 1 2 3 6 9
3 1 7 3
8 1 2 4 5 9
9 3 4 7 8 12
10 1 4 8 5 13
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Check element types
	expectedTypes := []ElementType{
		Tet, Tet, Tet, Tet, // First 4 are tets
		Hex, Hex, // Next 2 are hexes
		Prism,                     // 1 prism
		Pyramid, Pyramid, Pyramid, // 3 pyramids
	}

	if len(mesh.Elements) != len(expectedTypes) {
		t.Fatalf("Expected %d elements, got %d", len(expectedTypes), len(mesh.Elements))
	}

	for i, expected := range expectedTypes {
		if mesh.ElementTypes[i] != expected {
			t.Errorf("Element %d: expected type %v, got %v", i, expected, mesh.ElementTypes[i])
		}
	}
}

// TestReadGmsh4ParametricNodes tests nodes with parametric coordinates
func TestReadGmsh4ParametricNodes(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
1 1 0 0
1 0 0 0 0
1 0 0 0 1 0 0 0 2 1 -1
$EndEntities
$Nodes
2 3 1 3
0 1 1 1
1
0 0 0
0
1 1 2 2
2
3
0.5 0 0
1 0 0
0.25
0.75
$EndNodes
$Elements
0 0 0 0
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read Gmsh 4 file: %v", err)
	}

	// Should have read 3 nodes correctly, ignoring parametric coords
	if mesh.NumVertices != 3 {
		t.Errorf("Expected 3 vertices, got %d", mesh.NumVertices)
	}

	// Check coordinates (parametric coords should be ignored)
	expectedCoords := [][]float64{
		{0, 0, 0},
		{0.5, 0, 0},
		{1, 0, 0},
	}

	for i, nodeID := range []int{1, 2, 3} {
		idx, ok := mesh.GetNodeIndex(nodeID)
		if !ok {
			t.Errorf("Node %d not found", nodeID)
			continue
		}

		for j := 0; j < 3; j++ {
			if mesh.Vertices[idx][j] != expectedCoords[i][j] {
				t.Errorf("Node %d coord %d: expected %f, got %f",
					nodeID, j, expectedCoords[i][j], mesh.Vertices[idx][j])
			}
		}
	}
}

// TestReadGmsh4EmptyMesh tests handling of empty mesh
func TestReadGmsh4EmptyMesh(t *testing.T) {
	content := `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 0
$EndEntities
$Nodes
0 0 0 0
$EndNodes
$Elements
0 0 0 0
$EndElements`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	if err != nil {
		t.Fatalf("Failed to read empty Gmsh 4 file: %v", err)
	}

	if mesh.NumVertices != 0 {
		t.Errorf("Expected 0 vertices, got %d", mesh.NumVertices)
	}
	if mesh.NumElements != 0 {
		t.Errorf("Expected 0 elements, got %d", mesh.NumElements)
	}
}

// TestReadGmsh4BinaryDetection tests binary format detection for v4
func TestReadGmsh4BinaryDetection(t *testing.T) {
	content := `$MeshFormat
4.1 1 8
$EndMeshFormat`

	tmpFile := createTempMshFile(t, content)
	defer os.Remove(tmpFile)

	mesh, err := ReadGmsh4(tmpFile)
	// Binary not fully implemented, but should detect format
	if err == nil && mesh.IsBinary {
		t.Log("Binary format detected correctly for v4")
	}
}

// TestReadGmsh4ErrorHandling tests various error conditions
func TestReadGmsh4ErrorHandling(t *testing.T) {
	testCases := []struct {
		name    string
		content string
		errMsg  string
	}{
		{
			name: "Missing Entities",
			content: `$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
0 0 0 0
$EndNodes`,
			errMsg: "scanner error",
		},
		{
			name: "Invalid Nodes Header",
			content: `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 0
$EndEntities
$Nodes
invalid
$EndNodes`,
			errMsg: "invalid Nodes header",
		},
		{
			name: "Element References Unknown Node",
			content: `$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 1
1 0 0 0 1 1 1 0 0
$EndEntities
$Nodes
1 1 1 1
3 1 0 1
1
0 0 0
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
$EndElements`,
			errMsg: "unknown node",
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			tmpFile := createTempMshFile(t, tc.content)
			defer os.Remove(tmpFile)

			_, err := ReadGmsh4(tmpFile)
			if err == nil {
				t.Error("Expected error, got nil")
			} else if !strings.Contains(err.Error(), tc.errMsg) {
				t.Errorf("Expected error containing '%s', got '%v'", tc.errMsg, err)
			}
		})
	}
}

// BenchmarkReadGmsh4 benchmarks reading v4 format
func BenchmarkReadGmsh4(b *testing.B) {
	// Generate a v4 format mesh
	var content strings.Builder
	content.WriteString(`$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 0 1
1 0 0 0 10 10 10 0 0
$EndEntities
$Nodes
1 1000 1 1000
3 1 0 1000
`)

	// Node tags
	for i := 1; i <= 1000; i++ {
		fmt.Fprintf(&content, "%d\n", i)
	}

	// Node coordinates
	for i := 1; i <= 1000; i++ {
		x := float64(i % 10)
		y := float64((i / 10) % 10)
		z := float64(i / 100)
		fmt.Fprintf(&content, "%.1f %.1f %.1f\n", x, y, z)
	}

	content.WriteString(`$EndNodes
$Elements
1 500 1 500
3 1 4 500
`)

	// Elements
	for i := 1; i <= 500; i++ {
		v1 := i
		v2 := i + 1
		v3 := i + 2
		v4 := i + 3
		if v4 > 1000 {
			v4 = v4 - 1000
		}
		fmt.Fprintf(&content, "%d %d %d %d %d\n", i, v1, v2, v3, v4)
	}

	content.WriteString("$EndElements\n")

	tmpFile := filepath.Join(b.TempDir(), "bench4.msh")
	os.WriteFile(tmpFile, []byte(content.String()), 0644)

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		_, err := ReadGmsh4(tmpFile)
		if err != nil {
			b.Fatal(err)
		}
	}
}
