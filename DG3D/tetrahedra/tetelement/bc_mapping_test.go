package tetelement

import (
	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/utils"
	"testing"
)

// TestBCMapping_SingleTetAllBoundary tests BC mapping for a single tet with all boundary faces
// Following Unit Testing Principle: Start with fundamentals
func TestBCMapping_SingleTetAllBoundary(t *testing.T) {
	// Use test mesh helpers
	tm := utils.GetStandardTestMeshes()

	// Create a single tet mesh
	singleTetMesh := utils.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []utils.ElementSet{tm.SingleTet},
		Dimension: 3,
	}

	// Convert to mesh
	m := mesh.ConvertToMesh(singleTetMesh)

	// Add boundary elements for all faces
	m.BoundaryElements = map[string][]mesh.BoundaryElement{
		"inlet": {
			{ParentElement: 0, ParentFace: 0, ElementType: utils.Triangle},
		},
		"outlet": {
			{ParentElement: 0, ParentFace: 1, ElementType: utils.Triangle},
		},
		"wall": {
			{ParentElement: 0, ParentFace: 2, ElementType: utils.Triangle},
			{ParentElement: 0, ParentFace: 3, ElementType: utils.Triangle},
		},
	}

	// Create Element3D
	order := 2
	el, err := NewElement3DFromMesh(order, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Test 1: Verify BC maps were created
	if el.BCMaps == nil {
		t.Fatal("BC maps not initialized")
	}

	// Test 2: Check face BC mapping
	testCases := []struct {
		face     int
		expected utils.BCType
	}{
		{0, utils.BCInflow},
		{1, utils.BCOutflow},
		{2, utils.BCWall},
		{3, utils.BCWall},
	}

	for _, tc := range testCases {
		bcType := el.GetFaceBCType(0, tc.face)
		if bcType != tc.expected {
			t.Errorf("Face %d: expected BC type %s, got %s",
				tc.face, tc.expected, bcType)
		}
	}

	// Test 3: Verify all boundary nodes have BC assignments
	if len(el.BCMaps.NodeBC) == 0 {
		t.Error("No boundary nodes have BC assignments")
	}

	// Test 4: Check boundary statistics
	stats := el.GetBoundaryStatistics()
	expectedStats := map[utils.BCType]int{
		utils.BCInflow:  el.Nfp,     // One face
		utils.BCOutflow: el.Nfp,     // One face
		utils.BCWall:    2 * el.Nfp, // Two faces
	}

	for bcType, count := range expectedStats {
		if stats[bcType] != count {
			t.Errorf("BC type %s: expected %d nodes, got %d",
				bcType, count, stats[bcType])
		}
	}
}

// TestBCMapping_TwoTets tests BC mapping with interior and boundary faces
// Following Unit Testing Principle: Build systematically
func TestBCMapping_TwoTets(t *testing.T) {
	// Use test mesh helpers
	tm := utils.GetStandardTestMeshes()

	// Convert to mesh
	m := mesh.ConvertToMesh(tm.TwoTetMesh)

	// Add boundary elements (only exterior faces)
	m.BoundaryElements = map[string][]mesh.BoundaryElement{
		"wall": {
			// All exterior faces are walls
			{ParentElement: 0, ParentFace: 0, ElementType: utils.Triangle},
			{ParentElement: 0, ParentFace: 1, ElementType: utils.Triangle},
			{ParentElement: 0, ParentFace: 3, ElementType: utils.Triangle},
			{ParentElement: 1, ParentFace: 0, ElementType: utils.Triangle},
			{ParentElement: 1, ParentFace: 2, ElementType: utils.Triangle},
			{ParentElement: 1, ParentFace: 3, ElementType: utils.Triangle},
		},
	}

	// Create Element3D
	order := 1
	el, err := NewElement3DFromMesh(order, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Test 1: Interior faces should have BCNone
	// Face 2 of element 0 and face 1 of element 1 are shared
	if bc := el.GetFaceBCType(0, 2); bc != utils.BCNone {
		t.Errorf("Interior face should have BCNone, got %s", bc)
	}
	if bc := el.GetFaceBCType(1, 1); bc != utils.BCNone {
		t.Errorf("Interior face should have BCNone, got %s", bc)
	}

	// Test 2: Boundary faces should have BCWall
	boundaryFaces := []struct{ elem, face int }{
		{0, 0}, {0, 1}, {0, 3},
		{1, 0}, {1, 2}, {1, 3},
	}

	for _, bf := range boundaryFaces {
		if bc := el.GetFaceBCType(bf.elem, bf.face); bc != utils.BCWall {
			t.Errorf("Boundary face (%d,%d) should have BCWall, got %s",
				bf.elem, bf.face, bc)
		}
	}

	// Test 3: Number of boundary nodes should match MapB length
	// Note: MapB includes duplicate entries for nodes shared between faces
	// (e.g., nodes at edges and vertices), so we can't simply calculate 6 * Nfp
	stats := el.GetBoundaryStatistics()
	totalBoundaryNodes := 0
	for _, count := range stats {
		totalBoundaryNodes += count
	}

	// The total should match the length of MapB (which includes duplicates for shared nodes)
	if totalBoundaryNodes != len(el.DG3D.MapB) {
		t.Errorf("Total boundary nodes %d should match len(MapB) %d",
			totalBoundaryNodes, len(el.DG3D.MapB))
	}

	// All boundary nodes should be BCWall
	if wallCount, ok := stats[utils.BCWall]; !ok || wallCount != totalBoundaryNodes {
		t.Errorf("All %d boundary nodes should be BCWall", totalBoundaryNodes)
	}
}

// TestBCMapping_MixedBCTypes tests multiple BC types
// Following Unit Testing Principle: Test specific properties
func TestBCMapping_MixedBCTypes(t *testing.T) {
	// Use the two tet mesh for testing mixed BC types
	tm := utils.GetStandardTestMeshes()

	// Convert to mesh
	m := mesh.ConvertToMesh(tm.TwoTetMesh)

	// Mix of BC types - assign to ALL faces (interior ones will be ignored)
	m.BoundaryElements = map[string][]mesh.BoundaryElement{
		"velocity_inlet": {
			{ParentElement: 0, ParentFace: 0, ElementType: utils.Triangle},
		},
		"pressure_outlet": {
			{ParentElement: 1, ParentFace: 3, ElementType: utils.Triangle},
		},
		"slip_wall": {
			{ParentElement: 0, ParentFace: 1, ElementType: utils.Triangle},
			{ParentElement: 1, ParentFace: 0, ElementType: utils.Triangle},
		},
		"symmetry": {
			{ParentElement: 0, ParentFace: 3, ElementType: utils.Triangle},
			{ParentElement: 1, ParentFace: 2, ElementType: utils.Triangle},
		},
		"adiabatic": {
			{ParentElement: 0, ParentFace: 2, ElementType: utils.Triangle}, // This might be interior
			{ParentElement: 1, ParentFace: 1, ElementType: utils.Triangle}, // This might be interior
		},
	}

	order := 2
	el, err := NewElement3DFromMesh(order, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Test each BC type is correctly assigned
	testCases := []struct {
		elem     int
		face     int
		expected utils.BCType
	}{
		{0, 0, utils.BCVelocityInlet},
		{1, 3, utils.BCPressureOutlet},
		{0, 1, utils.BCSlipWall},
		{1, 0, utils.BCSlipWall},
		{0, 3, utils.BCSymmetry},
		{1, 2, utils.BCSymmetry},
		{0, 2, utils.BCAdiabatic},
		{1, 1, utils.BCAdiabatic},
	}

	// Add this loop to actually use testCases:
	for _, tc := range testCases {
		bcType := el.GetFaceBCType(tc.elem, tc.face)
		if bcType != tc.expected {
			t.Errorf("Element %d Face %d: expected %s, got %s",
				tc.elem, tc.face, tc.expected, bcType)
		}
	}

	// Verify statistics - we assigned 5 different BC types
	stats := el.GetBoundaryStatistics()
	if len(stats) != 5 {
		t.Errorf("Expected 5 different BC types, got %d", len(stats))
	}
}

// TestBCMapping_UnknownBCNames tests handling of unknown BC names
// Following Unit Testing Principle: Test edge cases
func TestBCMapping_UnknownBCNames(t *testing.T) {
	// Use test mesh helpers
	tm := utils.GetStandardTestMeshes()

	// Create a single tet mesh
	singleTetMesh := utils.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []utils.ElementSet{tm.SingleTet},
		Dimension: 3,
	}

	// Convert to mesh
	m := mesh.ConvertToMesh(singleTetMesh)

	// Use unknown BC names
	m.BoundaryElements = map[string][]mesh.BoundaryElement{
		"unknown_bc_1": {
			{ParentElement: 0, ParentFace: 0, ElementType: utils.Triangle},
		},
		"another_unknown": {
			{ParentElement: 0, ParentFace: 1, ElementType: utils.Triangle},
		},
	}

	order := 1
	el, err := NewElement3DFromMesh(order, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Unknown BC names should default to BCWall
	for face := 0; face < 2; face++ {
		bcType := el.GetFaceBCType(0, face)
		if bcType != utils.BCWall {
			t.Errorf("Unknown BC name should default to BCWall, got %s", bcType)
		}
	}
}
