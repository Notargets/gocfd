package utils

import (
	"testing"
)

func TestStandardTestMeshes(t *testing.T) {
	tm := GetStandardTestMeshes()

	t.Run("CubeNodes", func(t *testing.T) {
		// Verify we have 16 nodes
		if len(tm.CubeNodes.Nodes) != 16 {
			t.Errorf("Expected 16 cube nodes, got %d", len(tm.CubeNodes.Nodes))
		}

		// Verify node map consistency
		for name, idx := range tm.CubeNodes.NodeMap {
			if idx < 0 || idx >= len(tm.CubeNodes.Nodes) {
				t.Errorf("Invalid node index %d for %s", idx, name)
			}

			// Check ID map
			expectedID := idx + 1
			if tm.CubeNodes.NodeIDMap[name] != expectedID {
				t.Errorf("Node %s: ID map has %d, expected %d",
					name, tm.CubeNodes.NodeIDMap[name], expectedID)
			}
		}

		// Verify specific coordinates
		origin := tm.CubeNodes.Nodes[tm.CubeNodes.NodeMap["origin"]]
		if origin[0] != 0 || origin[1] != 0 || origin[2] != 0 {
			t.Errorf("Origin should be at (0,0,0), got (%f,%f,%f)",
				origin[0], origin[1], origin[2])
		}

		xyz := tm.CubeNodes.Nodes[tm.CubeNodes.NodeMap["xyz"]]
		if xyz[0] != 1 || xyz[1] != 1 || xyz[2] != 1 {
			t.Errorf("XYZ corner should be at (1,1,1), got (%f,%f,%f)",
				xyz[0], xyz[1], xyz[2])
		}
	})

	t.Run("ElementSets", func(t *testing.T) {
		// Test single tet
		if tm.SingleTet.Type != Tet {
			t.Errorf("SingleTet should be type Tet, got %v", tm.SingleTet.Type)
		}
		if len(tm.SingleTet.Elements) != 1 {
			t.Errorf("SingleTet should have 1 element, got %d", len(tm.SingleTet.Elements))
		}
		if len(tm.SingleTet.Elements[0]) != 4 {
			t.Errorf("Tet should have 4 nodes, got %d", len(tm.SingleTet.Elements[0]))
		}

		// Test single hex
		if tm.SingleHex.Type != Hex {
			t.Errorf("SingleHex should be type Hex, got %v", tm.SingleHex.Type)
		}
		if len(tm.SingleHex.Elements[0]) != 8 {
			t.Errorf("Hex should have 8 nodes, got %d", len(tm.SingleHex.Elements[0]))
		}

		// Test single prism
		if tm.SinglePrism.Type != Prism {
			t.Errorf("SinglePrism should be type Prism, got %v", tm.SinglePrism.Type)
		}
		if len(tm.SinglePrism.Elements[0]) != 6 {
			t.Errorf("Prism should have 6 nodes, got %d", len(tm.SinglePrism.Elements[0]))
		}

		// Test single pyramid
		if tm.SinglePyramid.Type != Pyramid {
			t.Errorf("SinglePyramid should be type Pyramid, got %v", tm.SinglePyramid.Type)
		}
		if len(tm.SinglePyramid.Elements[0]) != 5 {
			t.Errorf("Pyramid should have 5 nodes, got %d", len(tm.SinglePyramid.Elements[0]))
		}
	})

	t.Run("CompleteMeshes", func(t *testing.T) {
		// Test TwoTetMesh
		if len(tm.TwoTetMesh.Elements) != 1 {
			t.Errorf("TwoTetMesh should have 1 element set, got %d",
				len(tm.TwoTetMesh.Elements))
		}
		if len(tm.TwoTetMesh.Elements[0].Elements) != 2 {
			t.Errorf("TwoTetMesh should have 2 tets, got %d",
				len(tm.TwoTetMesh.Elements[0].Elements))
		}

		// Test MixedMesh
		if len(tm.MixedMesh.Elements) != 4 {
			t.Errorf("MixedMesh should have 4 element sets, got %d",
				len(tm.MixedMesh.Elements))
		}

		expectedTypes := []ElementType{Tet, Hex, Prism, Pyramid}
		for i, expectedType := range expectedTypes {
			if tm.MixedMesh.Elements[i].Type != expectedType {
				t.Errorf("MixedMesh element set %d: expected type %v, got %v",
					i, expectedType, tm.MixedMesh.Elements[i].Type)
			}
		}

		// Test CubeMesh
		if len(tm.CubeMesh.Elements) != 1 {
			t.Errorf("CubeMesh should have 1 element set, got %d",
				len(tm.CubeMesh.Elements))
		}
		if len(tm.CubeMesh.Elements[0].Elements) != 6 {
			t.Errorf("CubeMesh should have 6 tets, got %d",
				len(tm.CubeMesh.Elements[0].Elements))
		}
	})
}

func TestValidationHelpers(t *testing.T) {
	t.Run("ValidateNodeCoordinates", func(t *testing.T) {
		nodes := [][]float64{
			{0, 0, 0},
			{1, 0, 0},
			{0, 1, 0},
		}

		// Exact match should pass
		err := ValidateNodeCoordinates(nodes, nodes, 1e-10)
		if err != nil {
			t.Errorf("Exact match failed: %v", err)
		}

		// Small difference within tolerance should pass
		nodesWithSmallDiff := [][]float64{
			{0, 0, 0},
			{1.0000001, 0, 0},
			{0, 1, 0},
		}
		err = ValidateNodeCoordinates(nodes, nodesWithSmallDiff, 1e-6)
		if err != nil {
			t.Errorf("Small difference within tolerance failed: %v", err)
		}

		// Large difference should fail
		nodesWithLargeDiff := [][]float64{
			{0, 0, 0},
			{1.1, 0, 0},
			{0, 1, 0},
		}
		err = ValidateNodeCoordinates(nodes, nodesWithLargeDiff, 1e-6)
		if err == nil {
			t.Error("Large difference should have failed validation")
		}

		// Different count should fail
		fewerNodes := [][]float64{
			{0, 0, 0},
			{1, 0, 0},
		}
		err = ValidateNodeCoordinates(nodes, fewerNodes, 1e-6)
		if err == nil {
			t.Error("Different node count should have failed validation")
		}
	})

	t.Run("ValidateElementConnectivity", func(t *testing.T) {
		elements := [][]int{
			{0, 1, 2, 3},
			{1, 2, 3, 4},
		}

		// Exact match should pass
		err := ValidateElementConnectivity(elements, elements)
		if err != nil {
			t.Errorf("Exact match failed: %v", err)
		}

		// Different connectivity should fail
		differentElements := [][]int{
			{0, 1, 2, 3},
			{1, 2, 3, 5}, // Changed last node
		}
		err = ValidateElementConnectivity(elements, differentElements)
		if err == nil {
			t.Error("Different connectivity should have failed validation")
		}

		// Different node count should fail
		differentNodeCount := [][]int{
			{0, 1, 2, 3},
			{1, 2, 3}, // Missing a node
		}
		err = ValidateElementConnectivity(elements, differentNodeCount)
		if err == nil {
			t.Error("Different node count should have failed validation")
		}

		// Different element count should fail
		fewerElements := [][]int{
			{0, 1, 2, 3},
		}
		err = ValidateElementConnectivity(elements, fewerElements)
		if err == nil {
			t.Error("Different element count should have failed validation")
		}
	})
}
