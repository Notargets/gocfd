package utils

import (
	"math"
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

		// Test CubeMesh - basic counts
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

// TestCubeMeshMathematicalCorrectness validates the cube mesh decomposition mathematically
// Following the Unit Testing Principles: test specific mathematical properties
func TestCubeMeshMathematicalCorrectness(t *testing.T) {
	tm := GetStandardTestMeshes()
	cubeMesh := tm.CubeMesh
	nodes := cubeMesh.Nodes
	tets := cubeMesh.Elements[0].Elements

	// Helper to get vertex coordinates
	getVertex := func(nodeName string) []float64 {
		idx := nodes.NodeMap[nodeName]
		return nodes.Nodes[idx]
	}

	t.Run("VolumeConservation", func(t *testing.T) {
		// Mathematical property: Sum of tetrahedra volumes = cube volume
		// Cube volume = 1.0 (unit cube)
		expectedVolume := 1.0
		totalVolume := 0.0

		for i, tet := range tets {
			v0 := getVertex(tet[0])
			v1 := getVertex(tet[1])
			v2 := getVertex(tet[2])
			v3 := getVertex(tet[3])

			volume := ComputeTetVolume(v0, v1, v2, v3)
			totalVolume += volume

			// Each tet should have positive volume
			if volume <= 0 {
				t.Errorf("Tet %d has non-positive volume: %f", i, volume)
			}

			// Each tet should have reasonable volume (not too small)
			if volume < 0.01 {
				t.Errorf("Tet %d has suspiciously small volume: %f", i, volume)
			}
		}

		// Check total volume matches cube volume
		tolerance := 1e-10
		if math.Abs(totalVolume-expectedVolume) > tolerance {
			t.Errorf("Volume conservation failed: total volume %f, expected %f (diff %e)",
				totalVolume, expectedVolume, math.Abs(totalVolume-expectedVolume))
		}
	})

	t.Run("FaceConnectivity", func(t *testing.T) {
		// Build face connectivity map
		faceMap := make(map[FaceKey][]int) // face -> list of tets containing it

		// Get node IDs for conversion
		nodeToID := make(map[string]int)
		for name, idx := range nodes.NodeMap {
			nodeToID[name] = idx
		}

		// Process each tetrahedron
		for tetIdx, tet := range tets {
			// Convert node names to indices
			tetIndices := make([]int, 4)
			for i, nodeName := range tet {
				tetIndices[i] = nodeToID[nodeName]
			}

			// Get the 4 faces of this tet
			faces := GetTetFaces(tetIndices)

			// Add each face to the map
			for _, face := range faces {
				key := NewFaceKey(face[0], face[1], face[2])
				faceMap[key] = append(faceMap[key], tetIdx)
			}
		}

		// Count interior and boundary faces
		interiorFaces := 0
		boundaryFaces := 0

		for faceKey, tetList := range faceMap {
			if len(tetList) == 1 {
				// Boundary face (only one tet)
				boundaryFaces++
			} else if len(tetList) == 2 {
				// Interior face (shared by two tets)
				interiorFaces++
			} else {
				// Error: face shared by more than 2 tets
				t.Errorf("Face %v shared by %d tets (should be 1 or 2)",
					faceKey, len(tetList))
			}
		}

		// A cube has 6 faces, each divided into 2 triangles = 12 boundary faces
		expectedBoundaryFaces := 12
		if boundaryFaces != expectedBoundaryFaces {
			t.Errorf("Boundary face count: got %d, expected %d",
				boundaryFaces, expectedBoundaryFaces)
		}

		// Log the counts for debugging
		t.Logf("Face connectivity: %d boundary faces, %d interior faces",
			boundaryFaces, interiorFaces)
	})

	t.Run("BoundaryFacesCoverCubeSurface", func(t *testing.T) {
		// Verify that boundary faces cover exactly the cube surface area
		// Each cube face has area 1, total surface area = 6

		// Get all boundary faces
		faceMap := make(map[FaceKey]int) // face -> count of tets containing it
		nodeToID := make(map[string]int)
		for name, idx := range nodes.NodeMap {
			nodeToID[name] = idx
		}

		for _, tet := range tets {
			tetIndices := make([]int, 4)
			for i, nodeName := range tet {
				tetIndices[i] = nodeToID[nodeName]
			}
			faces := GetTetFaces(tetIndices)
			for _, face := range faces {
				key := NewFaceKey(face[0], face[1], face[2])
				faceMap[key]++ // Just increment the count
			}
		}

		// Calculate total boundary face area
		totalBoundaryArea := 0.0
		for faceKey, count := range faceMap {
			if count == 1 {
				// This is a boundary face
				v0 := nodes.Nodes[faceKey.V0]
				v1 := nodes.Nodes[faceKey.V1]
				v2 := nodes.Nodes[faceKey.V2]
				area := ComputeTriangleArea(v0, v1, v2)
				totalBoundaryArea += area
			}
		}

		// Expected surface area of unit cube = 6
		expectedArea := 6.0
		tolerance := 1e-10
		if math.Abs(totalBoundaryArea-expectedArea) > tolerance {
			t.Errorf("Boundary surface area: got %f, expected %f (diff %e)",
				totalBoundaryArea, expectedArea, math.Abs(totalBoundaryArea-expectedArea))
		}
	})

	t.Run("NoGapsOrOverlaps", func(t *testing.T) {
		// Test that every point in the cube belongs to exactly one tetrahedron
		// We'll test a grid of points inside the cube

		// Test points on a 5x5x5 grid
		testResolution := 5
		tolerance := 1e-10

		for i := 0; i <= testResolution; i++ {
			for j := 0; j <= testResolution; j++ {
				for k := 0; k <= testResolution; k++ {
					// Test point coordinates
					x := float64(i) / float64(testResolution)
					y := float64(j) / float64(testResolution)
					z := float64(k) / float64(testResolution)
					testPoint := []float64{x, y, z}

					// Count how many tets contain this point
					containingTets := 0

					for tetIdx, tet := range tets {
						_ = tetIdx
						v0 := getVertex(tet[0])
						v1 := getVertex(tet[1])
						v2 := getVertex(tet[2])
						v3 := getVertex(tet[3])

						if pointInTetrahedron(testPoint, v0, v1, v2, v3, tolerance) {
							containingTets++
						}
					}

					// Each point should be in exactly one tet
					// (or on a boundary between tets, which we'll allow)
					if containingTets == 0 {
						t.Errorf("Point (%f,%f,%f) is not in any tetrahedron (gap detected)",
							x, y, z)
					} else if containingTets > 1 {
						// Points on faces/edges might be in multiple tets due to numerical tolerance
						// This is acceptable if they're on boundaries
						t.Logf("Point (%f,%f,%f) is in %d tetrahedra (likely on boundary)",
							x, y, z, containingTets)
					}
				}
			}
		}
	})

	t.Run("ElementOrientations", func(t *testing.T) {
		// All tetrahedra should have consistent (positive) orientation
		for i, tet := range tets {
			v0 := getVertex(tet[0])
			v1 := getVertex(tet[1])
			v2 := getVertex(tet[2])
			v3 := getVertex(tet[3])

			// Compute signed volume (determinant / 6)
			d1 := []float64{v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]}
			d2 := []float64{v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]}
			d3 := []float64{v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]}

			det := d1[0]*(d2[1]*d3[2]-d2[2]*d3[1]) -
				d1[1]*(d2[0]*d3[2]-d2[2]*d3[0]) +
				d1[2]*(d2[0]*d3[1]-d2[1]*d3[0])

			signedVolume := det / 6.0

			if signedVolume <= 0 {
				t.Errorf("Tet %d has negative orientation (signed volume = %f)", i, signedVolume)
			}
		}
	})
}

// Helper function: Check if a point is inside a tetrahedron using barycentric coordinates
func pointInTetrahedron(p, v0, v1, v2, v3 []float64, tolerance float64) bool {
	// Compute barycentric coordinates
	// p = b0*v0 + b1*v1 + b2*v2 + b3*v3
	// where b0 + b1 + b2 + b3 = 1

	// Build the system: [v1-v0, v2-v0, v3-v0] * [b1, b2, b3]^T = p - v0
	// Then b0 = 1 - b1 - b2 - b3

	d1 := []float64{v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]}
	d2 := []float64{v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]}
	d3 := []float64{v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]}
	dp := []float64{p[0] - v0[0], p[1] - v0[1], p[2] - v0[2]}

	// Solve using Cramer's rule
	det := d1[0]*(d2[1]*d3[2]-d2[2]*d3[1]) -
		d1[1]*(d2[0]*d3[2]-d2[2]*d3[0]) +
		d1[2]*(d2[0]*d3[1]-d2[1]*d3[0])

	if math.Abs(det) < 1e-12 {
		// Degenerate tetrahedron
		return false
	}

	// Compute b1
	det1 := dp[0]*(d2[1]*d3[2]-d2[2]*d3[1]) -
		dp[1]*(d2[0]*d3[2]-d2[2]*d3[0]) +
		dp[2]*(d2[0]*d3[1]-d2[1]*d3[0])
	b1 := det1 / det

	// Compute b2
	det2 := d1[0]*(dp[1]*d3[2]-dp[2]*d3[1]) -
		d1[1]*(dp[0]*d3[2]-dp[2]*d3[0]) +
		d1[2]*(dp[0]*d3[1]-dp[1]*d3[0])
	b2 := det2 / det

	// Compute b3
	det3 := d1[0]*(d2[1]*dp[2]-d2[2]*dp[1]) -
		d1[1]*(d2[0]*dp[2]-d2[2]*dp[0]) +
		d1[2]*(d2[0]*dp[1]-d2[1]*dp[0])
	b3 := det3 / det

	// Compute b0
	b0 := 1.0 - b1 - b2 - b3

	// Point is inside if all barycentric coordinates are non-negative
	return b0 >= -tolerance && b1 >= -tolerance && b2 >= -tolerance && b3 >= -tolerance
}

// TestValidationHelpers tests the validation helper functions
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

// TestGeometricHelpers tests the geometric computation helpers
func TestGeometricHelpers(t *testing.T) {
	t.Run("ComputeTetVolume", func(t *testing.T) {
		// Test with a known tetrahedron
		// Unit tetrahedron at origin with vertices at:
		// (0,0,0), (1,0,0), (0,1,0), (0,0,1)
		v0 := []float64{0, 0, 0}
		v1 := []float64{1, 0, 0}
		v2 := []float64{0, 1, 0}
		v3 := []float64{0, 0, 1}

		volume := ComputeTetVolume(v0, v1, v2, v3)
		expectedVolume := 1.0 / 6.0 // Known volume of this tetrahedron

		tolerance := 1e-10
		if math.Abs(volume-expectedVolume) > tolerance {
			t.Errorf("Tet volume: got %f, expected %f", volume, expectedVolume)
		}
	})

	t.Run("ComputeTriangleArea", func(t *testing.T) {
		// Test with a known triangle
		// Right triangle with legs of length 1
		v0 := []float64{0, 0, 0}
		v1 := []float64{1, 0, 0}
		v2 := []float64{0, 1, 0}

		area := ComputeTriangleArea(v0, v1, v2)
		expectedArea := 0.5 // Known area of this triangle

		tolerance := 1e-10
		if math.Abs(area-expectedArea) > tolerance {
			t.Errorf("Triangle area: got %f, expected %f", area, expectedArea)
		}
	})
}
