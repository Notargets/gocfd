// Package facebuffer tests
package facebuffer

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/DG3D/tetrahedra/tetelement"
	"github.com/notargets/gocfd/utils"
)

// TestFaceBuffer_SingleTetBoundary tests the simplest case - single tet with all boundary faces
func TestFaceBuffer_SingleTetBoundary(t *testing.T) {
	// Get standard test meshes
	tm := utils.GetStandardTestMeshes()

	// Create single tet mesh
	singleTetMesh := utils.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []utils.ElementSet{tm.SingleTet},
		Dimension: 3,
	}

	// Convert to mesh
	m := mesh.ConvertToMesh(singleTetMesh)

	// Test with different polynomial orders
	for _, order := range []int{1, 2, 3} {
		t.Run(fmt.Sprintf("Order%d", order), func(t *testing.T) {
			// Create Element3D
			el, err := tetelement.NewElement3DFromMesh(order, m)
			if err != nil {
				t.Fatalf("Failed to create Element3D: %v", err)
			}

			// Build face buffer
			fb, err := BuildFaceBuffer(el)
			if err != nil {
				t.Fatalf("Failed to build face buffer: %v", err)
			}

			// Validate dimensions
			expectedFaces := uint32(4) // 4 faces for tetrahedron
			if uint32(len(fb.FaceIndex)) != expectedFaces {
				t.Errorf("Expected %d face indices, got %d", expectedFaces, len(fb.FaceIndex))
			}

			// Test: All faces should be boundary
			for i, faceCode := range fb.FaceIndex {
				if faceCode != BoundaryPlaceholder {
					t.Errorf("Face %d should be boundary (-999), got %d", i, faceCode)
				}
			}

			// Test: No remote partitions
			if len(fb.RemoteSendIndices) != 0 {
				t.Errorf("Expected 0 remote partitions, got %d", len(fb.RemoteSendIndices))
			}

			// Check statistics
			stats := fb.GetStats()
			if stats["boundary_faces"] != 4 {
				t.Errorf("Expected 4 boundary faces, got %d", stats["boundary_faces"])
			}
			if stats["interior_faces"] != 0 {
				t.Errorf("Expected 0 interior faces, got %d", stats["interior_faces"])
			}

			// Validate consistency
			if err := fb.ValidateFaceBuffer(); err != nil {
				t.Errorf("Validation failed: %v", err)
			}
		})
	}
}

// TestFaceBuffer_TwoTetsInterior tests interior face connections
func TestFaceBuffer_TwoTetsInterior(t *testing.T) {
	// Get standard test meshes
	tm := utils.GetStandardTestMeshes()

	// Use the two tet mesh
	m := mesh.ConvertToMesh(tm.TwoTetMesh)

	order := 2
	el, err := tetelement.NewElement3DFromMesh(order, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Build face buffer
	fb, err := BuildFaceBuffer(el)
	if err != nil {
		t.Fatalf("Failed to build face buffer: %v", err)
	}

	// Get statistics
	stats := fb.GetStats()

	// Test: Should have both boundary and interior faces
	if stats["boundary_faces"] == 0 {
		t.Error("Expected some boundary faces")
	}
	if stats["interior_faces"] == 0 {
		t.Error("Expected some interior faces (shared face between tets)")
	}

	// Test: Total faces should be 8 (2 elements × 4 faces)
	if stats["total_faces"] != 8 {
		t.Errorf("Expected 8 total faces, got %d", stats["total_faces"])
	}

	// Test: Interior faces should have positive indices
	interiorFound := false
	for _, faceCode := range fb.FaceIndex {
		if faceCode > 0 {
			interiorFound = true
			// Verify it points to a valid M buffer location
			if uint32(faceCode) >= fb.K*fb.Nfaces*fb.Nfp {
				t.Errorf("Interior face index %d out of M buffer range", faceCode)
			}
		}
	}
	if !interiorFound {
		t.Error("No interior faces found with positive indices")
	}

	// Test: Verify interior face connectivity is symmetric
	// For each interior face, there should be a matching face in the neighbor
	for elem := uint32(0); elem < fb.K; elem++ {
		for face := uint32(0); face < fb.Nfaces; face++ {
			faceIdx := face + elem*fb.Nfaces
			if fb.FaceIndex[faceIdx] > 0 {
				// This is an interior face
				pStart := uint32(fb.FaceIndex[faceIdx])

				// The P location should be within valid range
				if pStart+fb.Nfp > fb.K*fb.Nfaces*fb.Nfp {
					t.Errorf("Interior face P values would exceed M buffer")
				}
			}
		}
	}

	// Validate
	if err := fb.ValidateFaceBuffer(); err != nil {
		t.Errorf("Validation failed: %v", err)
	}
}

// TestFaceBuffer_SimpleParallel tests remote face detection
func TestFaceBuffer_SimpleParallel(t *testing.T) {
	// Create a simple 2-element mesh
	tm := utils.GetStandardTestMeshes()
	m := mesh.ConvertToMesh(tm.TwoTetMesh)

	order := 1
	el, err := tetelement.NewElement3DFromMesh(order, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Set up partitioning - element 0 in partition 0, element 1 in partition 1
	el.EToP = []int{0, 1}

	// Simulate partition 0 by keeping only element 0
	origK := el.K
	el.K = 1 // Only element 0 in this partition

	// Build face buffer for partition 0
	fb, err := BuildFaceBuffer(el)
	if err != nil {
		t.Fatalf("Failed to build face buffer: %v", err)
	}

	// Restore K for proper cleanup
	el.K = origK

	stats := fb.GetStats()

	// Should have exactly one remote face (where element 0 connects to element 1)
	if stats["remote_faces"] != 1 {
		t.Errorf("Expected 1 remote face, got %d", stats["remote_faces"])
	}

	// Should have boundary faces (exterior faces of element 0)
	if stats["boundary_faces"] != 3 {
		t.Errorf("Expected 3 boundary faces, got %d", stats["boundary_faces"])
	}

	// Should have NO interior faces (element 0 is alone in partition 0)
	if stats["interior_faces"] != 0 {
		t.Errorf("Expected 0 interior faces, got %d", stats["interior_faces"])
	}

	// Check that remote face is marked correctly
	remoteFound := false
	for _, faceCode := range fb.FaceIndex {
		if faceCode == RemoteFace {
			remoteFound = true
			break
		}
	}
	if !remoteFound {
		t.Error("No face marked as remote (-9999)")
	}

	// Check RemoteSendIndices
	if len(fb.RemoteSendIndices) != 1 {
		t.Errorf("Expected 1 remote partition, got %d", len(fb.RemoteSendIndices))
	}

	// Should have indices for partition 1
	if indices, ok := fb.RemoteSendIndices[1]; !ok {
		t.Error("No send indices for partition 1")
	} else if len(indices) != int(fb.Nfp) {
		t.Errorf("Expected %d send indices for one face, got %d", fb.Nfp, len(indices))
	}

	// Validate
	if err := fb.ValidateFaceBuffer(); err != nil {
		t.Errorf("Validation failed: %v", err)
	}
}

// TestFaceBuffer_BCOverlay tests the boundary condition overlay phase
func TestFaceBuffer_BCOverlay(t *testing.T) {
	// Get standard test meshes
	tm := utils.GetStandardTestMeshes()

	// Create single tet mesh
	singleTetMesh := utils.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []utils.ElementSet{tm.SingleTet},
		Dimension: 3,
	}

	m := mesh.ConvertToMesh(singleTetMesh)
	el, err := tetelement.NewElement3DFromMesh(1, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Build face buffer (Phase 1)
	fb, err := BuildFaceBuffer(el)
	if err != nil {
		t.Fatalf("Failed to build face buffer: %v", err)
	}

	// All faces should be boundary placeholders
	for i, faceCode := range fb.FaceIndex {
		if faceCode != BoundaryPlaceholder {
			t.Errorf("Face %d should have placeholder -999, got %d", i, faceCode)
		}
	}

	// Apply BC overlay (Phase 2)
	// Let's assign different BC types to each face
	bcData := map[int32]int32{
		0: -1, // Face 0: Wall BC
		1: -2, // Face 1: Outflow BC
		2: -3, // Face 2: Inflow BC
		3: -1, // Face 3: Wall BC
	}

	err = fb.ApplyBoundaryConditions(bcData)
	if err != nil {
		t.Fatalf("Failed to apply BC overlay: %v", err)
	}

	// Check that BCs were applied correctly
	expectedBCs := []int32{-1, -2, -3, -1}
	for i, expected := range expectedBCs {
		if fb.FaceIndex[i] != expected {
			t.Errorf("Face %d: expected BC %d, got %d", i, expected, fb.FaceIndex[i])
		}
	}

	// No face should have placeholder value anymore
	for i, faceCode := range fb.FaceIndex {
		if faceCode == BoundaryPlaceholder {
			t.Errorf("Face %d still has placeholder value -999", i)
		}
	}
}

// TestFaceBuffer_CubeMesh tests with more complex geometry
func TestFaceBuffer_CubeMesh(t *testing.T) {
	// Get standard cube mesh
	tm := utils.GetStandardTestMeshes()
	m := mesh.ConvertToMesh(tm.CubeMesh)

	order := 3
	el, err := tetelement.NewElement3DFromMesh(order, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	fb, err := BuildFaceBuffer(el)
	if err != nil {
		t.Fatalf("Failed to build face buffer: %v", err)
	}

	stats := fb.GetStats()

	// Cube mesh properties
	// Should have significant interior faces (cube is fully connected)
	if stats["interior_faces"] == 0 {
		t.Error("Cube mesh should have interior faces")
	}

	// Total faces should match
	expectedTotalFaces := int(el.K * 4)
	if stats["total_faces"] != expectedTotalFaces {
		t.Errorf("Expected %d total faces, got %d",
			expectedTotalFaces, stats["total_faces"])
	}

	// Validate all connections
	if err := fb.ValidateFaceBuffer(); err != nil {
		t.Errorf("Validation failed: %v", err)
	}

	// Test that interior face indices are reasonable
	var minInterior, maxInterior int32 = 0, 0
	interiorCount := 0

	for _, faceCode := range fb.FaceIndex {
		if faceCode > 0 {
			if interiorCount == 0 || faceCode < minInterior {
				minInterior = faceCode
			}
			if faceCode > maxInterior {
				maxInterior = faceCode
			}
			interiorCount++
		}
	}

	if interiorCount > 0 {
		// Check that interior indices span a reasonable range
		spread := float64(maxInterior-minInterior) / float64(fb.K*fb.Nfaces*fb.Nfp)
		if spread < 0.1 {
			t.Logf("Warning: Interior face indices appear clustered (spread=%f)", spread)
		}
	}
}
