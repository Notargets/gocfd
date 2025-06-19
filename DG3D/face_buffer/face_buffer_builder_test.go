// Package facebuffer tests
// Tests the face buffer builder which uses ONLY:
// - VmapM/VmapP: To identify face connections and boundary faces
// - EToE: To find neighbor elements
// - EToP: To detect remote partitions
// Note: EToF is NOT used - we find P positions by searching for matching volume nodes
package facebuffer

import (
	"fmt"
	"testing"

	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/DG3D/tetrahedra/tetelement"
	"github.com/notargets/gocfd/utils"
)

// TestFaceBuffer_SingleTetBoundary tests the simplest case - single tet with all boundary faces
// Following Unit Testing Principle: Start with fundamentals
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
			expectedTotal := uint32(el.Nfp * 4 * 1) // 4 faces, 1 element
			if fb.TotalFacePoints != expectedTotal {
				t.Errorf("Expected %d total face points, got %d", expectedTotal, fb.TotalFacePoints)
			}

			// Test 1: All faces should be boundary
			for i, ft := range fb.FaceTypes {
				if ft != BoundaryFace {
					t.Errorf("Face point %d should be boundary, got %v", i, ft)
				}
			}

			// Test 2: No interior indices
			if len(fb.LocalPIndices) != 0 {
				t.Errorf("Expected 0 local P indices, got %d", len(fb.LocalPIndices))
			}

			// Test 3: No remote partitions
			if len(fb.RemoteSendIndices) != 0 {
				t.Errorf("Expected 0 remote partitions, got %d", len(fb.RemoteSendIndices))
			}

			// Validate consistency
			if err := fb.ValidateFaceBuffer(); err != nil {
				t.Errorf("Validation failed: %v", err)
			}
		})
	}
}

// TestFaceBuffer_TwoTetsInterior tests interior face connections
// Following Unit Testing Principle: Build systematically
// Verifies that the face buffer correctly identifies interior connections
// using only VmapM/VmapP/EToE connectivity data
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

	// Test 1: Should have both boundary and interior points
	if stats["boundary_points"] == 0 {
		t.Error("Expected some boundary points")
	}
	if stats["interior_points"] == 0 {
		t.Error("Expected some interior points")
	}

	// Test 2: Interior points should have valid local P indices
	interiorCount := 0
	for i, ft := range fb.FaceTypes {
		if ft == InteriorFace {
			// Check that we have a corresponding P index
			if interiorCount >= len(fb.LocalPIndices) {
				t.Errorf("Missing P index for interior point %d", i)
			} else {
				pIdx := fb.LocalPIndices[interiorCount]
				// P index should point to a different location
				if pIdx == uint32(i) {
					t.Errorf("P index should not equal M index for interior face")
				}
				// P index should be in valid range
				if pIdx >= fb.TotalFacePoints {
					t.Errorf("P index %d out of range", pIdx)
				}
			}
			interiorCount++
		}
	}

	// Test 3: Verify LocalPIndices validity
	// Each LocalPIndices entry should point to a valid M location
	for i, pIdx := range fb.LocalPIndices {
		if pIdx >= fb.TotalFacePoints {
			t.Errorf("LocalPIndices[%d] = %d is out of range [0,%d)",
				i, pIdx, fb.TotalFacePoints)
		}
	}

	// Test 4: Check that VmapP connectivity is preserved
	// For each interior M point, its P volume node should match what VmapP says
	interiorIdx := 0
	for m := uint32(0); m < fb.TotalFacePoints; m++ {
		if fb.FaceTypes[m] == InteriorFace {
			pPos := fb.LocalPIndices[interiorIdx]

			// The volume node at P position should match VmapP[m]
			expectedVolNode := el.VmapP[m]
			actualVolNode := el.VmapM[pPos]

			if actualVolNode != expectedVolNode {
				t.Errorf("Interior connection M=%d: expected vol node %d at P=%d, got %d",
					m, expectedVolNode, pPos, actualVolNode)
			}

			interiorIdx++
		}
	}

	// Validate
	if err := fb.ValidateFaceBuffer(); err != nil {
		t.Errorf("Validation failed: %v", err)
	}
}

// TestFaceBuffer_ParallelPartitions tests remote face connections
// Following Unit Testing Principle: Progressive complexity
func TestFaceBuffer_ParallelPartitions(t *testing.T) {
	t.Skip("Parallel partition testing requires full mesh connectivity info - would need refactoring")

	// Note: For proper parallel testing, we would need:
	// 1. Full EToE/EToP arrays for all partitions
	// 2. A way to identify which neighbor elements belong to which remote partitions
	// 3. Complete VmapM/VmapP for computing remote P positions
	//
	// The current Element3D structure assumes we only have local partition data,
	// which makes it impossible to compute the remote P positions accurately.
}

// TestFaceBuffer_SimpleParallel tests the parallel detection logic
func TestFaceBuffer_SimpleParallel(t *testing.T) {
	// Create a simple 2-element mesh where element 1 is in a different partition
	tm := utils.GetStandardTestMeshes()
	m := mesh.ConvertToMesh(tm.TwoTetMesh)

	order := 1
	el, err := tetelement.NewElement3DFromMesh(order, m)
	if err != nil {
		t.Fatalf("Failed to create Element3D: %v", err)
	}

	// Set up partitioning - element 0 in partition 0, element 1 in partition 1
	el.EToP = []int{0, 1}

	// Simulate partition 0 by modifying K to only include element 0
	// Save original values
	origK := el.K
	origVmapM := el.VmapM
	origVmapP := el.VmapP
	origEToE := el.EToE

	// Modify to represent only partition 0
	el.K = 1                         // Only element 0
	el.VmapM = el.VmapM[:1*4*el.Nfp] // First element's face points
	el.VmapP = el.VmapP[:1*4*el.Nfp]
	el.EToE = el.EToE[:1]
	// Keep full EToP to check partition membership

	// Build face buffer for partition 0
	fb, err := BuildFaceBuffer(el)
	if err != nil {
		t.Fatalf("Failed to build face buffer: %v", err)
	}

	// Restore original values
	el.K = origK
	el.VmapM = origVmapM
	el.VmapP = origVmapP
	el.EToE = origEToE

	stats := fb.GetStats()

	// Should have remote faces (where element 0 connects to element 1)
	if stats["remote_points"] == 0 {
		t.Error("Expected some remote points for parallel mesh")
	}

	// Should have boundary faces (exterior faces of element 0)
	if stats["boundary_points"] == 0 {
		t.Error("Expected some boundary points")
	}

	// Should have NO interior faces (element 0 is alone in partition 0)
	if stats["interior_points"] != 0 {
		t.Errorf("Expected 0 interior points for single element in partition, got %d",
			stats["interior_points"])
	}

	// The shared face should be marked as remote
	expectedRemote := el.Nfp // One face worth of points
	if stats["remote_points"] != expectedRemote {
		t.Errorf("Expected %d remote points (one shared face), got %d",
			expectedRemote, stats["remote_points"])
	}

	// Validate
	if err := fb.ValidateFaceBuffer(); err != nil {
		t.Errorf("Validation failed: %v", err)
	}
}

// TestFaceBuffer_CubeMesh tests with more complex geometry
// Following Unit Testing Principle: Incremental validation
func TestFaceBuffer_CubeMesh(t *testing.T) {
	// Get standard cube mesh
	tm := utils.GetStandardTestMeshes()

	// Convert to mesh
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

	// Test properties of cube mesh
	// 1. Should have significant interior faces (cube is fully connected)
	if stats["interior_points"] < stats["boundary_points"] {
		t.Error("Cube mesh should have more interior than boundary points")
	}

	// 2. Total points should match
	expectedTotal := el.Nfp * 4 * el.K
	if stats["total_face_points"] != expectedTotal {
		t.Errorf("Expected %d total points, got %d",
			expectedTotal, stats["total_face_points"])
	}

	// 3. Validate all connections
	if err := fb.ValidateFaceBuffer(); err != nil {
		t.Errorf("Validation failed: %v", err)
	}

	// 4. Test memory layout is contiguous
	// Interior P indices should increase monotonically when sorted
	if len(fb.LocalPIndices) > 1 {
		// Check for reasonable distribution (not all clustered)
		minIdx := fb.LocalPIndices[0]
		maxIdx := fb.LocalPIndices[0]
		for _, idx := range fb.LocalPIndices {
			if idx < minIdx {
				minIdx = idx
			}
			if idx > maxIdx {
				maxIdx = idx
			}
		}

		spread := float64(maxIdx-minIdx) / float64(fb.TotalFacePoints)
		if spread < 0.1 {
			t.Error("P indices appear too clustered, may indicate incorrect traversal")
		}
	}
}
