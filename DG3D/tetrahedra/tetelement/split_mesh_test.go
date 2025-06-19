package tetelement

import (
	"fmt"
	mesh2 "github.com/notargets/gocfd/DG3D/mesh"
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"
)

// TestEtoPNormalization tests conversion to 0-based indexing
func TestEtoPNormalization(t *testing.T) {
	testCases := []struct {
		name     string
		input    []int
		expected []int
	}{
		{
			name:     "AlreadyZeroBased",
			input:    []int{0, 1, 2, 0, 1},
			expected: []int{0, 1, 2, 0, 1},
		},
		{
			name:     "OneBased",
			input:    []int{1, 2, 3, 1, 2},
			expected: []int{0, 1, 2, 0, 1},
		},
		{
			name:     "ArbitraryIDs",
			input:    []int{5, 7, 5, 9, 7},
			expected: []int{0, 1, 0, 2, 1},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			result := normalizeToZeroBased(tc.input)
			for i, v := range result {
				if v != tc.expected[i] {
					t.Errorf("Index %d: got %d, expected %d", i, v, tc.expected[i])
				}
			}
		})
	}
}

// TestBCMapping tests the primary functionality: BC transformation from global to local
func TestBCMapping(t *testing.T) {
	tm := utils.GetStandardTestMeshes()
	order := 2 // Polynomial order for DG3D

	testCases := []struct {
		name   string
		mesh   utils.CompleteMesh
		eToP   []int
		testBC map[string][]mesh2.BoundaryElement // Global BCs to test
	}{
		{
			name: "SingleTetSinglePartition",
			mesh: utils.CompleteMesh{
				Nodes:    tm.TetraNodes,
				Elements: []utils.ElementSet{tm.SingleTet},
			},
			eToP: []int{0}, // Single partition
			testBC: map[string][]mesh2.BoundaryElement{
				"wall": {
					{ElementType: utils.Triangle, ParentElement: 0, ParentFace: 0},
					{ElementType: utils.Triangle, ParentElement: 0, ParentFace: 1},
					{ElementType: utils.Triangle, ParentElement: 0, ParentFace: 2},
					{ElementType: utils.Triangle, ParentElement: 0, ParentFace: 3},
				},
			},
		},
		{
			name: "TwoTetsTwoPartitions",
			mesh: tm.TwoTetMesh,
			eToP: []int{0, 1}, // One element per partition
			testBC: map[string][]mesh2.BoundaryElement{
				"inlet": {
					{ElementType: utils.Triangle, ParentElement: 0, ParentFace: 0},
					{ElementType: utils.Triangle, ParentElement: 0, ParentFace: 1},
				},
				"outlet": {
					{ElementType: utils.Triangle, ParentElement: 1, ParentFace: 2},
					{ElementType: utils.Triangle, ParentElement: 1, ParentFace: 3},
				},
			},
		},
		{
			name: "CubeMeshThreePartitions",
			mesh: tm.CubeMesh,
			eToP: []int{0, 0, 1, 1, 2, 2}, // 6 tets in 3 partitions
			testBC: map[string][]mesh2.BoundaryElement{
				"bottom": {
					{ElementType: utils.Triangle, ParentElement: 0, ParentFace: 0},
					{ElementType: utils.Triangle, ParentElement: 2, ParentFace: 1},
				},
				"top": {
					{ElementType: utils.Triangle, ParentElement: 3, ParentFace: 2},
					{ElementType: utils.Triangle, ParentElement: 5, ParentFace: 3},
				},
				"side": {
					{ElementType: utils.Triangle, ParentElement: 1, ParentFace: 1},
					{ElementType: utils.Triangle, ParentElement: 4, ParentFace: 0},
				},
			},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Convert to mesh format
			globalMesh := mesh2.ConvertToMesh(tc.mesh)
			globalMesh.BoundaryElements = tc.testBC

			// Get vertex coordinates
			vx := make([]float64, globalMesh.NumVertices)
			vy := make([]float64, globalMesh.NumVertices)
			vz := make([]float64, globalMesh.NumVertices)
			for i, v := range globalMesh.Vertices {
				vx[i] = v[0]
				vy[i] = v[1]
				vz[i] = v[2]
			}

			// Create splitter and split
			globalMesh.EToP = tc.eToP
			el2, err := NewElement3DFromMesh(order, globalMesh)
			partElements, peToE, err := el2.Split, el2.PEToE, err
			if err != nil {
				t.Fatalf("SplitMesh failed: %v", err)
			}

			// Verify element structure
			for partID, el := range partElements {
				// Check Mesh is nil
				if el.Mesh != nil {
					t.Errorf("Partition %d: Mesh should be nil", partID)
				}

				// Check DG3D is initialized
				if el.DG3D == nil {
					t.Errorf("Partition %d: DG3D not initialized", partID)
				}

				// Check K matches number of elements
				expectedK := len(peToE[partID])
				if el.K != expectedK {
					t.Errorf("Partition %d: K=%d, expected %d", partID, el.K, expectedK)
				}

				// Check BCMaps is initialized
				if el.BCMaps == nil {
					t.Errorf("Partition %d: BCMaps not initialized", partID)
				}
			}

			// Verify ALL global BCs appear correctly in partitions
			for bcTag, globalBCs := range tc.testBC {
				t.Run(bcTag, func(t *testing.T) {
					foundCount := 0
					bcType := utils.ParseBCName(bcTag)

					for _, globalBC := range globalBCs {
						globalElemID := globalBC.ParentElement

						// Which partition should have this BC?
						partID := normalizeToZeroBased(tc.eToP)[globalElemID]
						el := partElements[partID]

						// Find local element ID
						localElemID := -1
						for i, gID := range peToE[partID] {
							if gID == globalElemID {
								localElemID = i
								break
							}
						}

						if localElemID == -1 {
							t.Errorf("Global element %d not found in partition %d",
								globalElemID, partID)
							continue
						}

						// Verify BC exists in BCMaps.FaceBC
						faceKey := FaceKey{
							Element: localElemID,
							Face:    globalBC.ParentFace,
						}

						if storedBCType, exists := el.BCMaps.FaceBC[faceKey]; exists {
							if storedBCType == bcType {
								foundCount++
							} else {
								t.Errorf("BC %s [elem:%d,face:%d] has wrong type: got %v, expected %v",
									bcTag, globalElemID, globalBC.ParentFace, storedBCType, bcType)
							}
						} else {
							t.Errorf("BC %s [elem:%d,face:%d] not found in partition %d FaceBC map",
								bcTag, globalElemID, globalBC.ParentFace, partID)
						}
					}

					if foundCount != len(globalBCs) {
						t.Errorf("BC %s: found %d/%d BCs", bcTag, foundCount, len(globalBCs))
					}
				})
			}

			// Verify no extra BCs were created
			for partID, el := range partElements {
				for faceKey, bcType := range el.BCMaps.FaceBC {
					localElemID := faceKey.Element
					faceID := faceKey.Face

					// Convert to global
					globalElemID := peToE[partID][localElemID]

					// This BC should exist in global mesh
					found := false
					for bcTag, globalBCs := range globalMesh.BoundaryElements {
						expectedBCType := utils.ParseBCName(bcTag)
						for _, globalBC := range globalBCs {
							if globalBC.ParentElement == globalElemID &&
								globalBC.ParentFace == faceID &&
								bcType == expectedBCType {
								found = true
								break
							}
						}
						if found {
							break
						}
					}

					if !found {
						t.Errorf("Partition %d has spurious BC type %v [local:%d,face:%d] = [global:%d,face:%d]",
							partID, bcType, localElemID, faceID, globalElemID, faceID)
					}
				}
			}
		})
	}
}

// TestVolumeConservation verifies total volume is preserved across partitions
func TestVolumeConservation(t *testing.T) {
	tm := utils.GetStandardTestMeshes()
	order := 1 // Low order is fine for volume tests

	testCases := []struct {
		name string
		mesh utils.CompleteMesh
		eToP []int
	}{
		{
			name: "CubeMeshTwoPartitions",
			mesh: tm.CubeMesh,
			eToP: []int{0, 0, 0, 1, 1, 1},
		},
		{
			name: "CubeMeshThreePartitions",
			mesh: tm.CubeMesh,
			eToP: []int{0, 0, 1, 1, 2, 2},
		},
		{
			name: "CubeMeshUneven",
			mesh: tm.CubeMesh,
			eToP: []int{0, 1, 1, 1, 2, 2},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Import mesh package functionality
			globalMesh := mesh2.ConvertToMesh(tc.mesh)

			// Extract vertex coordinates
			vx := make([]float64, globalMesh.NumVertices)
			vy := make([]float64, globalMesh.NumVertices)
			vz := make([]float64, globalMesh.NumVertices)
			for i, v := range globalMesh.Vertices {
				vx[i] = v[0]
				vy[i] = v[1]
				vz[i] = v[2]
			}

			// Create global Element3D to compute volume consistently
			globalEl, err := NewElement3DFromMesh(order, globalMesh)
			if err != nil {
				t.Fatalf("Failed to create global Element3D: %v", err)
			}

			// Compute global volume using average Jacobian
			// Without quadrature weights, use simple average
			globalVolume := 0.0
			for k := 0; k < globalEl.K; k++ {
				elemVolume := 0.0
				// Sum Jacobian values for this element
				for i := 0; i < globalEl.Np; i++ {
					idx := k*globalEl.Np + i
					elemVolume += globalEl.DG3D.J.DataP[idx]
				}
				// Average Jacobian times reference element volume (2/3 for tet)
				globalVolume += elemVolume * 2.0 / (3.0 * float64(globalEl.Np))
			}

			// Split mesh
			globalMesh.EToP = tc.eToP
			el2, err := NewElement3DFromMesh(order, globalMesh)
			partElements, _, err := el2.Split, el2.PEToE, err
			if err != nil {
				t.Fatalf("SplitMesh failed: %v", err)
			}

			// Compute partition volumes using same method
			partitionVolume := 0.0
			for _, el := range partElements {
				for k := 0; k < el.K; k++ {
					elemVolume := 0.0
					for i := 0; i < el.Np; i++ {
						idx := k*el.Np + i
						elemVolume += el.DG3D.J.DataP[idx]
					}
					partitionVolume += elemVolume * 2.0 / (3.0 * float64(el.Np))
				}
			}

			// Check conservation
			tolerance := 1e-12
			if math.Abs(globalVolume-partitionVolume) > tolerance {
				t.Errorf("Volume not conserved: global=%g, partitions=%g, diff=%g",
					globalVolume, partitionVolume, globalVolume-partitionVolume)
			}
		})
	}
}

// TestElementDistribution verifies each element appears in exactly one partition
func TestElementDistribution(t *testing.T) {
	tm := utils.GetStandardTestMeshes()
	globalMesh := mesh2.ConvertToMesh(tm.CubeMesh)
	order := 1

	// Extract actual vertex coordinates from the mesh
	vx := make([]float64, globalMesh.NumVertices)
	vy := make([]float64, globalMesh.NumVertices)
	vz := make([]float64, globalMesh.NumVertices)
	for i, v := range globalMesh.Vertices {
		vx[i] = v[0]
		vy[i] = v[1]
		vz[i] = v[2]
	}

	eToP := []int{0, 1, 0, 2, 1, 2}

	globalMesh.EToP = eToP
	el2, err := NewElement3DFromMesh(order, globalMesh)
	_, peToE, err := el2.Split, el2.PEToE, err
	if err != nil {
		t.Fatalf("SplitMesh failed: %v", err)
	}

	// Check each element appears exactly once
	seen := make(map[int]int) // globalElemID -> partID
	for partID, elems := range peToE {
		for _, globalElemID := range elems {
			if prevPartID, exists := seen[globalElemID]; exists {
				t.Errorf("Element %d appears in both partition %d and %d",
					globalElemID, prevPartID, partID)
			}
			seen[globalElemID] = partID
		}
	}

	// Check all elements are accounted for
	if len(seen) != globalMesh.NumElements {
		t.Errorf("Element count mismatch: seen %d, expected %d",
			len(seen), globalMesh.NumElements)
	}
}

// TestIncrementalPartitions tests systematic increase in partition count
// Following Unit Testing Principle: Test every intermediate level
func TestIncrementalPartitions(t *testing.T) {
	tm := utils.GetStandardTestMeshes()
	globalMesh := mesh2.ConvertToMesh(tm.CubeMesh)
	order := 2

	vx := make([]float64, globalMesh.NumVertices)
	vy := make([]float64, globalMesh.NumVertices)
	vz := make([]float64, globalMesh.NumVertices)
	for i, v := range globalMesh.Vertices {
		vx[i] = v[0]
		vy[i] = v[1]
		vz[i] = v[2]
	}

	// Test 1 through 6 partitions (6 elements in cube mesh)
	for numParts := 1; numParts <= 6; numParts++ {
		t.Run(fmt.Sprintf("%dPartitions", numParts), func(t *testing.T) {
			// Create partition assignment
			eToP := make([]int, 6)
			for i := 0; i < 6; i++ {
				eToP[i] = i % numParts
			}

			globalMesh.EToP = eToP
			el2, err := NewElement3DFromMesh(order, globalMesh)
			partElements, _, err := el2.Split, el2.PEToE, err
			if err != nil {
				t.Fatalf("SplitMesh failed: %v", err)
			}

			// Verify correct number of partitions
			if len(partElements) != numParts {
				t.Errorf("Expected %d partitions, got %d", numParts, len(partElements))
			}

			// Verify element distribution
			totalElems := 0
			for _, el := range partElements {
				totalElems += el.K
			}
			if totalElems != 6 {
				t.Errorf("Total elements %d != 6", totalElems)
			}

			// Verify each Element3D is properly initialized
			for partID, el := range partElements {
				if el.DG3D == nil {
					t.Errorf("Partition %d: DG3D not initialized", partID)
				}
				if el.Mesh != nil {
					t.Errorf("Partition %d: Mesh should be nil", partID)
				}
				if el.BCMaps == nil {
					t.Errorf("Partition %d: BCMaps not initialized", partID)
				}
			}
		})
	}
}

// TestDegeneratePartitions tests edge cases and degenerate configurations
// Following Unit Testing Principle: Test degeneracies and special cases
func TestDegeneratePartitions(t *testing.T) {
	tm := utils.GetStandardTestMeshes()
	order := 1

	testCases := []struct {
		name string
		mesh utils.CompleteMesh
		eToP []int
		desc string
	}{
		{
			name: "AllElementsInOnePartition",
			mesh: tm.CubeMesh,
			eToP: []int{0, 0, 0, 0, 0, 0}, // All elements in single partition
			desc: "All elements in first partition",
		},
		{
			name: "SingleElementPerPartition",
			mesh: tm.CubeMesh,
			eToP: []int{0, 1, 2, 3, 4, 5}, // Each element in its own partition
			desc: "Maximum partitioning",
		},
		{
			name: "UnbalancedPartitions",
			mesh: tm.CubeMesh,
			eToP: []int{0, 0, 0, 0, 0, 1}, // 5 elements in partition 0, 1 in partition 1
			desc: "Highly unbalanced",
		},
		{
			name: "NonSequentialPartitionIDs",
			mesh: tm.CubeMesh,
			eToP: []int{10, 20, 10, 30, 20, 30}, // Non-sequential IDs
			desc: "Tests EToP normalization",
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			globalMesh := mesh2.ConvertToMesh(tc.mesh)

			vx := make([]float64, globalMesh.NumVertices)
			vy := make([]float64, globalMesh.NumVertices)
			vz := make([]float64, globalMesh.NumVertices)
			for i, v := range globalMesh.Vertices {
				vx[i] = v[0]
				vy[i] = v[1]
				vz[i] = v[2]
			}

			globalMesh.EToP = tc.eToP
			el2, err := NewElement3DFromMesh(order, globalMesh)
			_, peToE, err := el2.Split, el2.PEToE, err
			if err != nil {
				t.Fatalf("SplitMesh failed for %s: %v", tc.desc, err)
			}

			// Verify all elements accounted for
			totalElems := 0
			for _, elems := range peToE {
				totalElems += len(elems)
			}
			if totalElems != globalMesh.NumElements {
				t.Errorf("%s: Lost elements, got %d expected %d",
					tc.name, totalElems, globalMesh.NumElements)
			}

			// Verify no duplicate elements
			seen := make(map[int]bool)
			for partID, elems := range peToE {
				for _, elem := range elems {
					if seen[elem] {
						t.Errorf("%s: Element %d duplicated in partition %d",
							tc.name, elem, partID)
					}
					seen[elem] = true
				}
			}
		})
	}
}

// Helper function
func tetVolume(v0, v1, v2, v3 []float64) float64 {
	// Volume = |det(v1-v0, v2-v0, v3-v0)| / 6
	d1 := []float64{v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]}
	d2 := []float64{v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]}
	d3 := []float64{v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]}

	det := d1[0]*(d2[1]*d3[2]-d2[2]*d3[1]) -
		d1[1]*(d2[0]*d3[2]-d2[2]*d3[0]) +
		d1[2]*(d2[0]*d3[1]-d2[1]*d3[0])

	return math.Abs(det) / 6.0
}
