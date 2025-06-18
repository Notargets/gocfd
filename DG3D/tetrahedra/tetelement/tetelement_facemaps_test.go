package tetelement

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/utils"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

// Level 0: Single tetrahedron - simplest possible case
func TestBuildMaps3D_SingleTet(t *testing.T) {
	// Create a single tetrahedron using test helpers
	tm := utils.GetStandardTestMeshes()
	singleTetMesh := &utils.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []utils.ElementSet{tm.SingleTet},
		Dimension: 3,
		BoundingBox: [2][3]float64{
			{0, 0, 0},
			{1, 1, 1},
		},
	}

	// Convert to mesh
	meshObj := mesh.ConvertToMesh(*singleTetMesh)

	// Test multiple polynomial orders
	orders := []int{1, 2, 3}

	for _, order := range orders {
		t.Run(fmt.Sprintf("Order_%d", order), func(t *testing.T) {
			// Create Element3D
			el, err := NewElement3DFromMesh(order, meshObj)
			require.NoError(t, err)
			require.NotNil(t, el)

			// Basic validations
			assert.Equal(t, 1, el.K, "Should have 1 element")
			assert.NotNil(t, el.VmapM, "VmapM should be initialized")
			assert.NotNil(t, el.VmapP, "VmapP should be initialized")

			// Expected sizes
			Nfp := el.Nfp
			Nfaces := 4
			expectedSize := Nfp * Nfaces * el.K

			assert.Equal(t, expectedSize, len(el.VmapM), "VmapM size")
			assert.Equal(t, expectedSize, len(el.VmapP), "VmapP size")

			// For a single tet, all faces are boundary
			// Therefore VmapP should equal VmapM everywhere
			boundaryCount := 0
			for i := 0; i < len(el.VmapM); i++ {
				if el.VmapP[i] == el.VmapM[i] {
					boundaryCount++
				} else {
					t.Errorf("Single tet: VmapP[%d]=%d should equal VmapM[%d]=%d",
						i, el.VmapP[i], i, el.VmapM[i])
				}
			}

			assert.Equal(t, expectedSize, boundaryCount,
				"All face points should be boundary")

			// Verify VmapM maps to valid volume nodes
			for i := 0; i < len(el.VmapM); i++ {
				vmapM := el.VmapM[i]
				elem := vmapM / el.Np
				node := vmapM % el.Np

				assert.Equal(t, 0, elem, "VmapM should map to element 0")
				assert.Less(t, node, el.Np, "Node index should be < Np")

				// Decode which face point this is
				face := (i / Nfp) % Nfaces
				point := i % Nfp

				// VmapM should map face nodes to volume nodes
				// The volume node should be one of the face nodes
				faceNodeIdx := el.Fmask[face][point]
				expectedVmapM := 0*el.Np + faceNodeIdx

				assert.Equal(t, expectedVmapM, vmapM,
					"VmapM[%d] (face %d, point %d) should map to volume node %d",
					i, face, point, faceNodeIdx)
			}

			// Verify MapB (boundary map) contains all face points
			assert.Equal(t, expectedSize, len(el.MapB),
				"All face points should be in boundary map")
		})
	}
}

// Level 1: Two tetrahedra sharing a face
func TestBuildMaps3D_TwoTetsSharedFace(t *testing.T) {
	// Use the standard two-tet mesh
	tm := utils.GetStandardTestMeshes()
	meshObj := mesh.ConvertToMesh(tm.TwoTetMesh)

	// Test with order 1 and 2
	orders := []int{1, 2}

	for _, order := range orders {
		t.Run(fmt.Sprintf("Order_%d", order), func(t *testing.T) {
			el, err := NewElement3DFromMesh(order, meshObj)
			require.NoError(t, err)

			// Basic validations
			assert.Equal(t, 2, el.K, "Should have 2 elements")
			assert.Equal(t, 5, el.Mesh.NumVertices, "Should have 5 vertices")

			// Find the shared face using EToE/EToF
			sharedFaces := make(map[string]bool)
			for elem := 0; elem < 2; elem++ {
				for face := 0; face < 4; face++ {
					neighbor := el.EToE[elem][face]
					if neighbor != elem && neighbor >= 0 {
						key := fmt.Sprintf("elem%d_face%d", elem, face)
						sharedFaces[key] = true
					}
				}
			}

			assert.Equal(t, 2, len(sharedFaces),
				"Should have exactly 2 shared face entries (one from each side)")

			// Count boundary vs interior face points
			boundaryCount := 0
			interiorCount := 0

			for i := 0; i < len(el.VmapM); i++ {
				if el.VmapP[i] == el.VmapM[i] {
					boundaryCount++
				} else {
					interiorCount++
				}
			}

			expectedBoundaryPoints := 6 * el.Nfp // 6 boundary faces
			expectedInteriorPoints := 2 * el.Nfp // 2 sides of shared face

			assert.Equal(t, expectedBoundaryPoints, boundaryCount,
				"Should have correct number of boundary points")
			assert.Equal(t, expectedInteriorPoints, interiorCount,
				"Should have correct number of interior points")

			// Verify reciprocal mappings for interior faces
			t.Run("ReciprocalMappings", func(t *testing.T) {
				verifyReciprocalMappings(t, el)
			})

			// Verify node coordinate matching
			t.Run("NodeCoordinateMatching", func(t *testing.T) {
				verifyNodeCoordinateMatching(t, el)
			})
		})
	}
}

// Level 2: Cube mesh with 6 tetrahedra
func TestBuildMaps3D_CubeMesh(t *testing.T) {
	tm := utils.GetStandardTestMeshes()
	meshObj := mesh.ConvertToMesh(tm.CubeMesh)

	// Test with order 1
	el, err := NewElement3DFromMesh(1, meshObj)
	require.NoError(t, err)

	assert.Equal(t, 6, el.K, "Cube should have 6 tetrahedra")

	// Verify connectivity is reasonable
	// Each tet in the cube mesh should have some interior faces
	interiorFaceCount := 0
	boundaryFaceCount := 0

	for elem := 0; elem < el.K; elem++ {
		for face := 0; face < 4; face++ {
			if el.EToE[elem][face] == elem {
				boundaryFaceCount++
			} else {
				interiorFaceCount++
			}
		}
	}

	// Cube has 12 exterior faces (2 per cube face * 6 cube faces)
	// Total faces = 6 tets * 4 faces = 24
	// Interior faces are counted twice
	assert.Equal(t, 12, boundaryFaceCount, "Should have 12 boundary faces")
	assert.Equal(t, 12, interiorFaceCount, "Should have 12 interior face connections")

	// Verify face mappings
	boundaryPoints := 0
	interiorPoints := 0
	invalidPoints := 0

	for i := 0; i < len(el.VmapP); i++ {
		if el.VmapP[i] == el.VmapM[i] {
			boundaryPoints++
		} else {
			// Check if VmapP points to a valid element
			elemP := el.VmapP[i] / el.Np
			if elemP >= el.K {
				invalidPoints++
				elem := i / (el.Nfp * 4)
				face := (i / el.Nfp) % 4
				point := i % el.Nfp
				t.Errorf("Invalid VmapP[%d]=%d (elem %d, face %d, point %d)",
					i, el.VmapP[i], elem, face, point)
			} else {
				interiorPoints++
			}
		}
	}

	assert.Equal(t, 0, invalidPoints, "Should have no invalid mappings")
	assert.Greater(t, boundaryPoints, 0, "Should have boundary points")
	assert.Greater(t, interiorPoints, 0, "Should have interior points")

	totalPoints := el.Nfp * 4 * el.K
	assert.Equal(t, totalPoints, boundaryPoints+interiorPoints,
		"All points should be classified")

	// Verify reciprocal mappings
	t.Run("ReciprocalMappings", func(t *testing.T) {
		verifyReciprocalMappings(t, el)
	})
}

// Helper function to verify reciprocal mappings
func verifyReciprocalMappings(t *testing.T, el *Element3D) {
	// For interior faces, mappings should be reciprocal
	// If face point A maps to face point B, then B should map back to A

	reciprocalFailures := 0
	checked := make(map[int]bool)

	for i := 0; i < len(el.VmapP); i++ {
		// Skip if already checked or if boundary
		if checked[i] || el.VmapP[i] == el.VmapM[i] {
			continue
		}

		// This is an interior face point
		// Find what it maps to
		vmapM_i := el.VmapM[i]
		vmapP_i := el.VmapP[i]
		mapP_i := el.MapP[i]

		// The reciprocal point should map back
		if mapP_i >= 0 && mapP_i < len(el.VmapP) {
			vmapM_j := el.VmapM[mapP_i]
			vmapP_j := el.VmapP[mapP_i]
			mapP_j := el.MapP[mapP_i]

			// Check reciprocity
			if vmapP_i != vmapM_j {
				reciprocalFailures++
				if reciprocalFailures <= 5 { // Limit output
					t.Errorf("VmapP reciprocity failed at %d: "+
						"VmapP[%d]=%d but VmapM[%d]=%d",
						i, i, vmapP_i, mapP_i, vmapM_j)
				}
			}

			if vmapP_j != vmapM_i {
				reciprocalFailures++
				if reciprocalFailures <= 5 {
					t.Errorf("VmapP reciprocity failed at %d: "+
						"VmapP[%d]=%d but VmapM[%d]=%d",
						mapP_i, mapP_i, vmapP_j, i, vmapM_i)
				}
			}

			if mapP_j != i {
				reciprocalFailures++
				if reciprocalFailures <= 5 {
					t.Errorf("MapP reciprocity failed: "+
						"MapP[%d]=%d but MapP[%d]=%d (should be %d)",
						i, mapP_i, mapP_i, mapP_j, i)
				}
			}

			checked[i] = true
			checked[mapP_i] = true
		}
	}

	assert.Equal(t, 0, reciprocalFailures,
		"Should have no reciprocal mapping failures")
}

// Helper function to verify node coordinates match
func verifyNodeCoordinateMatching(t *testing.T, el *Element3D) {
	NODETOL := 1e-10
	mismatchCount := 0

	// For each face point
	for i := 0; i < len(el.VmapM); i++ {
		// Skip boundary faces
		if el.VmapM[i] == el.VmapP[i] {
			continue
		}

		// Get M and P node coordinates
		vmapM := el.VmapM[i]
		vmapP := el.VmapP[i]

		elemM := vmapM / el.Np
		nodeM := vmapM % el.Np
		elemP := vmapP / el.Np
		nodeP := vmapP % el.Np

		// Get coordinates
		xM := el.X.At(nodeM, elemM)
		yM := el.Y.At(nodeM, elemM)
		zM := el.Z.At(nodeM, elemM)

		xP := el.X.At(nodeP, elemP)
		yP := el.Y.At(nodeP, elemP)
		zP := el.Z.At(nodeP, elemP)

		// Check distance
		dx := xM - xP
		dy := yM - yP
		dz := zM - zP
		dist := math.Sqrt(dx*dx + dy*dy + dz*dz)

		if dist > NODETOL {
			mismatchCount++
			if mismatchCount <= 5 { // Limit output
				elem := i / (el.Nfp * 4)
				face := (i / el.Nfp) % 4
				point := i % el.Nfp

				t.Errorf("Node coordinate mismatch at face point %d "+
					"(elem=%d,face=%d,point=%d): distance=%g",
					i, elem, face, point, dist)
			}
		}
	}

	assert.Equal(t, 0, mismatchCount,
		"Interior face nodes should have matching coordinates")
}

// Test specific known face mappings for order 1 tet
func TestBuildMaps3D_SpecificMappings_Order1(t *testing.T) {
	// Create a two-tet mesh with known connectivity
	// We'll verify specific node mappings
	tm := utils.GetStandardTestMeshes()
	meshObj := mesh.ConvertToMesh(tm.TwoTetMesh)

	el, err := NewElement3DFromMesh(1, meshObj)
	require.NoError(t, err)

	// For order 1, Nfp = 3 (3 nodes per face)
	assert.Equal(t, 3, el.Nfp)

	// The two tets share face with vertices {1,2,3}
	// Need to find which face this is for each tet

	// Tet 0 has vertices {0,1,2,3}
	// Face vertex mappings:
	// Face 0: {0,1,2}
	// Face 1: {0,1,3}
	// Face 2: {1,2,3} <- This is the shared face
	// Face 3: {0,2,3}

	// Tet 1 has vertices {1,2,3,4}
	// Face 0: {1,2,3} <- This is the shared face
	// Face 1: {1,2,4}
	// Face 2: {2,3,4}
	// Face 3: {1,3,4}

	// So Tet 0 face 2 should connect to Tet 1 face 0
	assert.Equal(t, 1, el.EToE[0][2], "Tet 0 face 2 should connect to Tet 1")
	assert.Equal(t, 0, el.EToF[0][2], "Tet 0 face 2 should connect to Tet 1 face 0")

	assert.Equal(t, 0, el.EToE[1][0], "Tet 1 face 0 should connect to Tet 0")
	assert.Equal(t, 2, el.EToF[1][0], "Tet 1 face 0 should connect to Tet 0 face 2")

	// Verify face point mappings for the shared face
	// Face points for Tet 0, face 2
	for p := 0; p < el.Nfp; p++ {
		idx0 := 0*4*el.Nfp + 2*el.Nfp + p // elem 0, face 2, point p

		// This should map to some point on Tet 1, face 0
		vmapP := el.VmapP[idx0]
		elemP := vmapP / el.Np

		assert.Equal(t, 1, elemP,
			"Tet 0 face 2 point %d should map to Tet 1", p)

		// Verify the mapping is reciprocal through MapP
		mapP := el.MapP[idx0]
		elem1 := mapP / (el.Nfp * 4)
		face1 := (mapP / el.Nfp) % 4

		assert.Equal(t, 1, elem1,
			"MapP should point to Tet 1")
		assert.Equal(t, 0, face1,
			"MapP should point to Tet 1 face 0")
	}
}

// Test that demonstrates debugging a specific face mapping issue
func TestBuildMaps3D_DebuggingExample(t *testing.T) {
	// This test shows how to debug face mapping issues
	// using a simple, known configuration

	tm := utils.GetStandardTestMeshes()
	meshObj := mesh.ConvertToMesh(tm.TwoTetMesh)

	el, err := NewElement3DFromMesh(1, meshObj)
	require.NoError(t, err)

	// Focus on a specific face connection
	elem := 0
	face := 2 // The shared face

	t.Logf("=== Debugging Element %d, Face %d ===", elem, face)
	t.Logf("Connects to Element %d, Face %d",
		el.EToE[elem][face], el.EToF[elem][face])

	// Print all mappings for this face
	for p := 0; p < el.Nfp; p++ {
		idx := elem*4*el.Nfp + face*el.Nfp + p

		vmapM := el.VmapM[idx]
		vmapP := el.VmapP[idx]
		mapP := el.MapP[idx]

		// Decode the mappings
		elemM := vmapM / el.Np
		nodeM := vmapM % el.Np
		elemP := vmapP / el.Np
		nodeP := vmapP % el.Np

		// Get coordinates
		xM := el.X.At(nodeM, elemM)
		yM := el.Y.At(nodeM, elemM)
		zM := el.Z.At(nodeM, elemM)

		var xP, yP, zP float64
		if elemP < el.K {
			xP = el.X.At(nodeP, elemP)
			yP = el.Y.At(nodeP, elemP)
			zP = el.Z.At(nodeP, elemP)
		}

		t.Logf("Point %d:", p)
		t.Logf("  VmapM[%d] = %d (elem %d, node %d) @ (%.3f, %.3f, %.3f)",
			idx, vmapM, elemM, nodeM, xM, yM, zM)
		t.Logf("  VmapP[%d] = %d (elem %d, node %d) @ (%.3f, %.3f, %.3f)",
			idx, vmapP, elemP, nodeP, xP, yP, zP)
		t.Logf("  MapP[%d] = %d", idx, mapP)

		// Verify coordinates match
		if elemP < el.K {
			dist := math.Sqrt(math.Pow(xM-xP, 2) +
				math.Pow(yM-yP, 2) + math.Pow(zM-zP, 2))
			t.Logf("  Distance between M and P: %g", dist)
		}
	}
}
