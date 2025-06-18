package gonudg

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/utils"
)

// Helper function to convert CompleteMesh to DG3D format
func meshToDG3DFormat(m utils.CompleteMesh) (VX, VY, VZ []float64, EToV [][]int) {
	// Extract vertex coordinates
	numNodes := len(m.Nodes.Nodes)
	VX = make([]float64, numNodes)
	VY = make([]float64, numNodes)
	VZ = make([]float64, numNodes)

	for i := 0; i < numNodes; i++ {
		VX[i] = m.Nodes.Nodes[i][0]
		VY[i] = m.Nodes.Nodes[i][1]
		VZ[i] = m.Nodes.Nodes[i][2]
	}

	// Extract element connectivity
	// Assuming we only have tetrahedral elements
	for _, elemSet := range m.Elements {
		if elemSet.Type == utils.Tet {
			EToV = make([][]int, len(elemSet.Elements))
			for i, elem := range elemSet.Elements {
				EToV[i] = make([]int, 4)
				for j, nodeName := range elem {
					// Convert node name to index (0-based)
					EToV[i][j] = m.Nodes.NodeMap[nodeName]
				}
			}
			break
		}
	}

	return
}

// TestBuildMaps3DBasic tests basic properties of the connectivity maps
func TestBuildMaps3DBasic(t *testing.T) {
	// Use the standard single tet from test helpers
	tm := utils.GetStandardTestMeshes()

	// Create a simple single tet mesh
	singleTetMesh := utils.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []utils.ElementSet{tm.SingleTet},
		Dimension: 3,
	}

	VX, VY, VZ, EToV := meshToDG3DFormat(singleTetMesh)

	for _, N := range []int{1, 2, 3} {
		t.Run(fmt.Sprintf("N=%d", N), func(t *testing.T) {
			dg, err := NewDG3D(N, VX, VY, VZ, EToV)
			if err != nil {
				t.Fatalf("Failed to create DG3D: %v", err)
			}

			// Check array dimensions
			expectedSize := dg.Nfp * dg.Nfaces * dg.K

			if len(dg.vmapM) != expectedSize {
				t.Errorf("vmapM size: got %d, want %d", len(dg.vmapM), expectedSize)
			}
			if len(dg.vmapP) != expectedSize {
				t.Errorf("vmapP size: got %d, want %d", len(dg.vmapP), expectedSize)
			}
			if len(dg.mapM) != expectedSize {
				t.Errorf("mapM size: got %d, want %d", len(dg.mapM), expectedSize)
			}
			if len(dg.mapP) != expectedSize {
				t.Errorf("mapP size: got %d, want %d", len(dg.mapP), expectedSize)
			}

			// Check that mapM contains sequential indices
			for i := 0; i < expectedSize; i++ {
				if dg.mapM[i] != i {
					t.Errorf("mapM[%d] = %d, expected %d", i, dg.mapM[i], i)
				}
			}

			// For a single element, all faces are boundaries
			// So mapP should equal mapM
			for i := 0; i < expectedSize; i++ {
				if dg.mapP[i] != dg.mapM[i] {
					t.Errorf("mapP[%d] = %d, expected %d (boundary)", i, dg.mapP[i], dg.mapM[i])
				}
			}

			// vmapP should equal vmapM for boundaries
			for i := 0; i < expectedSize; i++ {
				if dg.vmapP[i] != dg.vmapM[i] {
					t.Errorf("vmapP[%d] = %d, expected %d (boundary)", i, dg.vmapP[i], dg.vmapM[i])
				}
			}

			// Check vmapM points to valid volume nodes
			for i := 0; i < expectedSize; i++ {
				vid := dg.vmapM[i]
				if vid < 0 || vid >= dg.Np*dg.K {
					t.Errorf("vmapM[%d] = %d is out of range [0, %d)", i, vid, dg.Np*dg.K)
				}
			}

			// All nodes should be boundary nodes
			if len(dg.mapB) != expectedSize {
				t.Errorf("mapB size: got %d, want %d", len(dg.mapB), expectedSize)
			}
			if len(dg.vmapB) != expectedSize {
				t.Errorf("vmapB size: got %d, want %d", len(dg.vmapB), expectedSize)
			}
		})
	}
}

// TestBuildMaps3DFaceMapping verifies face node to volume node mapping
func TestBuildMaps3DFaceMapping(t *testing.T) {
	// Use standard single tet
	tm := utils.GetStandardTestMeshes()
	singleTetMesh := utils.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []utils.ElementSet{tm.SingleTet},
		Dimension: 3,
	}

	VX, VY, VZ, EToV := meshToDG3DFormat(singleTetMesh)

	N := 2
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// Verify that face nodes map to correct volume nodes
	// For each face, check that the mapped volume nodes satisfy the face equation
	tolerance := 1e-14

	// Face 0: t = -1
	for i := 0; i < dg.Nfp; i++ {
		faceIdx := 0*dg.Nfp + i
		volIdx := dg.vmapM[faceIdx]

		// Get volume node in element-local index
		localIdx := volIdx % dg.Np
		tVal := dg.t[localIdx]

		if math.Abs(tVal+1.0) > tolerance {
			t.Errorf("Face 0, node %d: t = %f, expected -1", i, tVal)
		}
	}

	// Face 1: s = -1
	for i := 0; i < dg.Nfp; i++ {
		faceIdx := 1*dg.Nfp + i
		volIdx := dg.vmapM[faceIdx]
		localIdx := volIdx % dg.Np
		sVal := dg.s[localIdx]

		if math.Abs(sVal+1.0) > tolerance {
			t.Errorf("Face 1, node %d: s = %f, expected -1", i, sVal)
		}
	}

	// Face 2: r+s+t = -1
	for i := 0; i < dg.Nfp; i++ {
		faceIdx := 2*dg.Nfp + i
		volIdx := dg.vmapM[faceIdx]
		localIdx := volIdx % dg.Np
		sum := dg.r[localIdx] + dg.s[localIdx] + dg.t[localIdx]

		if math.Abs(sum+1.0) > tolerance {
			t.Errorf("Face 2, node %d: r+s+t = %f, expected -1", i, sum)
		}
	}

	// Face 3: r = -1
	for i := 0; i < dg.Nfp; i++ {
		faceIdx := 3*dg.Nfp + i
		volIdx := dg.vmapM[faceIdx]
		localIdx := volIdx % dg.Np
		rVal := dg.r[localIdx]

		if math.Abs(rVal+1.0) > tolerance {
			t.Errorf("Face 3, node %d: r = %f, expected -1", i, rVal)
		}
	}
}

// TestBuildMaps3DTwoElements tests connectivity between two tetrahedra
func TestBuildMaps3DTwoElements(t *testing.T) {
	// Use the standard TwoTetMesh
	tm := utils.GetStandardTestMeshes()
	VX, VY, VZ, EToV := meshToDG3DFormat(tm.TwoTetMesh)

	N := 2
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// We need to set up connectivity manually since tiConnect3D is not implemented
	// For the TwoTetMesh, we need to determine which faces are connected
	// Convert to mesh format to get connectivity
	meshObj := mesh.ConvertToMesh(tm.TwoTetMesh)
	meshObj.BuildConnectivity()

	// Copy connectivity from mesh
	dg.EToE = meshObj.EToE
	dg.EToF = meshObj.EToF

	// Rebuild maps with connectivity
	dg.BuildMaps3D()

	// Test that shared face nodes are properly connected
	NF := dg.Nfp * dg.Nfaces

	// Count connected nodes (nodes that are not on boundaries)
	connectedCount := 0
	for i := 0; i < dg.K*NF; i++ {
		if dg.mapP[i] != i {
			connectedCount++
		}
	}

	// The TwoTetMesh has two tets sharing one face
	// Each tet has 4 faces, so total 8 faces
	// They share 1 face, so we have 6 boundary faces and 2 interior faces
	// Interior faces contribute Nfp nodes each, so 2*Nfp interior nodes
	expectedConnected := 2 * dg.Nfp // One shared face seen from both sides
	if connectedCount != expectedConnected {
		t.Logf("Nfp = %d, Nfaces = %d, K = %d", dg.Nfp, dg.Nfaces, dg.K)
		t.Logf("Total face nodes = %d", dg.K*NF)

		// Debug: print connectivity info
		for k := 0; k < dg.K; k++ {
			for f := 0; f < dg.Nfaces; f++ {
				neighbor := dg.EToE[k][f]
				neighborFace := dg.EToF[k][f]
				if neighbor != k || neighborFace != f {
					t.Logf("Element %d face %d connects to element %d face %d",
						k, f, neighbor, neighborFace)
				}
			}
		}
		t.Errorf("Connected nodes: got %d, expected %d", connectedCount, expectedConnected)
	}

	// Test that boundary nodes are correctly identified
	// The TwoTetMesh should have 6 boundary faces (4 faces per tet - 1 shared face per tet = 3 per tet, times 2 tets)
	expectedBoundaryNodes := 6 * dg.Nfp
	if len(dg.mapB) != expectedBoundaryNodes {
		t.Logf("mapB length = %d", len(dg.mapB))
		t.Logf("Expected %d boundary faces * %d Nfp = %d boundary nodes",
			6, dg.Nfp, expectedBoundaryNodes)
		t.Errorf("Boundary nodes: got %d, expected %d", len(dg.mapB), expectedBoundaryNodes)
	}
}

// TestBuildMaps3DNodeMatching tests that shared face nodes are correctly matched
func TestBuildMaps3DNodeMatching(t *testing.T) {
	// Use the cube mesh which has more complex connectivity
	tm := utils.GetStandardTestMeshes()
	VX, VY, VZ, EToV := meshToDG3DFormat(tm.CubeMesh)

	N := 3
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// Get connectivity from mesh
	meshObj := mesh.ConvertToMesh(tm.CubeMesh)
	meshObj.BuildConnectivity()
	dg.EToE = meshObj.EToE
	dg.EToF = meshObj.EToF

	dg.BuildMaps3D()

	// For connected nodes, verify that physical coordinates match
	NF := dg.Nfp * dg.Nfaces
	tolerance := math.Sqrt(dg.NODETOL)

	for k := 0; k < dg.K; k++ {
		for f := 0; f < dg.Nfaces; f++ {
			for i := 0; i < dg.Nfp; i++ {
				idx := k*NF + f*dg.Nfp + i

				if dg.mapP[idx] != idx {
					// This is a connected node
					vidM := dg.vmapM[idx]
					vidP := dg.vmapP[idx]

					// Extract physical coordinates
					// Global index = element * Np + local index
					kM := vidM / dg.Np
					iM := vidM % dg.Np
					kP := vidP / dg.Np
					iP := vidP % dg.Np

					xM := dg.x.At(iM, kM)
					yM := dg.y.At(iM, kM)
					zM := dg.z.At(iM, kM)

					xP := dg.x.At(iP, kP)
					yP := dg.y.At(iP, kP)
					zP := dg.z.At(iP, kP)

					// Check that coordinates match
					dist := math.Sqrt((xP-xM)*(xP-xM) + (yP-yM)*(yP-yM) + (zP-zM)*(zP-zM))
					if dist > tolerance {
						t.Errorf("Connected nodes don't match: distance = %e", dist)
						t.Errorf("  Node M: (%f, %f, %f)", xM, yM, zM)
						t.Errorf("  Node P: (%f, %f, %f)", xP, yP, zP)
					}
				}
			}
		}
	}
}

// TestBuildMaps3DBoundaryFaces tests the GetBoundaryFaces helper function
func TestBuildMaps3DBoundaryFaces(t *testing.T) {
	// Single tetrahedron - all faces are boundaries
	tm := utils.GetStandardTestMeshes()
	singleTetMesh := utils.CompleteMesh{
		Nodes:     tm.TetraNodes,
		Elements:  []utils.ElementSet{tm.SingleTet},
		Dimension: 3,
	}

	VX, VY, VZ, EToV := meshToDG3DFormat(singleTetMesh)

	N := 1
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// Need to initialize EToE for boundary detection
	dg.EToE = [][]int{{0, 0, 0, 0}}
	dg.EToF = [][]int{{0, 1, 2, 3}}

	elements, faces := dg.GetBoundaryFaces()

	// Should have 4 boundary faces
	if len(elements) != 4 {
		t.Errorf("Boundary faces: got %d, expected 4", len(elements))
	}

	// All should be from element 0
	for i, elem := range elements {
		if elem != 0 {
			t.Errorf("Boundary face %d: element = %d, expected 0", i, elem)
		}
	}

	// Should have faces 0, 1, 2, 3
	for i, face := range faces {
		if face != i {
			t.Errorf("Boundary face %d: face = %d, expected %d", i, face, i)
		}
	}
}

// TestBuildMaps3DSharedFaces tests the GetSharedFaces helper function
func TestBuildMaps3DSharedFaces(t *testing.T) {
	// Use the cube mesh which has many shared faces
	tm := utils.GetStandardTestMeshes()
	VX, VY, VZ, EToV := meshToDG3DFormat(tm.CubeMesh)

	N := 2
	dg, err := NewDG3D(N, VX, VY, VZ, EToV)
	if err != nil {
		t.Fatalf("Failed to create DG3D: %v", err)
	}

	// Get connectivity from mesh
	meshObj := mesh.ConvertToMesh(tm.CubeMesh)
	meshObj.BuildConnectivity()
	dg.EToE = meshObj.EToE
	dg.EToF = meshObj.EToF

	elem1, face1, elem2, face2 := dg.GetSharedFaces()

	// The cube mesh should have many interior faces
	if len(elem1) == 0 {
		t.Error("No shared faces found in cube mesh")
	}

	// Verify that each reported connection is reciprocal
	for i := range elem1 {
		// Check that elem2's face points back to elem1
		if dg.EToE[elem2[i]][face2[i]] != elem1[i] {
			t.Errorf("Non-reciprocal connection at index %d", i)
		}
		if dg.EToF[elem2[i]][face2[i]] != face1[i] {
			t.Errorf("Non-reciprocal face connection at index %d", i)
		}
	}

	// Each pair should be listed only once (elem1 < elem2)
	for i := range elem1 {
		if elem1[i] >= elem2[i] {
			t.Errorf("Shared face %d: elem1 (%d) should be < elem2 (%d)",
				i, elem1[i], elem2[i])
		}
	}
}
