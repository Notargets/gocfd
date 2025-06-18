package gonudg

import (
	"sort"
)

// tiConnect3D method for DG3D struct - calls the standalone tiConnect3D function
func (dg *DG3D) tiConnect3D() {
	dg.EToE, dg.EToF = tiConnect3D(dg.EToV)
}

// tiConnect3D implements tetrahedral face connectivity algorithm due to Toby Isaac
// Returns element-to-element (EToE) and element-to-face (EToF) connectivity arrays
//
// INDEXING NOTE: This is the 0-based version. The original C++ uses 1-based indexing.
// The algorithm:
//  1. Extracts all faces from all tetrahedra (4 faces per tet)
//  2. Sorts nodes within each face for canonical representation
//  3. Assigns unique IDs to faces based on sorted node numbers
//  4. Finds matching faces by comparing IDs
//  5. Updates connectivity arrays with neighbor information
//
// For boundary faces, elements remain self-connected (EToE[k][f] = k, EToF[k][f] = f)
func tiConnect3D(EToV [][]int) (EToE, EToF [][]int) {
	Nfaces := 4
	K := len(EToV)

	// Find maximum node number to determine array sizes
	Nnodes := 0
	for _, elem := range EToV {
		for _, v := range elem {
			if v+1 > Nnodes {
				Nnodes = v + 1
			}
		}
	}

	// Define face-to-node connectivity for tetrahedron
	// Each face is defined by 3 vertices
	faceVertices := [][]int{
		{0, 1, 2}, // Face 0
		{0, 1, 3}, // Face 1
		{1, 2, 3}, // Face 2
		{0, 2, 3}, // Face 3
	}

	// Create list of all faces
	type face struct {
		nodes    [3]int // sorted node indices
		elem     int    // element index
		faceNum  int    // local face number (0-3)
		globalID int64  // unique face identifier
	}

	faces := make([]face, K*Nfaces)

	// Extract all faces from all elements
	faceIdx := 0
	for k := 0; k < K; k++ {
		for f := 0; f < Nfaces; f++ {
			// Get the three vertices of this face
			var faceNodes [3]int
			for i := 0; i < 3; i++ {
				faceNodes[i] = EToV[k][faceVertices[f][i]]
			}

			// Sort nodes to create canonical ordering
			sortedNodes := faceNodes
			sort.Ints(sortedNodes[:])

			// Create unique ID for this face based on node numbers
			// Using the formula from C++: id = n1*Nnodes^2 + n2*Nnodes + n3
			globalID := int64(sortedNodes[0])*int64(Nnodes)*int64(Nnodes) +
				int64(sortedNodes[1])*int64(Nnodes) +
				int64(sortedNodes[2])

			faces[faceIdx] = face{
				nodes:    sortedNodes,
				elem:     k,
				faceNum:  f,
				globalID: globalID,
			}
			faceIdx++
		}
	}

	// Sort faces by their global ID
	sort.Slice(faces, func(i, j int) bool {
		return faces[i].globalID < faces[j].globalID
	})

	// Initialize connectivity arrays with self-connectivity
	EToE = make([][]int, K)
	EToF = make([][]int, K)
	for k := 0; k < K; k++ {
		EToE[k] = make([]int, Nfaces)
		EToF[k] = make([]int, Nfaces)
		for f := 0; f < Nfaces; f++ {
			EToE[k][f] = k
			EToF[k][f] = f
		}
	}

	// Find matching faces
	for i := 0; i < len(faces)-1; i++ {
		if faces[i].globalID == faces[i+1].globalID {
			// Found a match - two faces share the same nodes
			elem1 := faces[i].elem
			face1 := faces[i].faceNum
			elem2 := faces[i+1].elem
			face2 := faces[i+1].faceNum

			// Update connectivity arrays
			EToE[elem1][face1] = elem2
			EToF[elem1][face1] = face2
			EToE[elem2][face2] = elem1
			EToF[elem2][face2] = face1
		}
	}

	return EToE, EToF
}
