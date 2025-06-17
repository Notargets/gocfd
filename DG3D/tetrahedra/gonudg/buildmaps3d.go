package gonudg

// INDEXING NOTE: Original C++ code uses 1-based indexing to emulate Matlab behavior.
// This Go port uses standard 0-based indexing. Example conversions:
//   C++: sk = 1; V3D(All,sk) = ...    ->    Go: sk = 0; V3D.SetCol(sk, ...)
//   C++: Fmask[1] (first face)        ->    Go: Fmask[0] (first face)
// The indexing has been correctly translated throughout this port.

import (
	"math"
)

// BuildMaps3D builds connectivity and boundary tables for nodes
// This is the 0-based index version of the C++ BuildMaps3D function
// Returns vmapM, vmapP, mapB, vmapB
func BuildMaps3D(K, Np, Nfp int, Nfaces int, x, y, z []float64,
	EToE, EToF [][]int, Fmask [][]int) (vmapM, vmapP, mapB, vmapB []int) {

	NODETOL := 1e-7
	NF := Nfp * Nfaces

	// Initialize arrays
	vmapM = make([]int, Nfp*Nfaces*K)
	vmapP = make([]int, Nfp*Nfaces*K)
	mapM := make([]int, Nfp*Nfaces*K)
	mapP := make([]int, Nfp*Nfaces*K)

	// Initialize mapM and mapP as sequential indices
	for i := 0; i < len(mapM); i++ {
		mapM[i] = i
		mapP[i] = i
	}

	// Number volume nodes consecutively
	// Find index of face nodes with respect to volume node ordering
	for k1 := 0; k1 < K; k1++ {
		iL1 := k1 * NF
		for f := 0; f < Nfaces; f++ {
			for i := 0; i < Nfp; i++ {
				idsL := iL1 + f*Nfp + i
				// Fmask contains node indices within element (0-based)
				volumeNode := Fmask[f][i] + k1*Np
				vmapM[idsL] = volumeNode
			}
		}
	}

	// Initialize vmapP to vmapM (for boundary faces)
	copy(vmapP, vmapM)

	// Create face to face mapping
	for k1 := 0; k1 < K; k1++ {
		for f1 := 0; f1 < Nfaces; f1++ {
			// Find neighbor
			k2 := EToE[k1][f1]
			f2 := EToF[k1][f1]

			// Skip boundary faces (self-connected)
			if k2 == k1 && f2 == f1 {
				continue
			}

			// Skip invalid neighbors
			if k2 < 0 || k2 >= K {
				continue
			}

			skM := k1 * NF // offset to element k1
			skP := k2 * NF // offset to element k2

			// Build index lists for faces
			idsM := make([]int, Nfp)
			idsP := make([]int, Nfp)
			for i := 0; i < Nfp; i++ {
				idsM[i] = f1*Nfp + i + skM
				idsP[i] = f2*Nfp + i + skP
			}

			// Find volume node numbers of left and right nodes
			vidM := make([]int, Nfp)
			vidP := make([]int, Nfp)
			for i := 0; i < Nfp; i++ {
				vidM[i] = vmapM[idsM[i]]
				vidP[i] = vmapM[idsP[i]]
			}

			// Extract coordinates
			x1 := make([]float64, Nfp)
			y1 := make([]float64, Nfp)
			z1 := make([]float64, Nfp)
			x2 := make([]float64, Nfp)
			y2 := make([]float64, Nfp)
			z2 := make([]float64, Nfp)

			for i := 0; i < Nfp; i++ {
				// vidM[i] is the global node index
				x1[i] = x[vidM[i]]
				y1[i] = y[vidM[i]]
				z1[i] = z[vidM[i]]

				x2[i] = x[vidP[i]]
				y2[i] = y[vidP[i]]
				z2[i] = z[vidP[i]]
			}

			// Compute distance matrix and find matches
			for i := 0; i < Nfp; i++ {
				for j := 0; j < Nfp; j++ {
					dx := x1[i] - x2[j]
					dy := y1[i] - y2[j]
					dz := z1[i] - z2[j]
					dist := math.Sqrt(dx*dx + dy*dy + dz*dz)

					if dist < NODETOL {
						// Found matching nodes
						idM := idsM[i]
						vmapP[idM] = vidP[j] // Set external element node

						idP := idsP[j]
						mapP[idM] = idP // Set external face node
						break           // Only one match per node
					}
				}
			}
		}
	}

	// Create list of boundary nodes
	mapB = []int{}
	vmapB = []int{}

	for i := 0; i < len(vmapP); i++ {
		if vmapP[i] == vmapM[i] {
			mapB = append(mapB, i)
			vmapB = append(vmapB, vmapM[i])
		}
	}

	return vmapM, vmapP, mapB, vmapB
}

// Helper function to find matching nodes between two faces
func findFaceMatches(x1, y1, z1, x2, y2, z2 []float64, tol float64) [][]int {
	n1 := len(x1)
	n2 := len(x2)
	matches := [][]int{}

	for i := 0; i < n1; i++ {
		for j := 0; j < n2; j++ {
			dx := x1[i] - x2[j]
			dy := y1[i] - y2[j]
			dz := z1[i] - z2[j]
			dist := math.Sqrt(dx*dx + dy*dy + dz*dz)

			if dist < tol {
				matches = append(matches, []int{i, j})
			}
		}
	}

	return matches
}
