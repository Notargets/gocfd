package gonudg

// BuildMaps3D builds connectivity and boundary tables for nodes
// This populates vmapM, vmapP, mapM, mapP, mapB, and vmapB in the DG3D struct
func (dg *DG3D) BuildMaps3D() {
	// Initialize arrays
	NF := dg.Nfp * dg.Nfaces
	dg.vmapM = make([]int, NF*dg.K)
	dg.vmapP = make([]int, NF*dg.K)
	dg.mapM = make([]int, NF*dg.K)
	dg.mapP = make([]int, NF*dg.K)

	// Initialize mapM and mapP as sequential indices
	for i := 0; i < len(dg.mapM); i++ {
		dg.mapM[i] = i
		dg.mapP[i] = i
	}

	// Number volume nodes consecutively
	// Find index of face nodes with respect to volume node ordering
	for k1 := 0; k1 < dg.K; k1++ {
		iL1 := k1 * NF
		for f := 0; f < dg.Nfaces; f++ {
			for i := 0; i < dg.Nfp; i++ {
				idsL := iL1 + f*dg.Nfp + i
				// Fmask contains node indices within element (0-based)
				volumeNode := dg.Fmask[f][i] + k1*dg.Np
				dg.vmapM[idsL] = volumeNode
			}
		}
	}

	// Initialize vmapP to vmapM (for boundary faces)
	copy(dg.vmapP, dg.vmapM)

	// Create face to face mapping if connectivity information exists
	if dg.EToE != nil && dg.EToF != nil {
		dg.buildFaceConnectivity()
	}

	// Find boundary nodes
	dg.buildBoundaryMaps()
}

// buildFaceConnectivity creates face-to-face mappings for interior faces
func (dg *DG3D) buildFaceConnectivity() {
	NF := dg.Nfp * dg.Nfaces

	// Extract coordinate data into flat arrays for compatibility
	x := make([]float64, dg.Np*dg.K)
	y := make([]float64, dg.Np*dg.K)
	z := make([]float64, dg.Np*dg.K)

	for k := 0; k < dg.K; k++ {
		for i := 0; i < dg.Np; i++ {
			idx := k*dg.Np + i
			x[idx] = dg.x.At(i, k)
			y[idx] = dg.y.At(i, k)
			z[idx] = dg.z.At(i, k)
		}
	}

	// Create face to face mapping
	for k1 := 0; k1 < dg.K; k1++ {
		for f1 := 0; f1 < dg.Nfaces; f1++ {
			// Find neighbor
			k2 := dg.EToE[k1][f1]
			f2 := dg.EToF[k1][f1]

			// Skip boundary faces (self-connected)
			if k2 == k1 && f2 == f1 {
				continue
			}

			// Skip invalid neighbors
			if k2 < 0 || k2 >= dg.K {
				continue
			}

			skM := k1 * NF // offset to element k1
			skP := k2 * NF // offset to element k2

			// Build index lists for faces
			idsM := make([]int, dg.Nfp)
			idsP := make([]int, dg.Nfp)
			for i := 0; i < dg.Nfp; i++ {
				idsM[i] = f1*dg.Nfp + i + skM
				idsP[i] = f2*dg.Nfp + i + skP
			}

			// Find volume node numbers of left and right nodes
			vidM := make([]int, dg.Nfp)
			vidP := make([]int, dg.Nfp)
			for i := 0; i < dg.Nfp; i++ {
				vidM[i] = dg.vmapM[idsM[i]]
				vidP[i] = dg.vmapM[idsP[i]]
			}

			// Extract coordinates
			x1 := make([]float64, dg.Nfp)
			y1 := make([]float64, dg.Nfp)
			z1 := make([]float64, dg.Nfp)
			x2 := make([]float64, dg.Nfp)
			y2 := make([]float64, dg.Nfp)
			z2 := make([]float64, dg.Nfp)

			for i := 0; i < dg.Nfp; i++ {
				// vidM[i] is the global node index
				x1[i] = x[vidM[i]]
				y1[i] = y[vidM[i]]
				z1[i] = z[vidM[i]]

				x2[i] = x[vidP[i]]
				y2[i] = y[vidP[i]]
				z2[i] = z[vidP[i]]
			}

			// Compute distance matrix
			D := make([][]float64, dg.Nfp)
			for i := range D {
				D[i] = make([]float64, dg.Nfp)
			}

			for i := 0; i < dg.Nfp; i++ {
				for j := 0; j < dg.Nfp; j++ {
					dx := x2[j] - x1[i]
					dy := y2[j] - y1[i]
					dz := z2[j] - z1[i]
					D[i][j] = dx*dx + dy*dy + dz*dz
				}
			}

			// Find matching nodes
			for i := 0; i < dg.Nfp; i++ {
				// Find minimum distance
				minDist := D[i][0]
				minJ := 0
				for j := 1; j < dg.Nfp; j++ {
					if D[i][j] < minDist {
						minDist = D[i][j]
						minJ = j
					}
				}

				// Map to neighbor if within tolerance
				if minDist < dg.NODETOL*dg.NODETOL {
					dg.vmapP[idsM[i]] = vidP[minJ]
					dg.mapP[idsM[i]] = idsP[minJ]
				}
			}
		}
	}
}

// buildBoundaryMaps finds and stores boundary node information
func (dg *DG3D) buildBoundaryMaps() {
	// Create lists of boundary nodes
	mapB := make([]int, 0)
	vmapB := make([]int, 0)

	// Check each face node
	NF := dg.Nfp * dg.Nfaces
	for i := 0; i < dg.K*NF; i++ {
		if dg.mapP[i] == i {
			// This is a boundary node
			mapB = append(mapB, i)
			vmapB = append(vmapB, dg.vmapM[i])
		}
	}

	dg.mapB = mapB
	dg.vmapB = vmapB
}

// GetBoundaryFaces returns the element and face indices for all boundary faces
func (dg *DG3D) GetBoundaryFaces() (elements []int, faces []int) {
	elements = make([]int, 0)
	faces = make([]int, 0)

	if dg.EToE == nil {
		return
	}

	for k := 0; k < dg.K; k++ {
		for f := 0; f < dg.Nfaces; f++ {
			// Check if this is a boundary face
			if dg.EToE[k][f] == k && dg.EToF[k][f] == f {
				elements = append(elements, k)
				faces = append(faces, f)
			}
		}
	}

	return
}

// GetSharedFaces returns pairs of elements/faces that share a face
func (dg *DG3D) GetSharedFaces() (elem1 []int, face1 []int, elem2 []int, face2 []int) {
	elem1 = make([]int, 0)
	face1 = make([]int, 0)
	elem2 = make([]int, 0)
	face2 = make([]int, 0)

	if dg.EToE == nil {
		return
	}

	// To avoid duplicates, only record when k1 < k2
	for k1 := 0; k1 < dg.K; k1++ {
		for f1 := 0; f1 < dg.Nfaces; f1++ {
			k2 := dg.EToE[k1][f1]
			f2 := dg.EToF[k1][f1]

			// Skip boundary faces
			if k2 == k1 && f2 == f1 {
				continue
			}

			// Only record if k1 < k2 to avoid duplicates
			if k1 < k2 {
				elem1 = append(elem1, k1)
				face1 = append(face1, f1)
				elem2 = append(elem2, k2)
				face2 = append(face2, f2)
			}
		}
	}

	return
}
