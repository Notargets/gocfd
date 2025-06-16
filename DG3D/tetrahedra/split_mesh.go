package tetrahedra

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
)

// SplitByPartition splits the Element3D into separate Element3D instances for each partition
// This should be called after EToP is populated and before face_buffer or other parallel processing
func (el *Element3D) SplitByPartition() error {
	if el.EToP == nil {
		// Nothing to split - non-partitioned mesh
		return nil
	}

	// Count elements per partition and create partition mapping
	partitionCount := make(map[int]int)
	partitionElements := make(map[int][]int) // partition -> list of element indices

	for elemID, partID := range el.EToP {
		partitionCount[partID]++
		partitionElements[partID] = append(partitionElements[partID], elemID)
	}

	// Create split Element3D for each partition
	el.SplitElement3D = make([]*Element3D, 0, len(partitionCount))

	for partID, elemList := range partitionElements {
		splitEl, err := el.createPartitionElement(partID, elemList)
		if err != nil {
			return fmt.Errorf("failed to create partition %d: %v", partID, err)
		}
		el.SplitElement3D = append(el.SplitElement3D, splitEl)
	}

	return nil
}

// createPartitionElement creates a new Element3D for a specific partition
func (el *Element3D) createPartitionElement(partitionID int, elemIndices []int) (*Element3D, error) {
	localK := len(elemIndices)

	// Create global to local element mapping
	globalToLocal := make(map[int]int)
	for localIdx, globalIdx := range elemIndices {
		globalToLocal[globalIdx] = localIdx
	}

	// Create new Element3D for this partition
	partEl := &Element3D{
		K:        localK,
		TetBasis: el.TetBasis, // Share the basis
		Mesh:     el.Mesh,     // Reference to original mesh
	}

	// Allocate vertex coordinates for local elements
	partEl.VX = utils.NewVector(localK * 4)
	partEl.VY = utils.NewVector(localK * 4)
	partEl.VZ = utils.NewVector(localK * 4)

	// Copy vertex coordinates for local elements
	for localIdx, globalIdx := range elemIndices {
		for v := 0; v < 4; v++ {
			partEl.VX.Set(localIdx*4+v, el.VX.At(globalIdx*4+v))
			partEl.VY.Set(localIdx*4+v, el.VY.At(globalIdx*4+v))
			partEl.VZ.Set(localIdx*4+v, el.VZ.At(globalIdx*4+v))
		}
	}

	// Allocate and copy physical coordinates
	partEl.X = utils.NewMatrix(el.Np, localK)
	partEl.Y = utils.NewMatrix(el.Np, localK)
	partEl.Z = utils.NewMatrix(el.Np, localK)

	for localIdx, globalIdx := range elemIndices {
		for n := 0; n < el.Np; n++ {
			partEl.X.Set(n, localIdx, el.X.At(n, globalIdx))
			partEl.Y.Set(n, localIdx, el.Y.At(n, globalIdx))
			partEl.Z.Set(n, localIdx, el.Z.At(n, globalIdx))
		}
	}

	// Copy face coordinates
	if el.Fx != nil {
		partEl.Fx = make([]utils.Matrix, len(el.Fx))
		partEl.Fy = make([]utils.Matrix, len(el.Fy))
		partEl.Fz = make([]utils.Matrix, len(el.Fz))

		for f := 0; f < len(el.Fx); f++ {
			partEl.Fx[f] = utils.NewMatrix(el.Nfp, localK)
			partEl.Fy[f] = utils.NewMatrix(el.Nfp, localK)
			partEl.Fz[f] = utils.NewMatrix(el.Nfp, localK)

			for localIdx, globalIdx := range elemIndices {
				for i := 0; i < el.Nfp; i++ {
					partEl.Fx[f].Set(i, localIdx, el.Fx[f].At(i, globalIdx))
					partEl.Fy[f].Set(i, localIdx, el.Fy[f].At(i, globalIdx))
					partEl.Fz[f].Set(i, localIdx, el.Fz[f].At(i, globalIdx))
				}
			}
		}
	}

	// Copy and remap connectivity
	partEl.EToV = make([][]int, localK)
	for localIdx, globalIdx := range elemIndices {
		partEl.EToV[localIdx] = make([]int, len(el.EToV[globalIdx]))
		copy(partEl.EToV[localIdx], el.EToV[globalIdx])
	}

	// Set EToP - all elements in this partition have the same partition ID
	partEl.EToP = make([]int, localK)
	for i := range partEl.EToP {
		partEl.EToP[i] = partitionID
	}

	// Copy and remap connectivity arrays
	if el.ConnectivityArrays != nil {
		partEl.ConnectivityArrays = &ConnectivityArrays{
			EToE: make([][]int, localK),
			EToF: make([][]int, localK),
		}

		for localIdx, globalIdx := range elemIndices {
			partEl.EToE[localIdx] = make([]int, len(el.EToE[globalIdx]))
			partEl.EToF[localIdx] = make([]int, len(el.EToF[globalIdx]))

			for f := 0; f < len(el.EToE[globalIdx]); f++ {
				neighborGlobal := el.EToE[globalIdx][f]

				if neighborGlobal == -1 {
					// Boundary face
					partEl.EToE[localIdx][f] = -1
					partEl.EToF[localIdx][f] = el.EToF[globalIdx][f]
				} else if localNeighbor, isLocal := globalToLocal[neighborGlobal]; isLocal {
					// Local neighbor within this partition
					partEl.EToE[localIdx][f] = localNeighbor
					partEl.EToF[localIdx][f] = el.EToF[globalIdx][f]
				} else {
					// Remote neighbor in different partition
					// Mark as boundary for now - face_buffer will handle remote connections
					partEl.EToE[localIdx][f] = -1
					partEl.EToF[localIdx][f] = el.EToF[globalIdx][f]
				}
			}
		}
	}

	// Copy boundary conditions
	if el.BCType != nil {
		nfaces := 4 // Tetrahedra
		partEl.BCType = make([]int, localK*nfaces)

		for localIdx, globalIdx := range elemIndices {
			for f := 0; f < nfaces; f++ {
				if globalIdx*nfaces+f < len(el.BCType) {
					partEl.BCType[localIdx*nfaces+f] = el.BCType[globalIdx*nfaces+f]
				}
			}
		}
	}

	// Copy geometric factors
	if el.GeometricFactors != nil {
		partEl.GeometricFactors = &GeometricFactors{
			Rx: utils.NewMatrix(el.Np, localK),
			Ry: utils.NewMatrix(el.Np, localK),
			Rz: utils.NewMatrix(el.Np, localK),
			Sx: utils.NewMatrix(el.Np, localK),
			Sy: utils.NewMatrix(el.Np, localK),
			Sz: utils.NewMatrix(el.Np, localK),
			Tx: utils.NewMatrix(el.Np, localK),
			Ty: utils.NewMatrix(el.Np, localK),
			Tz: utils.NewMatrix(el.Np, localK),
			J:  utils.NewMatrix(el.Np, localK),
		}

		for localIdx, globalIdx := range elemIndices {
			for n := 0; n < el.Np; n++ {
				partEl.Rx.Set(n, localIdx, el.Rx.At(n, globalIdx))
				partEl.Ry.Set(n, localIdx, el.Ry.At(n, globalIdx))
				partEl.Rz.Set(n, localIdx, el.Rz.At(n, globalIdx))
				partEl.Sx.Set(n, localIdx, el.Sx.At(n, globalIdx))
				partEl.Sy.Set(n, localIdx, el.Sy.At(n, globalIdx))
				partEl.Sz.Set(n, localIdx, el.Sz.At(n, globalIdx))
				partEl.Tx.Set(n, localIdx, el.Tx.At(n, globalIdx))
				partEl.Ty.Set(n, localIdx, el.Ty.At(n, globalIdx))
				partEl.Tz.Set(n, localIdx, el.Tz.At(n, globalIdx))
				partEl.J.Set(n, localIdx, el.J.At(n, globalIdx))
			}
		}
	}

	// Copy face geometric factors
	if el.FaceGeometricFactors != nil {
		nfp := el.Nfp
		nfaces := 4
		totalFaceNodes := nfp * nfaces

		partEl.FaceGeometricFactors = &FaceGeometricFactors{
			Nx:     utils.NewMatrix(totalFaceNodes, localK),
			Ny:     utils.NewMatrix(totalFaceNodes, localK),
			Nz:     utils.NewMatrix(totalFaceNodes, localK),
			SJ:     utils.NewMatrix(totalFaceNodes, localK),
			Fscale: utils.NewMatrix(totalFaceNodes, localK),
		}

		for localIdx, globalIdx := range elemIndices {
			for i := 0; i < totalFaceNodes; i++ {
				partEl.Nx.Set(i, localIdx, el.Nx.At(i, globalIdx))
				partEl.Ny.Set(i, localIdx, el.Ny.At(i, globalIdx))
				partEl.Nz.Set(i, localIdx, el.Nz.At(i, globalIdx))
				partEl.SJ.Set(i, localIdx, el.SJ.At(i, globalIdx))
				partEl.Fscale.Set(i, localIdx, el.Fscale.At(i, globalIdx))
			}
		}
	}

	// Rebuild face mappings for the local partition
	partEl.BuildMaps3D()

	// Store original element indices for reference
	partEl.tetIndices = make([]int, localK)
	copy(partEl.tetIndices, elemIndices)

	return partEl, nil
}

// GetPartition returns the Element3D for a specific partition ID
func (el *Element3D) GetPartition(partitionID int) (*Element3D, error) {
	if el.SplitElement3D == nil {
		return nil, fmt.Errorf("mesh has not been split by partition")
	}

	// Find partition with matching ID
	for _, partEl := range el.SplitElement3D {
		if partEl.EToP != nil && len(partEl.EToP) > 0 && partEl.EToP[0] == partitionID {
			return partEl, nil
		}
	}

	return nil, fmt.Errorf("partition %d not found", partitionID)
}

// GetOriginalElementIndex returns the global element index for a local element in a partition
func (el *Element3D) GetOriginalElementIndex(localIdx int) (int, error) {
	if el.tetIndices == nil {
		return -1, fmt.Errorf("this is not a partitioned Element3D")
	}

	if localIdx < 0 || localIdx >= len(el.tetIndices) {
		return -1, fmt.Errorf("local index %d out of range [0, %d)", localIdx, len(el.tetIndices))
	}

	return el.tetIndices[localIdx], nil
}
