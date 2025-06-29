package gonudg

import (
	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/utils"
)

// Lift3D computes the 3D surface to volume lift operator used in DG formulation
// Purpose: Compute 3D surface to volume lift operator used in DG formulation
func (dg *DG3D) Lift3D() error {
	Np := dg.Np
	Nfp := dg.Nfp
	Nfaces := dg.Nfaces

	// Initialize Emat - this will hold face mass matrices
	Emat := utils.NewMatrix(Np, Nfaces*Nfp)

	// Process each of the 4 faces
	for face := 0; face < Nfaces; face++ {
		// Extract the appropriate 2D coordinates for this face
		var faceR, faceS utils.Vector

		switch face {
		case 0: // Face 1: T = -1, use (R, S)
			faceR = utils.NewVector(len(dg.Fmask[face]))
			faceS = utils.NewVector(len(dg.Fmask[face]))
			for i, idx := range dg.Fmask[face] {
				faceR.Set(i, dg.R[idx])
				faceS.Set(i, dg.S[idx])
			}

		case 1: // Face 2: S = -1, use (R, T)
			faceR = utils.NewVector(len(dg.Fmask[face]))
			faceS = utils.NewVector(len(dg.Fmask[face]))
			for i, idx := range dg.Fmask[face] {
				faceR.Set(i, dg.R[idx])
				faceS.Set(i, dg.T[idx])
			}

		case 2: // Face 3: R+S+T = -1, use (S, T)
			faceR = utils.NewVector(len(dg.Fmask[face]))
			faceS = utils.NewVector(len(dg.Fmask[face]))
			for i, idx := range dg.Fmask[face] {
				faceR.Set(i, dg.S[idx])
				faceS.Set(i, dg.T[idx])
			}

		case 3: // Face 4: R = -1, use (S, T)
			faceR = utils.NewVector(len(dg.Fmask[face]))
			faceS = utils.NewVector(len(dg.Fmask[face]))
			for i, idx := range dg.Fmask[face] {
				faceR.Set(i, dg.S[idx])
				faceS.Set(i, dg.T[idx])
			}
		}

		// Compute 2D Vandermonde matrix for the face
		VFace := DG2D.Vandermonde2D(dg.N, faceR, faceS)

		// Compute face mass matrix: massFace = inv(VFace * VFace^T)
		massFace := VFace.Mul(VFace.Transpose()).InverseWithCheck()

		// Place face mass matrix into Emat at the appropriate location
		// The C++ code does: Emat(idr, JJ) = massFace
		// This sets a block of Emat where:
		// - rows are the face node indices (Fmask[face])
		// - columns are [face*Nfp : (face+1)*Nfp)
		for i, nodeIdx := range dg.Fmask[face] {
			for j := 0; j < Nfp; j++ {
				colIdx := face*Nfp + j
				Emat.Set(nodeIdx, colIdx, massFace.At(i, j))
			}
		}
	}

	// Compute LIFT = V * (V^T * Emat)
	VtE := dg.V.Transpose().Mul(Emat)
	dg.LIFT = dg.V.Mul(VtE)

	return nil
}
