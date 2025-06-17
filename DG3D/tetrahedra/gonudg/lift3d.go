package gonudg

import (
	"github.com/notargets/gocfd/DG2D"
	"github.com/notargets/gocfd/utils"
)

// Lift3D computes the 3D surface to volume lift operator used in DG formulation
// Purpose: Compute 3D surface to volume lift operator used in DG formulation
func Lift3D(N int, r, s, t []float64, V utils.Matrix, Fmask [][]int) utils.Matrix {
	Np := len(r)
	Nfp := (N + 1) * (N + 2) / 2
	Nfaces := 4

	// Initialize Emat - this will hold face mass matrices
	Emat := utils.NewMatrix(Np, Nfaces*Nfp)

	// Process each of the 4 faces
	for face := 0; face < Nfaces; face++ {
		// Extract the appropriate 2D coordinates for this face
		var faceR, faceS utils.Vector

		switch face {
		case 0: // Face 1: t = -1, use (r, s)
			faceR = utils.NewVector(len(Fmask[face]))
			faceS = utils.NewVector(len(Fmask[face]))
			for i, idx := range Fmask[face] {
				faceR.Set(i, r[idx])
				faceS.Set(i, s[idx])
			}

		case 1: // Face 2: s = -1, use (r, t)
			faceR = utils.NewVector(len(Fmask[face]))
			faceS = utils.NewVector(len(Fmask[face]))
			for i, idx := range Fmask[face] {
				faceR.Set(i, r[idx])
				faceS.Set(i, t[idx])
			}

		case 2: // Face 3: r+s+t = -1, use (s, t)
			faceR = utils.NewVector(len(Fmask[face]))
			faceS = utils.NewVector(len(Fmask[face]))
			for i, idx := range Fmask[face] {
				faceR.Set(i, s[idx])
				faceS.Set(i, t[idx])
			}

		case 3: // Face 4: r = -1, use (s, t)
			faceR = utils.NewVector(len(Fmask[face]))
			faceS = utils.NewVector(len(Fmask[face]))
			for i, idx := range Fmask[face] {
				faceR.Set(i, s[idx])
				faceS.Set(i, t[idx])
			}
		}

		// Compute 2D Vandermonde matrix for the face
		VFace := DG2D.Vandermonde2D(N, faceR, faceS)

		// Compute face mass matrix: massFace = inv(VFace * VFace^T)
		massFace := VFace.Mul(VFace.Transpose()).InverseWithCheck()

		// Place face mass matrix into Emat at the appropriate location
		// The C++ code does: Emat(idr, JJ) = massFace
		// This sets a block of Emat where:
		// - rows are the face node indices (Fmask[face])
		// - columns are [face*Nfp : (face+1)*Nfp)
		for i, nodeIdx := range Fmask[face] {
			for j := 0; j < Nfp; j++ {
				colIdx := face*Nfp + j
				Emat.Set(nodeIdx, colIdx, massFace.At(i, j))
			}
		}
	}

	// Compute LIFT = V * (V^T * Emat)
	VtE := V.Transpose().Mul(Emat)
	LIFT := V.Mul(VtE)

	return LIFT
}
