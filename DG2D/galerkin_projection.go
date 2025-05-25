package DG2D

import (
	"github.com/notargets/gocfd/utils"
)

type GalerkinProjection struct {
	FromBasis *JacobiBasis2D
	ToBasis   *JacobiBasis2D
	Cub       *Cubature
	A         utils.Matrix // Precomputed projection matrix
}

func NewGalerkinProjection(FromBasis, ToBasis *JacobiBasis2D) (gp *GalerkinProjection) {
	gp = &GalerkinProjection{
		FromBasis: FromBasis,
		ToBasis:   ToBasis,
		Cub:       NewCubature(FromBasis.P),
	}
	var (
		Nq = gp.Cub.Nq
	)
	// Compose the MassMatrix, which is re-used for projections
	/*
		for i := 0; i < N_b; i++ {
			for j := 0; j < N_b; j++ {
				M[i][j] = 0.0
				for q := 0; q < N_c; q++ {
					psi_i := ψ_i(R[q], S[q])
					psi_j := ψ_j(R[q], S[q])
					M[i][j] += psi_i * psi_j * W[q]
				}
			}
		}
	*/
	// An interpolation matrix for the projected low order basis P is:
	// P = ILowAtRS x Minv x ILowQ^T x W x IHighQ
	// We store A = Minv x ILowQ^T x W x IHighQ
	MassMatrix := utils.NewMatrix(ToBasis.Np, ToBasis.Np)
	for i := 0; i < ToBasis.Np; i++ {
		for j := 0; j < ToBasis.Np; j++ {
			for q := 0; q < gp.Cub.Nq; q++ {
				r_q, s_q, w_q := gp.Cub.R.AtVec(q), gp.Cub.S.AtVec(q), gp.Cub.W.AtVec(q)
				b_i := gp.ToBasis.GetOrthogonalPolynomialAtJ(r_q, s_q, i)
				b_j := gp.ToBasis.GetOrthogonalPolynomialAtJ(r_q, s_q, j)
				mm := MassMatrix.At(i, j) + b_i*b_j*w_q
				MassMatrix.Set(i, j, mm)
			}
		}
	}
	Minv := MassMatrix.InverseWithCheck()
	IHighQ := gp.FromBasis.GetInterpMatrix(gp.Cub.R, gp.Cub.S)
	ILowQ := gp.ToBasis.GetInterpMatrix(gp.Cub.R, gp.Cub.S)
	W := utils.NewDiagMatrix(Nq, gp.Cub.W.DataP)
	gp.A = Minv.Mul(ILowQ.Transpose().Mul(W).Mul(IHighQ))
	return
}

func (gp *GalerkinProjection) GetProjectedInterpolationMatrix(R, S utils.Vector) (I utils.Matrix) {
	ILowRS := gp.ToBasis.GetInterpMatrix(R, S)
	I = ILowRS.Mul(gp.A)
	return
}
