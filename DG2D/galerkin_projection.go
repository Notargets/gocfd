package DG2D

import "github.com/notargets/gocfd/utils"

type GalerkinProjection struct {
	FromBasis  *JacobiBasis2D
	InterpFrom utils.Matrix // Interpolation matrix to quadrature points from
	// FromBasis
	ToBasis          *JacobiBasis2D
	Cub              *Cubature
	MassMatrix, Minv utils.Matrix
}

func NewGalerkinProjection(FromBasis, ToBasis *JacobiBasis2D) (gp *GalerkinProjection) {
	gp = &GalerkinProjection{
		FromBasis: FromBasis,
		ToBasis:   ToBasis,
		Cub:       NewCubature(FromBasis.P),
	}
	gp.InterpFrom = FromBasis.GetInterpMatrix(gp.Cub.R, gp.Cub.S)

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
	gp.MassMatrix = utils.NewMatrix(ToBasis.Np, ToBasis.Np)
	for i := 0; i < ToBasis.Np; i++ {
		for j := 0; j < ToBasis.Np; j++ {
			for q := 0; q < gp.Cub.Nq; q++ {
				r_q, s_q, w_q := gp.Cub.R.AtVec(q), gp.Cub.S.AtVec(q), gp.Cub.W.AtVec(q)
				b_i := gp.ToBasis.GetOrthogonalPolynomialAtJ(r_q, s_q, i)
				b_j := gp.ToBasis.GetOrthogonalPolynomialAtJ(r_q, s_q, j)
				mm := gp.MassMatrix.At(i, j) + b_i*b_j*w_q
				gp.MassMatrix.Set(i, j, mm)
			}
		}
	}
	gp.Minv = gp.MassMatrix.InverseWithCheck()
	return
}

func (gp *GalerkinProjection) GetCoefficients(Uh utils.Matrix) (Coeffs []float64) {
	// Coeffs are the coefficients of the target polynomial
	var (
		Np, _ = Uh.Dims()
	)
	// Uh is the solution value vector from the From Basis
	if Np != gp.FromBasis.Np {
		panic("Uh should have length of total node points of FromBasis")
	}
	UhAtQ := gp.InterpFrom.Mul(Uh)
	B := utils.NewMatrix(gp.ToBasis.Np, 1)
	for i := 0; i < gp.ToBasis.Np; i++ {
		for q := 0; q < gp.Cub.Nq; q++ {
			r_q, s_q, w_q := gp.Cub.R.AtVec(q), gp.Cub.S.AtVec(q), gp.Cub.W.AtVec(q)
			b_i := gp.ToBasis.GetOrthogonalPolynomialAtJ(r_q, s_q, i)
			uh := UhAtQ.At(q, 0)
			bb := B.At(i, 0) + b_i*uh*w_q
			B.Set(i, 0, bb)
		}
	}
	Coeffs = gp.Minv.Mul(B).DataP
	return
}

func (gp *GalerkinProjection) GetProjectedValues(Uh utils.Matrix, r, s []float64, v *[]float64) {
	if len(*v) != len(s) {
		panic("value length should match input points length")
	}
	Coeffs := gp.GetCoefficients(Uh)
	for ii, rr := range r {
		ss := s[ii]
		(*v)[ii] = 0.
		for j := 0; j < gp.ToBasis.Np; j++ {
			cc := Coeffs[j]
			(*v)[ii] += cc * gp.ToBasis.GetOrthogonalPolynomialAtJ(rr, ss, j)
		}
	}
	return
}
