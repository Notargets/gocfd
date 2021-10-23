package DG2D

import (
	"fmt"

	"github.com/notargets/gocfd/utils"
)

type LagrangeElement2D struct {
	N, Nfp, Np, NFaces  int
	R, S                utils.Vector
	Dr, Ds              utils.Matrix
	V, Vinv, MassMatrix utils.Matrix
	Cub                 *Cubature
}

type Cubature struct {
	r, s, w                 utils.Vector
	W                       utils.Matrix
	V, Dr, Ds, VT, DrT, DsT utils.Matrix
	x, y, rx, sx, ry, sy, J utils.Matrix
	mm, mmCHOL              utils.Matrix
}

type NodeType string

const (
	Epsilon   = "Epsilon"
	Hesthaven = "Hesthaven"
)

func NewLagrangeElement2D(N int, nodeType NodeType) (el *LagrangeElement2D) {
	var (
		err error
	)
	el = &LagrangeElement2D{
		N:      N,
		Np:     (N + 1) * (N + 2) / 2,
		NFaces: 3,
	}
	if N < 1 {
		panic(fmt.Errorf("Polynomial order must be >= 1, have %d", N))
	}
	el.Nfp = el.N + 1
	el.Np = (el.N + 1) * (el.N + 2) / 2
	el.NFaces = 3
	// Compute nodal set
	switch nodeType {
	case Epsilon:
		el.R, el.S = NodesEpsilon(el.N)
	case Hesthaven:
		el.R, el.S = XYtoRS(Nodes2D(el.N))
	}
	// Build reference element matrices
	el.V = Vandermonde2D(el.N, el.R, el.S)
	if el.Vinv, err = el.V.Inverse(); err != nil {
		panic(err)
	}
	el.MassMatrix = el.Vinv.Transpose().Mul(el.Vinv)
	// Initialize the (r,s) differentiation matrices on the simplex, evaluated at (r,s) at order N
	/*
		Vr, Vs := GradVandermonde2D(el.N, el.R, el.S)
		el.Dr = Vr.Mul(el.Vinv)
		el.Ds = Vs.Mul(el.Vinv)
	*/
	el.Dr, el.Ds = el.GetDerivativeMatrices(el.R, el.S)
	// Mark fields read only
	el.V.SetReadOnly("V")
	el.Vinv.SetReadOnly("Vinv")
	el.MassMatrix.SetReadOnly("MassMatrix")
	el.Dr.SetReadOnly("Dr")
	el.Ds.SetReadOnly("Ds")
	return
}

func (el *LagrangeElement2D) GetDerivativeMatrices(R, S utils.Vector) (Dr, Ds utils.Matrix) {
	Vr, Vs := GradVandermonde2D(el.N, R, S)
	Dr, Ds = Vr.Mul(el.Vinv), Vs.Mul(el.Vinv)
	return
}

func (el *LagrangeElement2D) Simplex2DInterpolatingPolyMatrixJacobi(R, S utils.Vector) (polyMatrix utils.Matrix) {
	/*
		Uses Jacobi polynomials as the basis function

		Compose a matrix of interpolating polynomials where each row represents one [r,s] location to be interpolated
		This matrix can then be multiplied by a single vector of function values at the polynomial nodes to produce a
		vector of interpolated values, one for each interpolation location
	*/
	var (
		N  = el.N
		Np = el.Np
	)
	// First compute polynomial terms, used by all polynomials
	polyTerms := make([]float64, R.Len()*Np)
	var sk int
	for ii, r := range R.DataP {
		s := S.DataP[ii]
		for i := 0; i <= N; i++ {
			for j := 0; j <= (N - i); j++ {
				polyTerms[sk] = Simplex2DPTerm(r, s, i, j)
				sk++
			}
		}
	}
	ptV := utils.NewMatrix(R.Len(), Np, polyTerms).Transpose()
	polyMatrix = el.Vinv.Transpose().Mul(ptV).Transpose()
	return
}

func (el *LagrangeElement2D) NewCube2D(COrder int) {
	// function [cubR,cubS,cubW, Ncub] = Cubature2D(COrder)
	// Purpose: provide multidimensional quadrature (i.e. cubature)
	//          rules to integrate up to COrder polynomials

	if COrder > 28 {
		COrder = 28
	}

	if COrder <= 28 {
		cub2d := getCub(COrder)
		nr := len(cub2d) / 3
		cubMat := utils.NewMatrix(nr, 3, cub2d)
		el.Cub = &Cubature{
			r: cubMat.Col(0),
			s: cubMat.Col(1),
			w: cubMat.Col(2),
		}
	} else {
		err := fmt.Errorf("Cubature2D(%d): COrder > 28 not yet tested\n", COrder)
		panic(err)
		/*
		   DVec cuba,cubwa, cubb,cubwb
		   DMat cubA, cubB, cubR, cubS, cubW, tA,tB

		   int cubNA=(int)ceil((COrder+1.0)/2.0)
		   int cubNB=(int)ceil((COrder+1.0)/2.0)


		   JacobiGQ(1.0, 0.0, cubNB-1,  cubb,cubwb)

		   cubA = outer( ones(cubNB), cuba )
		   cubB = outer( cubb, ones(cubNA) )

		   tA = 1.0+cubA
		   tB = 1.0-cubB
		   cubR = 0.5 * tA.dm(tB) - 1.0
		   cubS = cubB
		   cubW = 0.5 * outer(cubwb, cubwa)

		   cub.r = cubR
		   cub.s = cubS
		   cub.w = cubW
		   cub.Ncub = cub.r.size()
		*/
	}
	return
}
