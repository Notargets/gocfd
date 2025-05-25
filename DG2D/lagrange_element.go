package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
)

type LagrangeElement2D struct {
	N, Nfp, Np, NFaces int
	R, S               utils.Vector
	Dr, Ds             utils.Matrix
	MassMatrix         utils.Matrix
	Cub                *Cubature
	JB2D               *JacobiBasis2D
}

type Cubature struct {
	R, S, W utils.Vector
	Nq      int
}

func NewCubature(P int) (cb *Cubature) {
	cub2d := getCub(P)
	Nq := len(cub2d) / 3
	cubMat := utils.NewMatrix(Nq, 3, cub2d)
	cb = &Cubature{
		R:  cubMat.Col(0),
		S:  cubMat.Col(1),
		W:  cubMat.Col(2),
		Nq: Nq,
	}
	return
}

type NodeType string

const (
	// Only the WSJ distribution has Quadrature colocated,
	// we will assume their use throughout
	WSJ       = NodeType("WSJ")
	Epsilon   = NodeType("Epsilon")
	Hesthaven = NodeType("Hesthaven")
	Uniform   = NodeType("Uniform")
)

func NewLagrangeElement2D(N int, nodeType NodeType) (el *LagrangeElement2D) {
	el = &LagrangeElement2D{
		N:      N,
		Np:     (N + 1) * (N + 2) / 2,
		NFaces: 3,
	}
	if N < 0 {
		panic(fmt.Errorf("Polynomial order must be >= 0, have %d", N))
	}
	el.Nfp = el.N + 1
	el.Np = (el.N + 1) * (el.N + 2) / 2
	el.NFaces = 3
	// Only the WSJ distribution has Quadrature colocated,
	// we will assume their use throughout
	switch nodeType {
	case WSJ:
		el.R, el.S = MakeRSFromPoints(WilliamsShunnJameson(el.N))
	case Epsilon:
		el.R, el.S = NodesEpsilon(el.N)
	case Hesthaven:
		el.R, el.S = XYtoRS(Nodes2D(el.N))
	}

	// Build reference element matrices
	el.JB2D = NewJacobiBasis2D(el.N, el.R, el.S, 0, 0)

	w := WilliamsShunnJamesonWeights(el.N)
	W := utils.NewDiagMatrix(len(w), w)
	el.MassMatrix = el.JB2D.V.Transpose().Mul(W).Mul(el.JB2D.V) // el.MassMatrix.Print("M")

	el.Dr, el.Ds = el.GetDerivativeMatrices(el.R, el.S)
	// Mark fields read only
	el.MassMatrix.SetReadOnly("MassMatrix")
	el.Dr.SetReadOnly("Dr")
	el.Ds.SetReadOnly("Ds")
	return
}

func (el *LagrangeElement2D) CheckBasisOrthogonality() {
	var (
		N   = el.N
		Np  = el.Np
		V   = el.JB2D.V
		W   = utils.NewDiagMatrix(Np, WilliamsShunnJamesonWeights(N))
		tol = 1.e-12
	)
	D := V.Transpose().Mul(W).Mul(V)
	// D.Print("D at order:" + strconv.Itoa(N))
	// Compute Frobenius norm of D - I
	var frobNorm float64
	for i := 0; i < Np; i++ {
		for j := 0; j < Np; j++ {
			delta := 0.0
			if i == j {
				delta = 1.0
			}
			diff := D.At(i, j) - delta
			frobNorm += diff * diff
		}
	}
	frobNorm = math.Sqrt(frobNorm)
	// fmt.Printf("Frobenius norm ||D - I||_F = %.10e\n", frobNorm)
	if frobNorm > tol {
		panic("The basis is not orthonormalized properly")
	}
}

func (el *LagrangeElement2D) GetDerivativeMatrices(R, S utils.Vector) (Dr, Ds utils.Matrix) {
	Vr, Vs := el.JB2D.GradVandermonde2D(el.N, R, S)
	Dr, Ds = Vr.Mul(el.JB2D.Vinv), Vs.Mul(el.JB2D.Vinv)
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
			R: cubMat.Col(0),
			S: cubMat.Col(1),
			W: cubMat.Col(2),
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

		   cub.R = cubR
		   cub.S = cubS
		   cub.W = cubW
		   cub.Ncub = cub.R.size()
		*/
	}
	return
}
