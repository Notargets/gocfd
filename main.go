package main

import (
	"github.com/notargets/gofluids/DG1D"
	"github.com/notargets/gofluids/utils"
	"gonum.org/v1/gonum/mat"
)

const (
	NODETOL = 1.e-12
	K       = 10
	N       = 8
	Nfp     = 1
	NFaces  = 2
	Np      = N + 1
)

func main() {
	Startup1D()
}

func SimpleMesh1D(xmin, xmax float64, K int) (VX *mat.VecDense, EToV *mat.Dense) {
	// K is the number of elements, there are K+1 vertices
	var (
		x             = make([]float64, K+1)
		elementVertex = make([]float64, K*2)
	)
	for i := 0; i < K+1; i++ {
		x[i] = (xmax-xmin)*float64(i)/float64(K) + xmin
	}
	var iter int
	for i := 0; i < K; i++ {
		elementVertex[iter] = float64(i)
		elementVertex[iter+1] = float64(i + 1)
		iter += 2
	}
	return mat.NewVecDense(K+1, x), mat.NewDense(K, 2, elementVertex)
}

func Startup1D() {
	VX, EToV := SimpleMesh1D(0, 2, K)

	_, R, W := DG1D.JacobiGL(0, 0, N)
	V := DG1D.Vandermonde1D(N, R)
	Vinv := mat.NewDense(Np, Np, nil)
	if err := Vinv.Inverse(V); err != nil {
		panic("error inverting V")
	}
	Vr := DG1D.GradVandermonde1D(R, N)
	Dr := mat.NewDense(Np, Np, nil)
	Dr.Product(Vr, Vinv)
	LIFT := DG1D.Lift1D(V, Np, NFaces, Nfp)

	NX := DG1D.Normals1D(NFaces, Nfp, K)

	//fmt.Printf("LIFT = \n%v\n", mat.Formatted(LIFT, mat.Squeeze()))
	va := EToV.ColView(0)
	vb := EToV.ColView(1)
	sT := mat.NewVecDense(va.Len(), nil)
	sT.SubVec(utils.VecSub(VX, vb), utils.VecSub(VX, va))

	// x = ones(Np)*VX(va) + 0.5*(r+1.)*sT(vc);
	ones := utils.VecConst(1, Np)
	mm := mat.NewDense(Np, K, nil)
	mm.Mul(ones, utils.VecSub(VX, va).T())

	rr := utils.VecScalarAdd(mat.VecDenseCopyOf(R), 1)
	rr.ScaleVec(0.5, rr)
	X := mat.NewDense(Np, K, nil)
	X.Mul(rr, sT.T())
	X.Add(X, mm)

	J, Rx := DG1D.GeometricFactors1D(Dr, X)

	fmask1 := utils.VecFind(utils.VecScalarAdd(R, 1), utils.Less, NODETOL, true)
	fmask2 := utils.VecFind(utils.VecScalarAdd(R, -1), utils.Less, NODETOL, true)
	FMask := utils.VecConcat(fmask1, fmask2)
	Fx := utils.MatSubRow(X, FMask)
	JJ := utils.MatSubRow(J, FMask)
	FScale := utils.MatElementInvert(JJ)

	EToE, EToF := DG1D.Connect1D(EToV)

	mapM, mapP, vmapM, vmapP, mapB, vmapB :=
		DG1D.BuildMaps1D(VX, FMask,
			EToV, EToE, EToF,
			K, Np, Nfp, NFaces,
			NODETOL)
	_, _, _, _, _, _, _, _, _, _, _, _, _ = VX, EToV, J, W, LIFT, NX, X, JJ, Rx, Fx, FScale, EToE, EToF
	_, _, _, _, _, _ = mapM, mapP, vmapM, vmapP, mapB, vmapB
}
