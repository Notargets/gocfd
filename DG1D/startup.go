package DG1D

import "github.com/notargets/gocfd/utils"

func SimpleMesh1D(xmin, xmax float64, K int) (VX utils.Vector, EToV utils.Matrix) {
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
	return utils.NewVector(K+1, x), utils.NewMatrix(K, 2, elementVertex)
}

func Startup1D(K, N, NFaces, Nfp int) (X utils.Matrix) {
	var (
		Np   = N + 1
		err  error
		Vinv utils.Matrix
	)
	VX, EToV := SimpleMesh1D(0, 2, K)

	R, W := JacobiGL(0, 0, N)
	V := Vandermonde1D(N, R)
	if Vinv, err = V.Inverse(); err != nil {
		panic("error inverting V")
	}
	Vr := GradVandermonde1D(R, N)

	Dr := Vr.Mul(Vinv)

	LIFT := Lift1D(V, Np, NFaces, Nfp)

	NX := Normals1D(NFaces, Nfp, K)

	va := EToV.Col(0).ToIndex()
	vb := EToV.Col(1).ToIndex()
	sT := VX.Subset(vb).Subtract(VX.Subset(va))
	// x = ones(Np)*VX(va) + 0.5*(r+1.)*sT(vc);
	mm := utils.NewVector(Np).Set(1).Mul(VX.Subset(va))
	X = R.Copy().AddScalar(1).Scale(0.5).Mul(sT).Add(mm)
	//fmt.Printf("VX = \n%v\n", mat.Formatted(VX, mat.Squeeze()))
	//fmt.Printf("X = \n%v\n", mat.Formatted(X, mat.Squeeze()))

	J, Rx := GeometricFactors1D(Dr, X)

	fmask1 := R.Copy().AddScalar(1).Find(utils.Less, utils.NODETOL, true)
	fmask2 := R.Copy().AddScalar(-1).Find(utils.Less, utils.NODETOL, true)
	FMask := fmask1.Concat(fmask2)
	Fx := X.SliceRows(FMask.ToIndex())
	FScale := J.SliceRows(FMask.ToIndex()).POW(-1)

	EToE, EToF := Connect1D(EToV)

	vmapM, vmapP, mapB, vmapB, mapI, vmapI, mapO, vmapO :=
		BuildMaps1D(VX, FMask,
			X, EToV, EToE, EToF,
			K, Np, Nfp, NFaces,
			utils.NODETOL)
	_, _, _, _, _, _ = W, LIFT, NX, Rx, Fx, FScale
	_, _, _, _, _, _, _, _ = vmapM, vmapP, mapB, vmapB, mapI, vmapI, mapO, vmapO
	return
}
