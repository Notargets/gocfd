package DG1D

import (
	"github.com/notargets/gocfd/utils"
	"gonum.org/v1/gonum/mat"
)

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

func Startup1D(K, N, NFaces, Nfp int) (X *mat.Dense) {
	var (
		Np = N + 1
	)
	VX, EToV := SimpleMesh1D(0, 2, K)

	_, R, W := JacobiGL(0, 0, N)
	V := Vandermonde1D(N, R)
	Vinv := mat.NewDense(Np, Np, nil)
	if err := Vinv.Inverse(V); err != nil {
		panic("error inverting V")
	}
	Vr := GradVandermonde1D(R, N)
	Dr := mat.NewDense(Np, Np, nil)
	Dr.Product(Vr, Vinv)
	LIFT := Lift1D(V, Np, NFaces, Nfp)

	NX := Normals1D(NFaces, Nfp, K)

	va := EToV.Col(0).ToIndex()
	vb := EToV.Col(1).ToIndex()
	sT := VX.Subset(vb).Subtract(VX.Subset(va))

	// x = ones(Np)*VX(va) + 0.5*(r+1.)*sT(vc);
	mm := utils.NewVector(Np).Set(1).Mul(VX.Subset(va))
	r := utils.Vector{mat.VecDenseCopyOf(R)}
	X = r.AddScalar(1).Scale(0.5).Mul(sT).Add(mm).M

	J, Rx := GeometricFactors1D(Dr, X)

	r = utils.Vector{R}
	fmask1 := r.Copy().AddScalar(1).Find(utils.Less, utils.NODETOL, true)
	fmask2 := r.Copy().AddScalar(-1).Find(utils.Less, utils.NODETOL, true)
	FMask := fmask1.Concat(fmask2)
	Xm := utils.Matrix{X}
	Fx := Xm.SliceRows(FMask.ToIndex())
	Jm := utils.Matrix{J}
	JJ := Jm.SliceRows(FMask.ToIndex())
	FScale := JJ.POW(-1)

	EToE, EToF := Connect1D(EToV.M)

	vmapM, vmapP, mapB, vmapB, mapI, vmapI, mapO, vmapO :=
		BuildMaps1D(VX.V, FMask.V,
			X, EToV.M, EToE, EToF,
			K, Np, Nfp, NFaces,
			utils.NODETOL)
	_, _, _, _, _, _ = W, LIFT, NX, Rx, Fx, FScale
	_, _, _, _, _, _, _, _ = vmapM, vmapP, mapB, vmapB, mapI, vmapI, mapO, vmapO
	return
}
