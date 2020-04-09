package DG1D

import "github.com/notargets/gocfd/utils"

func (el *Elements1D) Startup1D() {
	var (
		err  error
		Vinv utils.Matrix
		N    = el.Np - 1
	)
	R := JacobiGL(0, 0, N)
	V := Vandermonde1D(N, R)
	if Vinv, err = V.Inverse(); err != nil {
		panic("error inverting V")
	}
	Vr := GradVandermonde1D(R, N)

	el.Dr = Vr.Mul(Vinv)

	el.LIFT = Lift1D(V, el.Np, el.NFaces, el.Nfp)

	el.NX = Normals1D(el.NFaces, el.Nfp, el.K)

	va := el.EToV.Col(0).ToIndex()
	vb := el.EToV.Col(1).ToIndex()
	sT := el.VX.Subset(vb).Subtract(el.VX.Subset(va))
	// x = ones(Np)*VX(va) + 0.5*(r+1.)*sT(vc);
	mm := utils.NewVector(el.Np).Set(1).Mul(el.VX.Subset(va))
	el.X = R.Copy().AddScalar(1).Scale(0.5).Mul(sT).Add(mm)

	var J utils.Matrix
	J, el.Rx = GeometricFactors1D(el.Dr, el.X)

	fmask1 := R.Copy().AddScalar(1).Find(utils.Less, utils.NODETOL, true)
	fmask2 := R.Copy().AddScalar(-1).Find(utils.Less, utils.NODETOL, true)
	el.FMask = fmask1.Concat(fmask2)
	el.FScale = J.SliceRows(el.FMask.ToIndex()).POW(-1)

	el.EToE, el.EToF = Connect1D(el.EToV)

	el.vmapM, el.vmapP, el.mapB, el.vmapB, el.mapI, el.vmapI, el.mapO, el.vmapO =
		BuildMaps1D(el.VX, el.FMask,
			el.X, el.EToV, el.EToE, el.EToF,
			el.K, el.Np, el.Nfp, el.NFaces)
	return
}
