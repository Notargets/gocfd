package DG1D

import (
	"fmt"
	"math"

	"github.com/james-bowman/sparse"
	"github.com/notargets/gocfd/utils"
)

func (el *Elements1D) Startup1D() {
	var (
		err  error
		Vinv utils.Matrix
		N    = el.Np - 1
	)
	R := JacobiGL(0, 0, N)
	V := Vandermonde1D(N, R)
	if Vinv, err = V.Inverse(); err != nil {
		fmt.Println(err)
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
	el.Connect1D()
	el.BuildMaps1D()
	return
}

func (el *Elements1D) Connect1D() {
	var (
		NFaces     = 2
		K, _       = el.EToV.Dims()
		Nv         = K + 1
		TotalFaces = NFaces * K
		vn         = utils.NewVector(2, []float64{0, 1}) // local face to vertex connections
	)
	SpFToV_Tmp := sparse.NewDOK(TotalFaces, Nv)
	var sk int
	for k := 0; k < K; k++ {
		for face := 0; face < NFaces; face++ {
			col := int(vn.AtVec(face))
			SpFToV_Tmp.Set(sk, int(el.EToV.At(k, col)), 1)
			sk++
		}
	}
	SpFToF := sparse.NewCSR(TotalFaces, TotalFaces, nil, nil, nil)
	SpFToV := SpFToV_Tmp.ToCSR()
	SpFToF.Mul(SpFToV, SpFToV.T())
	for i := 0; i < TotalFaces; i++ {
		v := SpFToF.At(i, i)
		SpFToF.Set(i, i, v-2)
	}
	FacesIndex := utils.MatFind(SpFToF, utils.Equal, 1)

	element1 := FacesIndex.RI.Copy().Apply(func(val int) int { return val / NFaces })
	face1 := FacesIndex.RI.Copy().Apply(func(val int) int { return int(math.Mod(float64(val), float64(NFaces))) })

	element2 := FacesIndex.CI.Copy().Apply(func(val int) int { return val / NFaces })
	face2 := FacesIndex.CI.Copy().Apply(func(val int) int { return int(math.Mod(float64(val), float64(NFaces))) })

	// Rearrange into Nelements x Nfaces sized arrays
	el.EToE = utils.NewRangeOffset(1, K).Outer(utils.NewOnes(NFaces))
	el.EToF = utils.NewOnes(K).Outer(utils.NewRangeOffset(1, NFaces))
	var I2D utils.Index2D
	var err error
	nr, nc := el.EToE.Dims()
	if I2D, err = utils.NewIndex2D(nr, nc, element1, face1); err != nil {
		panic(err)
	}
	if err = el.EToE.IndexedAssign(I2D, element2); err != nil {
		panic(err)
	}
	if err = el.EToF.IndexedAssign(I2D, face2); err != nil {
		panic(err)
	}
	return
}

func (el *Elements1D) BuildMaps1D() {
	var (
		NF = el.Nfp * el.NFaces
	)
	// number volume nodes consecutively
	nodeids := utils.NewRangeOffset(1, el.Np*el.K)

	// find index of face nodes with respect to volume node ordering
	el.VmapM = utils.NewIndex(el.Nfp * el.NFaces * el.K)
	idsR := utils.NewFromFloat(el.FMask.RawVector().Data)
	for k1 := 0; k1 < el.K; k1++ {
		iL1 := k1 * NF
		iL2 := iL1 + NF
		idsL := utils.NewRangeOffset(iL1+1, iL2) // sequential indices for element k1
		if err := el.VmapM.IndexedAssign(idsL, nodeids.Subset(idsR)); err != nil {
			panic(err)
		}
		idsR.Add(el.Np)
	}

	//var one = utils.NewVecConst(Nfp, 1)
	var one = utils.NewVector(el.Nfp).Set(1)
	el.VmapP = utils.NewIndex(el.Nfp * el.NFaces * el.K)
	//fmt.Printf("X = \n%v\n", mat.Formatted(X, mat.Squeeze()))
	for k1 := 0; k1 < el.K; k1++ {
		for f1 := 0; f1 < el.NFaces; f1++ {
			k2 := int(el.EToE.At(k1, f1))
			f2 := int(el.EToF.At(k1, f1))
			skM := k1 * NF
			skP := k2 * NF
			idsM := utils.NewRangeOffset(1+f1*el.Nfp+skM, (f1+1)*el.Nfp+skM)
			idsP := utils.NewRangeOffset(1+f2*el.Nfp+skP, (f2+1)*el.Nfp+skP)
			vidM := el.VmapM.Subset(idsM)
			vidP := el.VmapM.Subset(idsP)
			x1 := el.X.SubsetVector(vidM)
			x2 := el.X.SubsetVector(vidP)
			X1 := x1.Outer(one)
			X2 := x2.Outer(one)
			D := X1.Copy().Subtract(X2.Transpose()).POW(2).Apply(math.Sqrt).Apply(math.Abs)
			v1 := int(el.EToV.At(k1, f1))
			v2 := int(el.EToV.At(k1, (f1+1)%el.NFaces))
			refd := math.Sqrt(utils.POW(el.VX.AtVec(v1)-el.VX.AtVec(v2), 2))
			idMP := D.Find(utils.Less, utils.NODETOL*refd)
			idM := idMP.RI
			idP := idMP.CI
			if err := el.VmapP.IndexedAssign(idM.Copy().Add(f1*el.Nfp+skM), vidP.Subset(idP)); err != nil {
				panic(err)
			}
		}
	}
	// Create list of boundary nodes
	el.MapB = el.VmapP.Compare(utils.Equal, el.VmapM)
	el.VmapB = el.VmapM.Subset(el.MapB)
	el.MapI = utils.NewIndex(1)
	el.MapO = utils.NewIndex(1).Copy().Add(el.K*el.NFaces - 1)
	el.VmapI = utils.NewIndex(1)
	el.VmapO = utils.NewIndex(1).Copy().Add(el.K*el.Np - 1)
	return
}
