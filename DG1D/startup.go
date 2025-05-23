package DG1D

import (
	"fmt"
	"math"

	"github.com/james-bowman/sparse"
	"github.com/notargets/gocfd/utils"
)

func (el *Elements1D) Startup1D(nt NODE_TYPE) {
	var (
		err error
		N   = el.Np - 1
	)
	switch nt {
	case GAUSS:
		el.R, _ = JacobiGQ(1, 1, N)
	case GAUSS_LOBATO:
		el.R = JacobiGL(0, 0, N)
	}
	el.V = Vandermonde1D(N, el.R)
	if el.Vinv, err = el.V.Inverse(); err != nil {
		fmt.Println(err)
		panic("error inverting V")
	}
	Vr := GradVandermonde1D(N, el.R)

	el.Dr = Vr.Mul(el.Vinv)

	el.LIFT = Lift1D(el.V, el.Np, el.NFaces, el.Nfp)

	el.NX = Normals1D(el.NFaces, el.Nfp, el.K)

	va := el.EToV.Col(0).ToIndex()
	vb := el.EToV.Col(1).ToIndex()
	sT := el.VX.SubsetIndex(vb).Subtract(el.VX.SubsetIndex(va))
	// x = ones(Np)*VX(va) + 0.5*(r+1.)*sT(vc);
	mm := utils.NewVector(el.Np).SetScalar(1).Mul(el.VX.SubsetIndex(va))
	el.X = el.R.Copy().AddScalar(1).Scale(0.5).Mul(sT).Add(mm)

	var J utils.Matrix
	J, el.Rx = GeometricFactors1D(el.Dr, el.X)

	if nt == GAUSS_LOBATO { // We need FScale for Galerkin on GL nodes, not for Gauss with DFR
		fmask1 := el.R.Copy().AddScalar(1).Find(utils.Less, utils.NODETOL, true)
		fmask2 := el.R.Copy().AddScalar(-1).Find(utils.Less, utils.NODETOL, true)
		el.FMask = fmask1.Concat(fmask2)
		el.FScale = J.SliceRows(el.FMask.ToIndex()).POW(-1)
	}
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
		vn         = []int{0, 1} // local face to vertex connections
	)
	/*
		EToE: We will find the neighbor elements for each element
		Example: K = 4, each element has a left and right neighbor:
		EToE =
		⎡0.0000  1.0000⎤	Boundary element, left node is itself, right node is 1
		⎢0.0000  2.0000⎥	Node 1, internal, left is 0, right is 2
		⎢1.0000  3.0000⎥
		⎣2.0000  3.0000⎦

		EToF: We will find the neighbor element face number for each element
		EToF =
		⎡0.0000  0.0000⎤    Boundary element, face number is arbitrary and the same for R and L
		⎢1.0000  0.0000⎥	Left neighbor's face number 1, right neighbor's face number 0
		⎢1.0000  0.0000⎥
		⎣1.0000  1.0000⎦
	*/
	SpFToV_Tmp := sparse.NewDOK(TotalFaces, Nv)
	var sk int
	for k := 0; k < K; k++ {
		for face := 0; face < NFaces; face++ {
			SpFToV_Tmp.Set(sk, int(el.EToV.At(k, vn[face])), 1)
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
	// Find where faces share a vertex with another face
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
	el.EToE.Assign(I2D.ToIndex(), element2)
	el.EToF.Assign(I2D.ToIndex(), face2)
	return
}

func (el *Elements1D) BuildMaps1D() {
	el.VmapM, el.VmapP = el.FaceMap()
	// Create list of boundary nodes
	el.MapB = el.VmapP.Compare(utils.Equal, el.VmapM)
	el.VmapB = el.VmapM.Subset(el.MapB)
	el.MapI = utils.NewIndex(1)
	el.MapO = utils.NewIndex(1).Add(el.K*el.NFaces - 1)
	el.VmapI = utils.NewIndex(1)
	el.VmapO = utils.NewIndex(1).Add(el.K*el.Np - 1)
	el.VmapI = utils.NewIndex(1)
	el.VmapO = utils.NewIndex(1).Add(el.K*el.Np - 1)
	return
}

func (el *Elements1D) FaceMap() (VmapM, VmapP utils.Index) {
	/*
				We need a map of left / right face points for each element
		    	VmapM(NFaces, K) and VmapP(NFaces, K)
				The mapping needs to be in column-major form to implement an (NFaces, K) matrix
	*/
	var (
		RangerNpK = utils.NewR2(el.Np, el.K)
	)
	// Left and right ends for each element
	indLeft := RangerNpK.Range(0, ":")
	indRight := RangerNpK.Range(-1, ":")
	VmapM = make(utils.Index, 2*el.K)
	var ind int
	for i, val := range indLeft {
		VmapM[ind] = val
		VmapM[ind+el.K] = indRight[i]
		ind++
	}
	VmapP = make(utils.Index, 2*el.K)
	ind = 0
	for k := 0; k < el.K; k++ {
		nLeft := el.Np - 1
		nRight := 0
		if k == 0 {
			nLeft = 0
		}
		if k == el.K-1 {
			nRight = el.Np - 1
		}
		kLeft := k - 1
		kRight := k + 1
		kLeft = int(math.Max(0, float64(kLeft)))
		kRight = int(math.Min(float64(el.K-1), float64(kRight)))
		/*
			fmt.Printf("nL, kL, nR, kR = %d, %d, %d, %d\n", nLeft, kLeft, nRight, kRight)
			fmt.Printf("Range Left = %v\n", RangerNpK.Range(nLeft, kLeft))
			fmt.Printf("Range Right = %v\n", RangerNpK.Range(nRight, kRight))
		*/
		VmapP[ind] = RangerNpK.Range(nLeft, kLeft)[0]
		VmapP[ind+el.K] = RangerNpK.Range(nRight, kRight)[0]
		ind++
	}
	return
}
