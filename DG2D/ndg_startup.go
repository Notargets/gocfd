package DG2D

import "C"
import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

type NDG2D struct {
	Element                           *LagrangeElement2D
	K                                 int
	NODETOL                           float64
	VX, VY, VZ                        utils.Vector
	FMask, Fx, Fy                     utils.Matrix
	EToV, EToE, EToF                  utils.Matrix
	BCType                            utils.Matrix
	X, Y, FScale, LIFT                utils.Matrix
	Rx, Ry, Sx, Sy                    utils.Matrix
	xs, xr, ys, yr                    utils.Matrix
	NX, NY                            utils.Matrix
	J, sJ                             utils.Matrix
	VmapM, VmapP, VmapB, VmapI, VmapO utils.Index
	MapB, MapI, MapO                  utils.Index
	FMaskI                            utils.Index
}

func NewNDG2D(N int, meshFile string, plotMesh bool) (ndg *NDG2D) {
	ndg = &NDG2D{
		Element: NewLagrangeElement2D(N, Hesthaven),
		NODETOL: 1.e-6,
	}
	if N < 1 {
		panic(fmt.Errorf("Polynomial order must be >= 1, have %d", N))
	}
	ndg.K, ndg.VX, ndg.VY, ndg.EToV, ndg.BCType = ReadGambit2d(meshFile)
	ndg.Startup2D()
	if plotMesh {
		PlotMesh(ndg.VX, ndg.VY, ndg.EToV, ndg.BCType, ndg.X, ndg.Y, true)
		utils.SleepFor(50000)
	}
	return
}

func (ndg *NDG2D) Startup2D() {
	var (
		R, S = ndg.Element.R, ndg.Element.S
	)
	// build coordinates of all the nodes
	va, vb, vc := ndg.EToV.Col(0), ndg.EToV.Col(1), ndg.EToV.Col(2)
	ndg.X = R.Copy().Add(S).Scale(-1).Outer(ndg.VX.SubsetIndex(va.ToIndex())).Add(
		R.Copy().AddScalar(1).Outer(ndg.VX.SubsetIndex(vb.ToIndex()))).Add(
		S.Copy().AddScalar(1).Outer(ndg.VX.SubsetIndex(vc.ToIndex()))).Scale(0.5)
	ndg.Y = R.Copy().Add(S).Scale(-1).Outer(ndg.VY.SubsetIndex(va.ToIndex())).Add(
		R.Copy().AddScalar(1).Outer(ndg.VY.SubsetIndex(vb.ToIndex()))).Add(
		S.Copy().AddScalar(1).Outer(ndg.VY.SubsetIndex(vc.ToIndex()))).Scale(0.5)
	fmask1 := S.Copy().AddScalar(1).Find(utils.Less, ndg.NODETOL, true)
	fmask2 := S.Copy().Add(R).Find(utils.Less, ndg.NODETOL, true)
	fmask3 := R.Copy().AddScalar(1).Find(utils.Less, ndg.NODETOL, true)
	if fmask1.Len() != 0 {
		ndg.FMask = utils.NewMatrix(ndg.Element.Nfp, 3)
		ndg.FMask.SetCol(0, fmask1.Data())
		ndg.FMask.SetCol(1, fmask2.Data())
		ndg.FMask.SetCol(2, fmask3.Data())
		ndg.FMaskI = utils.NewIndex(len(ndg.FMask.Data()), ndg.FMask.Data())
		ndg.Fx = utils.NewMatrix(3*ndg.Element.Nfp, ndg.K)
		for fp, val := range ndg.FMask.Data() {
			ind := int(val)
			ndg.Fx.M.SetRow(fp, ndg.X.M.RawRowView(ind))
		}
		ndg.Fy = utils.NewMatrix(3*ndg.Element.Nfp, ndg.K)
		for fp, val := range ndg.FMask.Data() {
			ind := int(val)
			ndg.Fy.M.SetRow(fp, ndg.Y.M.RawRowView(ind))
		}
		ndg.Lift2D()
	}
	ndg.GeometricFactors2D()
	ndg.Normals2D()
	ndg.FScale = ndg.sJ.ElDiv(ndg.J.Subset(ndg.GetFaces()))
	// Build connectivity matrices
	ndg.EToE, ndg.EToF = Connect2D(ndg.K, ndg.Element.NFaces, ndg.VX.Len(), ndg.EToV)

	// Mark fields read only
	ndg.LIFT.SetReadOnly("LIFT")
	ndg.X.SetReadOnly("X")
	ndg.Y.SetReadOnly("Y")
	ndg.Fx.SetReadOnly("Fx")
	ndg.Fy.SetReadOnly("Fy")
	ndg.FMask.SetReadOnly("FMask")
	ndg.NX.SetReadOnly("NX")
	ndg.NY.SetReadOnly("NY")
	ndg.FScale.SetReadOnly("FScale")
	ndg.EToE.SetReadOnly("EToE")
	ndg.EToF.SetReadOnly("EToF")
	return
}

func (ndg *NDG2D) GetFaces() (aI utils.Index, NFacePts, K int) {
	var (
		err      error
		allFaces utils.Index2D
	)
	NFacePts = ndg.Element.Nfp * ndg.Element.NFaces
	K = ndg.K
	allK := utils.NewRangeOffset(1, ndg.K)
	if allFaces, err = utils.NewIndex2D(ndg.Element.Np, ndg.K, ndg.FMaskI, allK, true); err != nil {
		panic(err)
	}
	aI = allFaces.ToIndex()
	return
}

func (ndg *NDG2D) Normals2D() {
	var (
		f1, f2, f3 utils.Index2D
		err        error
	)
	allK := utils.NewRangeOffset(1, ndg.K)
	aI, NFacePts, _ := ndg.GetFaces()
	// interpolate geometric factors to face nodes
	fxr := ndg.xr.Subset(aI, NFacePts, ndg.K)
	fxs := ndg.xs.Subset(aI, NFacePts, ndg.K)
	fyr := ndg.yr.Subset(aI, NFacePts, ndg.K)
	fys := ndg.ys.Subset(aI, NFacePts, ndg.K)
	// build normals
	faces1 := utils.NewRangeOffset(1, ndg.Element.Nfp)
	faces2 := utils.NewRangeOffset(1+ndg.Element.Nfp, 2*ndg.Element.Nfp)
	faces3 := utils.NewRangeOffset(1+2*ndg.Element.Nfp, 3*ndg.Element.Nfp)
	if f1, err = utils.NewIndex2D(NFacePts, ndg.K, faces1, allK, true); err != nil {
		panic(err)
	}
	if f2, err = utils.NewIndex2D(NFacePts, ndg.K, faces2, allK, true); err != nil {
		panic(err)
	}
	if f3, err = utils.NewIndex2D(NFacePts, ndg.K, faces3, allK, true); err != nil {
		panic(err)
	}
	ndg.NX, ndg.NY = utils.NewMatrix(NFacePts, ndg.K), utils.NewMatrix(NFacePts, ndg.K)
	// Face 1
	ndg.NX.Assign(f1.ToIndex(), fyr.Subset(f1.ToIndex(), f1.Len, 1))
	ndg.NY.Assign(f1.ToIndex(), fxr.Subset(f1.ToIndex(), f1.Len, 1).Scale(-1))
	// Face 2
	ndg.NX.Assign(f2.ToIndex(), fys.Subset(f2.ToIndex(), f2.Len, 1).Subtract(fyr.Subset(f2.ToIndex(), f2.Len, 1)))
	ndg.NY.Assign(f2.ToIndex(), fxs.Subset(f2.ToIndex(), f2.Len, 1).Scale(-1).Add(fxr.Subset(f2.ToIndex(), f2.Len, 1)))
	// Face 3
	ndg.NX.Assign(f3.ToIndex(), fys.Subset(f3.ToIndex(), f3.Len, 1).Scale(-1))
	ndg.NY.Assign(f3.ToIndex(), fxs.Subset(f3.ToIndex(), f3.Len, 1))
	ndg.sJ = ndg.NX.Copy().POW(2).Add(ndg.NY.Copy().POW(2)).Apply(func(val float64) (res float64) {
		res = math.Sqrt(val)
		return
	})
	ndg.NX.ElDiv(ndg.sJ)
	ndg.NY.ElDiv(ndg.sJ)
}

func (ndg *NDG2D) GeometricFactors2D() {
	/*
	  // function [rx,sx,ry,sy,J] = GeometricFactors2D(x,y,Dr,Ds)
	  // Purpose  : Compute the metric elements for the local
	  //            mappings of the elements
	  DMat xr=Dr*x
	       xs=Ds*x
	       yr=Dr*y
	       ys=Ds*y;
	  J  =  xr.dm(ys) - xs.dm(yr);
	  rx = ys.dd(J); sx = -yr.dd(J); ry = -xs.dd(J); sy = xr.dd(J);
	*/
	// Calculate geometric factors
	ndg.xr, ndg.xs = ndg.Element.Dr.Mul(ndg.X), ndg.Element.Ds.Mul(ndg.X)
	ndg.yr, ndg.ys = ndg.Element.Dr.Mul(ndg.Y), ndg.Element.Ds.Mul(ndg.Y)
	ndg.xr.SetReadOnly("xr")
	ndg.xs.SetReadOnly("xs")
	ndg.yr.SetReadOnly("yr")
	ndg.ys.SetReadOnly("ys")
	ndg.J = ndg.xr.Copy().ElMul(ndg.ys).Subtract(ndg.xs.Copy().ElMul(ndg.yr))
	ndg.Rx = ndg.ys.Copy().ElDiv(ndg.J)
	ndg.Sx = ndg.yr.Copy().ElDiv(ndg.J).Scale(-1)
	ndg.Ry = ndg.xs.Copy().ElDiv(ndg.J).Scale(-1)
	ndg.Sy = ndg.xr.Copy().ElDiv(ndg.J)
}

func (ndg *NDG2D) Lift2D() {
	var (
		err      error
		I2       utils.Index2D
		massEdge utils.Matrix
		V1D      utils.Matrix
		Emat     utils.Matrix
		R, S     = ndg.Element.R, ndg.Element.S
	)
	Emat = utils.NewMatrix(ndg.Element.Np, ndg.Element.NFaces*ndg.Element.Nfp)
	faceMap := func(basis utils.Vector, faceNum int, Ind utils.Index) {
		faceBasis := basis.SubsetIndex(ndg.FMask.Col(faceNum).ToIndex())
		V1D = DG1D.Vandermonde1D(ndg.Element.N, faceBasis)
		if massEdge, err = V1D.Mul(V1D.Transpose()).Inverse(); err != nil {
			panic(err)
		}
		if I2, err = utils.NewIndex2D(ndg.Element.Np, ndg.Element.NFaces*ndg.Element.Nfp, ndg.FMask.Col(faceNum).ToIndex(), Ind, true); err != nil {
			panic(err)
		}
		Emat.Assign(I2.ToIndex(), massEdge)
	}
	// face 1
	faceMap(R, 0, utils.NewRangeOffset(1, ndg.Element.Nfp))
	// face 2
	faceMap(R, 1, utils.NewRangeOffset(ndg.Element.Nfp+1, 2*ndg.Element.Nfp))
	// face 3
	faceMap(S, 2, utils.NewRangeOffset(2*ndg.Element.Nfp+1, 3*ndg.Element.Nfp))
	// inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
	ndg.LIFT = ndg.Element.V.Mul(ndg.Element.V.Transpose().Mul(Emat))
	return
}

func (ndg *NDG2D) BuildMaps2D() {
	return
}

func (ndg *NDG2D) NewCube2D(COrder int) {
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
		ndg.Element.Cub = &Cubature{
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
