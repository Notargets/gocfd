package DG2D

import "C"
import (
	"fmt"
	"image/color"
	"math"

	"github.com/notargets/avs/chart2d"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

type NDG_Elements2D struct {
	K, N, Nfp, Np, NFaces             int
	NODETOL                           float64
	R, S, VX, VY, VZ                  utils.Vector
	FMask, Fx, Fy                     utils.Matrix
	EToV, EToE, EToF                  utils.Matrix
	BCType                            utils.Matrix
	X, Y, Dr, Ds, FScale, LIFT        utils.Matrix
	Rx, Ry, Sx, Sy                    utils.Matrix
	xs, xr, ys, yr                    utils.Matrix
	V, Vinv, MassMatrix               utils.Matrix
	NX, NY                            utils.Matrix
	J, sJ                             utils.Matrix
	VmapM, VmapP, VmapB, VmapI, VmapO utils.Index
	MapB, MapI, MapO                  utils.Index
	Cub                               *Cubature
	FMaskI                            utils.Index
	RT                                *RTElement
}

type Cubature struct {
	r, s, w                 utils.Vector
	W                       utils.Matrix
	V, Dr, Ds, VT, DrT, DsT utils.Matrix
	x, y, rx, sx, ry, sy, J utils.Matrix
	mm, mmCHOL              utils.Matrix
}

func NewNDG2D(N int, meshFile string, plotMesh bool) (el *NDG_Elements2D) {
	el = &NDG_Elements2D{
		N:      N,
		Np:     (N + 1) * (N + 2) / 2,
		NFaces: 3,
	}
	if N < 1 {
		panic(fmt.Errorf("Polynomial order must be >= 1, have %d", N))
	}
	el.ReadGambit2d(meshFile, plotMesh)
	el.Startup2D()
	el.Startup2DDFR()
	if plotMesh {
		var (
			chart *chart2d.Chart2D
		)
		white := color.RGBA{
			R: 255,
			G: 255,
			B: 255,
			A: 0,
		}
		blue := color.RGBA{
			R: 50,
			G: 0,
			B: 255,
			A: 0,
		}
		red := color.RGBA{
			R: 255,
			G: 0,
			B: 50,
			A: 0,
		}
		_, _, _ = white, red, blue
		chart = PlotMesh(el.VX, el.VY, el.EToV, el.BCType, el.X, el.Y, true)
		_ = chart
		utils.SleepFor(50000)
	}
	return
}

func (el *NDG_Elements2D) Simplex2DInterpolate(r, s float64, f []float64) (value float64) {
	var (
		N  = el.N
		Np = el.Np
	)
	if len(f) != Np {
		panic(fmt.Errorf("not enough function values for the 2D basis, have %d, need %d", len(f), Np))
	}
	// First compute polynomial terms, used by all polynomials
	polyTerms := make([]float64, Np)
	var sk int
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			polyTerms[sk] = Simplex2DPTerm(r, s, i, j)
			sk++
		}
	}
	ptV := utils.NewMatrix(Np, 1, polyTerms)
	pAtR := el.Vinv.Transpose().Mul(ptV)
	fV := utils.NewVector(Np, f)
	value = pAtR.Col(0).Dot(fV)
	return
}

func (el *NDG_Elements2D) Simplex2DInterpolatingPolyMatrix(R, S utils.Vector) (polyMatrix utils.Matrix) {
	/*
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
	for ii, r := range R.Data() {
		s := S.Data()[ii]
		for i := 0; i <= N; i++ {
			for j := 0; j <= (N - i); j++ {
				polyTerms[sk] = Simplex2DPTerm(r, s, i, j)
				sk++
			}
		}
	}
	ptV := utils.NewMatrix(Np, R.Len(), polyTerms).Transpose()
	polyMatrix = el.Vinv.Transpose().Mul(ptV)
	return
}

func (el *NDG_Elements2D) Startup2DDFR() {
	// Build reference element matrices
	/*
			We build the mixed elements for the DFR scheme with:

			Solution Points: We use points within a reference triangle, excluding the edges, for a Lagrangian element
			of O(K) to store the solution. If we need derivatives, or interpolated quantities (Flux), we use the
			solution points element.

			Flux Points: We use a customized Raviart-Thomas (RT) vector element of O(K+1) to store the vector Flux function
		    computed from the solution values. The RT element is of order O(K+1) and is a combination of the points from
			the solution element for the interior, and points along the three triangle edges. The custom RT basis is
			established using a procedure outlined in: "Ainv Direct Flux Reconstruction Scheme for Advection-Diffusion
			Problems on Triangular Grids" by Romero, Witherden and Jameson (2017). Ainv complete RT basis, [ B ], is used
			together with unit basis vectors, [ w ], to satisfy the following:
					[ B_j(r_i) dot w_i ] [ C ] = [ delta_i_j ]
					=> solve for [ C ], the coefficients defining the custom RT basis

			[ C ] is the vector of coefficients defining the basis using the basis vectors [ w ] and [ B ].

			The [ w ] directions of the custom RT element basis are defined such that:
				w([r]) = w(edge_locations) = unit normals on each of three edges
				w([r]) = w(interior) = unit normals in the two primary geometry2D directions (r and s)

			For order K there are:
				- (K+1) locations on each edge, for a total of 3(K+1) edge basis functions.
				- (K)(K+1) locations in the interior, half for the w_r direction and half for the w_s direction
				- Total: (K+3)(K+1) basis functions for the custom RT_K element

			Notes:
				1) The number of interior points matches the Lagrangian element in 2D at order (K-1). Ainv Lagrange element
				at order (K) has N_p = (K+1)(K+2)/2 degrees of freedom, so an order (K-1) element has (K)(K+1)/2 DOF.
				Considering that we need a term for each of the two interior directions at each interior point, we need
				exactly 2*N_p DOF at order (K-1) for the interior of the custom RT element, resulting in (K)(K+1) terms.
				2) Note (1) confirms that the custom element requires exactly the same number of interior points
				(K)(K+1)/2 as a Lagrange element of order (K-1), which means we can use the custom RT element for the
				DFR approach, which needs to provide a O(K+1) element to preserve the gradient at O(K). We will use the
				solution points from the Lagrange element at O(K) to construct the interior of the O(K+1) RT element
				without requiring interpolation of the solution points, as they already reside at the same geometric
				locations.
				(3) To create the custom RT element, we initialize the solution element, then define the custom RT element
				from the interior point locations of the solution element to ensure that they are colocated.
				(4) To use the custom RT element:
				a) calculate the solution, calculate the flux vector field from the solution at the solution points
				b) transfer the flux vector field values to the DFR element interior
				c) interpolate flux values at from the interior of the RT element to the locations on the triangle edges
				d) use the method of characteristics to calculate the corrected flux using the neighbor element's edge
				flux combined with the edge flux from this element
				e) calculate the gradient of the vector flux field using the custom RT element
				f) transfer the gradient values from the RT element to the solution element for use in advancing the
				solution in differential form (directly)

			By calculating the flux gradient in a way that yields an O(K) polynomial on the solution points, we can use
			the differential form of the equations directly for the solution, rather than using the traditional Galerkin
			approach of repeated integration by parts to obtain an equation with only first derivatives. This simplifies
			the solution process, resulting in a more efficient computational approach, in addition to making it easier
			to solve more complex equations with the identical formulation.
	*/
	// Compute nodal set
	R, S := NodesEpsilon(el.N)
	el.RT = NewRTElement(el.N+1, R, S)
}

func (el *NDG_Elements2D) Startup2D() {
	var (
		err error
	)
	el.Nfp = el.N + 1
	el.Np = (el.N + 1) * (el.N + 2) / 2
	el.NFaces = 3
	el.NODETOL = 1.e-12
	// Compute nodal set
	//el.R, el.S = NodesEpsilon(el.N)
	el.R, el.S = XYtoRS(Nodes2D(el.N))
	// Build reference element matrices
	el.V = Vandermonde2D(el.N, el.R, el.S)
	if el.Vinv, err = el.V.Inverse(); err != nil {
		panic(err)
	}
	el.MassMatrix = el.Vinv.Transpose().Mul(el.Vinv)
	// Initialize the (r,s) differentiation matrices on the simplex, evaluated at (r,s) at order N
	Vr, Vs := GradVandermonde2D(el.N, el.R, el.S)
	el.Dr = Vr.Mul(el.Vinv)
	el.Ds = Vs.Mul(el.Vinv)

	// build coordinates of all the nodes
	va, vb, vc := el.EToV.Col(0), el.EToV.Col(1), el.EToV.Col(2)
	el.X = el.R.Copy().Add(el.S).Scale(-1).Outer(el.VX.SubsetIndex(va.ToIndex())).Add(
		el.R.Copy().AddScalar(1).Outer(el.VX.SubsetIndex(vb.ToIndex()))).Add(
		el.S.Copy().AddScalar(1).Outer(el.VX.SubsetIndex(vc.ToIndex()))).Scale(0.5)
	el.Y = el.R.Copy().Add(el.S).Scale(-1).Outer(el.VY.SubsetIndex(va.ToIndex())).Add(
		el.R.Copy().AddScalar(1).Outer(el.VY.SubsetIndex(vb.ToIndex()))).Add(
		el.S.Copy().AddScalar(1).Outer(el.VY.SubsetIndex(vc.ToIndex()))).Scale(0.5)
	fmask1 := el.S.Copy().AddScalar(1).Find(utils.Less, el.NODETOL, true)
	fmask2 := el.S.Copy().Add(el.R).Find(utils.Less, el.NODETOL, true)
	fmask3 := el.R.Copy().AddScalar(1).Find(utils.Less, el.NODETOL, true)
	if fmask1.Len() != 0 {

		el.FMask = utils.NewMatrix(el.Nfp, 3)
		el.FMask.SetCol(0, fmask1.Data())
		el.FMask.SetCol(1, fmask2.Data())
		el.FMask.SetCol(2, fmask3.Data())
		el.FMaskI = utils.NewIndex(len(el.FMask.Data()), el.FMask.Data())
		el.Fx = utils.NewMatrix(3*el.Nfp, el.K)
		for fp, val := range el.FMask.Data() {
			ind := int(val)
			el.Fx.M.SetRow(fp, el.X.M.RawRowView(ind))
		}
		el.Fy = utils.NewMatrix(3*el.Nfp, el.K)
		for fp, val := range el.FMask.Data() {
			ind := int(val)
			el.Fy.M.SetRow(fp, el.Y.M.RawRowView(ind))
		}
		el.Lift2D()
	}
	el.GeometricFactors2D()
	el.Normals2D()
	el.FScale = el.sJ.ElDiv(el.J.Subset(el.GetFaces()))
	// Build connectivity matrix
	el.Connect2D()

	// Mark fields read only
	el.Dr.SetReadOnly("Dr")
	el.Ds.SetReadOnly("Ds")
	el.LIFT.SetReadOnly("LIFT")
	el.X.SetReadOnly("X")
	el.Y.SetReadOnly("Y")
	el.Fx.SetReadOnly("Fx")
	el.Fy.SetReadOnly("Fy")
	el.FMask.SetReadOnly("FMask")
	el.MassMatrix.SetReadOnly("MassMatrix")
	el.V.SetReadOnly("V")
	el.Vinv.SetReadOnly("Vinv")
	el.NX.SetReadOnly("NX")
	el.NY.SetReadOnly("NY")
	el.FScale.SetReadOnly("FScale")
	el.EToE.SetReadOnly("EToE")
	el.EToF.SetReadOnly("EToF")
	return
}

/*
	Startup2D
  // Build connectivity maps
  BuildMaps2D();
  // Compute weak operators (could be done in preprocessing to save time)
  DMat Vr,Vs;  GradVandermonde2D(N, r, s, Vr, Vs);
  VVT = V*trans(V);
  Drw = (V*trans(Vr))/VVT;  Dsw = (V*trans(Vs))/VVT;
*/
func (el *NDG_Elements2D) Connect2D() {
	var (
		Nv         = el.VX.Len()
		TotalFaces = el.NFaces * el.K
	)
	SpFToVDOK := utils.NewDOK(TotalFaces, Nv)
	faces := utils.NewMatrix(3, 2, []float64{
		0, 1,
		1, 2,
		0, 2,
	})
	var sk int
	for k := 0; k < el.K; k++ {
		for face := 0; face < el.NFaces; face++ {
			edge := faces.Range(face, ":")
			//fmt.Println("Nv, TotalFaces, k, face, edge, range = ", Nv, TotalFaces, k, face, edge, el.EToV.Range(k, edge))
			SpFToVDOK.Equate(1, sk, el.EToV.Range(k, edge))
			sk++
		}
	}
	// Build global face to global face sparse array
	SpFToV := SpFToVDOK.ToCSR()
	SpFToF := utils.NewCSR(TotalFaces, TotalFaces)
	SpFToF.M.Mul(SpFToV, SpFToV.T())
	for i := 0; i < TotalFaces; i++ {
		SpFToF.M.Set(i, i, SpFToF.At(i, i)-2)
	}
	// Find complete face to face connections
	F12 := utils.MatFind(SpFToF, utils.Equal, 2)

	element1 := F12.RI.Copy().Apply(func(val int) int { return val / el.NFaces })
	face1 := F12.RI.Copy().Apply(func(val int) int { return int(math.Mod(float64(val), float64(el.NFaces))) })

	element2 := F12.CI.Copy().Apply(func(val int) int { return val / el.NFaces })
	face2 := F12.CI.Copy().Apply(func(val int) int { return int(math.Mod(float64(val), float64(el.NFaces))) })

	// Rearrange into Nelements x Nfaces sized arrays
	el.EToE = utils.NewRangeOffset(1, el.K).Outer(utils.NewOnes(el.NFaces))
	el.EToF = utils.NewOnes(el.K).Outer(utils.NewRangeOffset(1, el.NFaces))
	var I2D utils.Index2D
	var err error
	nr, nc := el.EToE.Dims()
	if I2D, err = utils.NewIndex2D(nr, nc, element1, face1); err != nil {
		panic(err)
	}
	el.EToE.Assign(I2D.ToIndex(), element2)
	el.EToF.Assign(I2D.ToIndex(), face2)
}

func (el *NDG_Elements2D) GetFaces() (aI utils.Index, NFacePts, K int) {
	var (
		err      error
		allFaces utils.Index2D
	)
	NFacePts = el.Nfp * el.NFaces
	K = el.K
	allK := utils.NewRangeOffset(1, el.K)
	if allFaces, err = utils.NewIndex2D(el.Np, el.K, el.FMaskI, allK, true); err != nil {
		panic(err)
	}
	aI = allFaces.ToIndex()
	return
}

func (el *NDG_Elements2D) Normals2D() {
	var (
		f1, f2, f3 utils.Index2D
		err        error
	)
	allK := utils.NewRangeOffset(1, el.K)
	aI, NFacePts, _ := el.GetFaces()
	// interpolate geometric factors to face nodes
	fxr := el.xr.Subset(aI, NFacePts, el.K)
	fxs := el.xs.Subset(aI, NFacePts, el.K)
	fyr := el.yr.Subset(aI, NFacePts, el.K)
	fys := el.ys.Subset(aI, NFacePts, el.K)
	// build normals
	faces1 := utils.NewRangeOffset(1, el.Nfp)
	faces2 := utils.NewRangeOffset(1+el.Nfp, 2*el.Nfp)
	faces3 := utils.NewRangeOffset(1+2*el.Nfp, 3*el.Nfp)
	if f1, err = utils.NewIndex2D(NFacePts, el.K, faces1, allK, true); err != nil {
		panic(err)
	}
	if f2, err = utils.NewIndex2D(NFacePts, el.K, faces2, allK, true); err != nil {
		panic(err)
	}
	if f3, err = utils.NewIndex2D(NFacePts, el.K, faces3, allK, true); err != nil {
		panic(err)
	}
	el.NX, el.NY = utils.NewMatrix(NFacePts, el.K), utils.NewMatrix(NFacePts, el.K)
	// Face 1
	el.NX.Assign(f1.ToIndex(), fyr.Subset(f1.ToIndex(), f1.Len, 1))
	el.NY.Assign(f1.ToIndex(), fxr.Subset(f1.ToIndex(), f1.Len, 1).Scale(-1))
	// Face 2
	el.NX.Assign(f2.ToIndex(), fys.Subset(f2.ToIndex(), f2.Len, 1).Subtract(fyr.Subset(f2.ToIndex(), f2.Len, 1)))
	el.NY.Assign(f2.ToIndex(), fxs.Subset(f2.ToIndex(), f2.Len, 1).Scale(-1).Add(fxr.Subset(f2.ToIndex(), f2.Len, 1)))
	// Face 3
	el.NX.Assign(f3.ToIndex(), fys.Subset(f3.ToIndex(), f3.Len, 1).Scale(-1))
	el.NY.Assign(f3.ToIndex(), fxs.Subset(f3.ToIndex(), f3.Len, 1))
	el.sJ = el.NX.Copy().POW(2).Add(el.NY.Copy().POW(2)).Apply(func(val float64) (res float64) {
		res = math.Sqrt(val)
		return
	})
	el.NX.ElDiv(el.sJ)
	el.NY.ElDiv(el.sJ)
}

func (el *NDG_Elements2D) GeometricFactors2D() {
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
	el.xr, el.xs = el.Dr.Mul(el.X), el.Ds.Mul(el.X)
	el.yr, el.ys = el.Dr.Mul(el.Y), el.Ds.Mul(el.Y)
	el.xr.SetReadOnly("xr")
	el.xs.SetReadOnly("xs")
	el.yr.SetReadOnly("yr")
	el.ys.SetReadOnly("ys")
	el.J = el.xr.Copy().ElMul(el.ys).Subtract(el.xs.Copy().ElMul(el.yr))
	el.Rx = el.ys.Copy().ElDiv(el.J)
	el.Sx = el.yr.Copy().ElDiv(el.J).Scale(-1)
	el.Ry = el.xs.Copy().ElDiv(el.J).Scale(-1)
	el.Sy = el.xr.Copy().ElDiv(el.J)
}

func (el *NDG_Elements2D) Lift2D() {
	var (
		err      error
		I2       utils.Index2D
		massEdge utils.Matrix
		V1D      utils.Matrix
		Emat     utils.Matrix
	)
	Emat = utils.NewMatrix(el.Np, el.NFaces*el.Nfp)
	faceMap := func(basis utils.Vector, faceNum int, Ind utils.Index) {
		faceBasis := basis.SubsetIndex(el.FMask.Col(faceNum).ToIndex())
		V1D = DG1D.Vandermonde1D(el.N, faceBasis)
		if massEdge, err = V1D.Mul(V1D.Transpose()).Inverse(); err != nil {
			panic(err)
		}
		if I2, err = utils.NewIndex2D(el.Np, el.NFaces*el.Nfp, el.FMask.Col(faceNum).ToIndex(), Ind, true); err != nil {
			panic(err)
		}
		Emat.Assign(I2.ToIndex(), massEdge)
	}
	// face 1
	faceMap(el.R, 0, utils.NewRangeOffset(1, el.Nfp))
	// face 2
	faceMap(el.R, 1, utils.NewRangeOffset(el.Nfp+1, 2*el.Nfp))
	// face 3
	faceMap(el.S, 2, utils.NewRangeOffset(2*el.Nfp+1, 3*el.Nfp))
	// inv(mass matrix)*\I_n (L_i,L_j)_{edge_n}
	el.LIFT = el.V.Mul(el.V.Transpose().Mul(Emat))
	return
}

func (el *NDG_Elements2D) BuildMaps2D() {
	return
}

func (el *NDG_Elements2D) NewCube2D(COrder int) {
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
