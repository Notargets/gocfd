package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

type Elements2D struct {
	K, N, Nfp, Np, NFaces             int
	NODETOL                           float64
	R, S, VX, VY, VZ                  utils.Vector
	FMask, Fx, Fy                     utils.Matrix
	EToV, EToE, EToF                  utils.Matrix
	BCType                            utils.Matrix
	X, Y, Dr, Ds, FScale, NX, LIFT    utils.Matrix
	Rx, Ry, Sx, Sy                    utils.Matrix
	V, Vinv, MassMatrix               utils.Matrix
	J                                 utils.Matrix
	VmapM, VmapP, VmapB, VmapI, VmapO utils.Index
	MapB, MapI, MapO                  utils.Index
	Cub                               *Cubature
}

type Cubature struct {
	r, s, w                 utils.Vector
	W                       utils.Matrix
	V, Dr, Ds, VT, DrT, DsT utils.Matrix
	x, y, rx, sx, ry, sy, J utils.Matrix
	mm, mmCHOL              utils.Matrix
}

func NewElements2D(N int, meshFile string, plotMesh bool) (el *Elements2D) {
	var (
		// choose order to integrate exactly
		CubatureOrder = int(math.Floor(2.0 * float64(N+1) * 3.0 / 2.0))
		NGauss        = int(math.Floor(2.0 * float64(N+1)))
	)
	if N < 1 {
		N = 1
	}
	el = &Elements2D{
		N:      N,
		Np:     (N + 1) * (N + 2) / 2,
		NFaces: 3,
	}
	el.ReadGambit2d(meshFile, plotMesh)
	el.NewCube2D(CubatureOrder)
	el.Startup2D()
	_ = NGauss
	/*
	  // build cubature node data for all elements
	  CubatureVolumeMesh2D(CubatureOrder);

	  // build Gauss node data for all element faces
	  GaussFaceMesh2D(NGauss);

	  Resize_cub();           // resize cubature arrays
	  MapGaussFaceData();     // {nx = gauss.nx}, etc.
	  PreCalcBdryData();      // gmapB = concat(mapI, mapO), etc.
	*/
	// N is the polynomial degree, Np is the number of interpolant points = N+1

	return
}

func (el *Elements2D) InterpMatrix2D() {
	/*
	   //---------------------------------------------------------
	   void NDG2D::InterpMatrix2D(Cub2D& cub)
	   //---------------------------------------------------------
	   {
	   // compute Vandermonde at (rout,sout)
	   DMat Vout = Vandermonde2D(this->N, cub.r, cub.s);
	   // build interpolation matrix
	   cub.V = Vout * this->invV;
	   // store transpose
	   cub.VT = trans(cub.V);
	   }
	   }
	*/
}

func Vandermonde2D(N int, r, s utils.Vector) (V2D utils.Matrix) {
	V2D = utils.NewMatrix(r.Len(), (N+1)*(N+2)/2)
	a, b := RStoAB(r, s)
	var sk int
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			V2D.SetCol(sk, Simplex2DP(a, b, i, j))
			sk++
		}
	}
	return
}

func Simplex2DP(a, b utils.Vector, i, j int) (P []float64) {
	var (
		Np = a.Len()
		bd = b.Data()
	)
	h1 := DG1D.JacobiP(a, 0, 0, i)
	h2 := DG1D.JacobiP(b, float64(2*i+1), 0, j)
	P = make([]float64, Np)
	sq2 := math.Sqrt(2)
	for ii := range h1 {
		tv1 := sq2 * h1[ii] * h2[ii]
		tv2 := utils.POW(1-bd[ii], i)
		P[ii] = tv1 * tv2
	}
	return
}

/*
	  // evaluate generalized Vandermonde of Lagrange interpolant functions at cubature nodes
	  InterpMatrix2D(m_cub);

	  // evaluate local derivatives of Lagrange interpolants at cubature nodes
	  Dmatrices2D(this->N, m_cub);

	  // evaluate the geometric factors at the cubature nodes
	  GeometricFactors2D(m_cub);

	  // custom mass matrix per element
	  DMat mmk; DMat_Diag D; DVec d;
	  m_cub.mmCHOL.resize(Np*Np, K);
	  m_cub.mm    .resize(Np*Np, K);

	  for (int k=1; k<=K; ++k) {
	    d=m_cub.J(All,k); d*=m_cub.w; D.diag(d);  // weighted diagonal
	    mmk = m_cub.VT * D * m_cub.V;     // mass matrix for element k
	    m_cub.mm(All,k)     = mmk;        // store mass matrix
	    m_cub.mmCHOL(All,k) = chol(mmk);  // store Cholesky factorization
	  }

	  // incorporate weights and Jacobian
	  m_cub.W = outer(m_cub.w, ones(K));
	  m_cub.W.mult_element(m_cub.J);

	  // compute coordinates of cubature nodes
	  m_cub.x = m_cub.V * this->x;
	  m_cub.y = m_cub.V * this->y;

	  return m_cub;
	}
*/

// Purpose  : Compute (x,y) nodes in equilateral triangle for
//            polynomial of order N
func Nodes2D(N int) (x, y utils.Vector) {
	var (
		alpha                                                               float64
		Np                                                                  = (N + 1) * (N + 2) / 2
		L1, L2, L3                                                          utils.Vector
		blend1, blend2, blend3, warp1, warp2, warp3, warpf1, warpf2, warpf3 []float64
	)
	L1, L2, L3, x, y =
		utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np)
	l1d, l2d, l3d, xd, yd := L1.Data(), L2.Data(), L3.Data(), x.Data(), y.Data()
	blend1, blend2, blend3, warp1, warp2, warp3 =
		make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np), make([]float64, Np)

	alpopt := [15]float64{
		0.0000, 0.0000, 1.4152, 0.1001, 0.2751,
		0.9800, 1.0999, 1.2832, 1.3648, 1.4773,
		1.4959, 1.5743, 1.5770, 1.6223, 1.6258,
	}
	if N < 16 {
		alpha = alpopt[N]
	} else {
		alpha = 5. / 3.
	}
	// Create equidistributed nodes on equilateral triangle
	fn := 1. / float64(N)
	var sk int
	for n := 0; n < N+1; n++ {
		for m := 0; m < (N + 1 - n); m++ {
			l1d[sk] = float64(n) * fn
			l3d[sk] = float64(m) * fn
			l2d[sk] = 1 - l1d[sk] - l3d[sk]
			xd[sk] = l3d[sk] - l2d[sk]
			yd[sk] = (-l3d[sk] - l2d[sk] + 2*l1d[sk]) / math.Sqrt(3)
			// Compute blending function at each node for each edge
			blend1[sk] = 4 * l3d[sk] * l2d[sk]
			blend2[sk] = 4 * l1d[sk] * l3d[sk]
			blend3[sk] = 4 * l2d[sk] * l1d[sk]
			sk++
		}
	}
	// Amount of warp for each node, for each edge
	warpf1 = Warpfactor(N, L3.Copy().Subtract(L2))
	warpf2 = Warpfactor(N, L1.Copy().Subtract(L3))
	warpf3 = Warpfactor(N, L2.Copy().Subtract(L1))
	// Combine blend & warp
	for i := range warpf1 {
		warp1[i] = blend1[i] * warpf1[i] * (1 + math.Sqrt(alpha*l1d[i]))
		warp2[i] = blend2[i] * warpf2[i] * (1 + math.Sqrt(alpha*l2d[i]))
		warp3[i] = blend3[i] * warpf3[i] * (1 + math.Sqrt(alpha*l3d[i]))
	}
	// Accumulate deformations associated with each edge
	for i := range xd {
		xd[i] += (warp1[i]+math.Cos(2*math.Pi/3))*warp2[i] + math.Cos(4*math.Pi/3)*warp3[i]
		yd[i] += math.Sin(2*math.Pi/3)*warp2[i] + math.Sin(4*math.Pi/3)*warp3[i]
	}
	return
}

func Warpfactor(N int, rout utils.Vector) (warpF []float64) {
	var (
		Nr   = rout.Len()
		Pmat = utils.NewMatrix(N+1, Nr)
	)
	// Compute LGL and equidistant node distribution
	LGLr := DG1D.JacobiGL(0, 0, N)
	req := utils.NewVector(N+1).Linspace(-1, 1)
	Veq := DG1D.Vandermonde1D(N, req)
	// Evaluate Lagrange polynomial at rout
	for i := 0; i < (N + 1); i++ {
		Pmat.M.SetRow(i, DG1D.JacobiP(rout, 0, 0, i))
	}
	Lmat := Veq.Transpose().LUSolve(Pmat)
	// Compute warp factor
	warp := Lmat.Transpose().Mul(LGLr.Subtract(req).ToMatrix())
	// Scale factor
	zerof := rout.Apply(func(val float64) (res float64) {
		if math.Abs(val) < 1.0 {
			res = 1
		} else {
			res = val
		}
		res -= 1.e-10
		return
	})
	sf := zerof.Copy().ElMul(rout).Apply(func(val float64) (res float64) {
		res = 1 - val*val
		return
	})
	warp.ElDiv(sf.ToMatrix()).ElMul(zerof.AddScalar(-1).ToMatrix())
	warpF = warp.Data()
	return
}

func (el *Elements2D) Startup2D() {
	var (
		err error
	)
	el.Nfp = el.N + 1
	el.Np = (el.N + 1) * (el.N + 2) / 2
	el.NFaces = 3
	el.NODETOL = 1.e-12
	// Compute nodal set
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
	el.FMask = utils.NewMatrix(el.Nfp, 3)
	el.FMask.SetCol(0, fmask1.Data())
	el.FMask.SetCol(1, fmask2.Data())
	el.FMask.SetCol(2, fmask3.Data())
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
	el.GeometricFactors2D()
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
	return
}

func (el *Elements2D) GeometricFactors2D() {
	var (
		xr, xs = el.Dr.Mul(el.X), el.Ds.Mul(el.X)
		yr, ys = el.Dr.Mul(el.Y), el.Ds.Mul(el.Y)
	)
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
	el.J = xr.Copy().ElMul(ys).Subtract(xs.Copy().ElMul(yr))
	el.Rx = ys.Copy().ElDiv(el.J)
	el.Sx = yr.Copy().ElDiv(el.J).Scale(-1)
	el.Ry = xs.Copy().ElDiv(el.J).Scale(-1)
	el.Sy = xr.Copy().ElDiv(el.J)
}

func (el *Elements2D) Lift2D() {
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

/*
	Startup2D
  // Create surface integral terms
  Lift2D();

  // calculate geometric factors
  ::GeometricFactors2D(x,y,Dr,Ds,  rx,sx,ry,sy,J);

  // calculate geometric factors
  Normals2D();
  Fscale = sJ.dd(J(Fmask,All));


#if (0)
  OutputNodes(false); // volume nodes
  OutputNodes(true);  // face nodes
  umERROR("Exiting early", "Check {volume,face} nodes");
#endif

  // Build connectivity matrix
  tiConnect2D(EToV, EToE,EToF);

  // Build connectivity maps
  BuildMaps2D();

  // Compute weak operators (could be done in preprocessing to save time)
  DMat Vr,Vs;  GradVandermonde2D(N, r, s, Vr, Vs);
  VVT = V*trans(V);
  Drw = (V*trans(Vr))/VVT;  Dsw = (V*trans(Vs))/VVT;
*/

func GradVandermonde2D(N int, r, s utils.Vector) (V2Dr, V2Ds utils.Matrix) {
	var (
		a, b = RStoAB(r, s)
		Np   = (N + 1) * (N + 2) / 2
		Nr   = r.Len()
	)
	V2Dr, V2Ds = utils.NewMatrix(Nr, Np), utils.NewMatrix(Nr, Np)
	var sk int
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			ddr, dds := GradSimplex2DP(a, b, i, j)
			V2Dr.M.SetCol(sk, ddr)
			V2Ds.M.SetCol(sk, dds)
			sk++
		}
	}
	return
}

func GradSimplex2DP(a, b utils.Vector, id, jd int) (ddr, dds []float64) {
	var (
		ad, bd = a.Data(), b.Data()
	)
	_ = ad
	fa := DG1D.JacobiP(a, 0, 0, id)
	dfa := DG1D.GradJacobiP(a, 0, 0, id)
	gb := DG1D.JacobiP(b, 2*float64(id)+1, 0, jd)
	dgb := DG1D.GradJacobiP(b, 2*float64(id)+1, 0, jd)
	// r-derivative
	// d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
	ddr = make([]float64, len(gb))
	for i := range ddr {
		ddr[i] = dfa[i] * gb[i]
		if id > 0 {
			ddr[i] *= utils.POW(0.5*(1-bd[i]), id-1)
		}
		// Normalize
		ddr[i] *= math.Pow(2, float64(id)+0.5)
	}
	// s-derivative
	// d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
	dds = make([]float64, len(gb))
	for i := range dds {
		dds[i] = 0.5 * dfa[i] * gb[i] * (1 + ad[i])
		if id > 0 {
			dds[i] *= utils.POW(0.5*(1-bd[i]), id-1)
		}
		tmp := dgb[i] * utils.POW(0.5*(1-bd[i]), id)
		if id > 0 {
			tmp -= 0.5 * float64(id) * gb[i] * utils.POW(0.5*(1-bd[i]), id-1)
		}
		dds[i] += fa[i] * tmp
		// Normalize
		dds[i] *= math.Pow(2, float64(id)+0.5)
	}
	return
}

func (el *Elements2D) Connect2D() {
	return
}

func (el *Elements2D) BuildMaps2D() {
	return
}

func (el *Elements2D) NewCube2D(COrder int) {
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

func RStoAB(r, s utils.Vector) (a, b utils.Vector) {
	var (
		Np     = r.Len()
		rd, sd = r.Data(), s.Data()
	)
	ad, bd := make([]float64, Np), make([]float64, Np)
	for n, sval := range sd {
		if sval != 1 {
			ad[n] = 2*(1+rd[n])/(1-sval) - 1
		} else {
			ad[n] = -1
		}
		bd[n] = sval
	}
	a, b = utils.NewVector(Np, ad), utils.NewVector(Np, bd)
	return
}

// function [r,s] = xytors(x,y)
// Purpose : Transfer from (x,y) in equilateral triangle
//           to (r,s) coordinates in standard triangle
func XYtoRS(x, y utils.Vector) (r, s utils.Vector) {
	r, s = utils.NewVector(x.Len()), utils.NewVector(x.Len())
	var (
		xd, yd = x.Data(), y.Data()
		rd, sd = r.Data(), s.Data()
	)
	sr3 := math.Sqrt(3)
	for i := range xd {
		l1 := (sr3*yd[i] + 1) / 3
		l2 := (-3*xd[i] - sr3*yd[i] + 2) / 6
		l3 := (3*xd[i] - sr3*yd[i] + 2) / 6
		rd[i] = -l2 + l3 - l1
		sd[i] = -l2 - l3 + l1
	}
	return
}
