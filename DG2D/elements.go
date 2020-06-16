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
	R, VX, VY, VZ, FMask              utils.Vector
	EToV, EToE, EToF                  utils.Matrix
	BCType                            utils.Matrix
	X, Dr, Rx, FScale, NX, LIFT       utils.Matrix
	V, Vinv                           utils.Matrix
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

func NewElements2D(N, K int, meshFile string, plotMesh bool) (el *Elements2D) {
	var (
		// choose order to integrate exactly
		CubatureOrder = int(math.Floor(2.0 * float64(N+1) * 3.0 / 2.0))
		NGauss        = int(math.Floor(2.0 * float64(N+1)))
	)
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
	el = &Elements2D{
		N:      N,
		K:      K,
		Np:     (N + 1) * (N + 2) / 2,
		NFaces: 3,
	}
	el.ReadGambit2d(meshFile, plotMesh)
	el.NewCube2D(CubatureOrder)
	el.Startup2D()
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
		}
		sk++
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
	tv1, tv2 := make([]float64, Np), make([]float64, Np)
	P = make([]float64, Np)
	for i, h1Val := range h1 {
		tv1[i] = h1Val * h2[i] * math.Sqrt(2)
		tv2[i] = utils.POW(1-bd[i], i)
		P[i] = tv1[i] * tv2[i]
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
	//DVec blend1,blend2,blend3,warp1,warp2,warp3,warpf1,warpf2,warpf3;
	var (
		alpha                                                                           float64
		Np                                                                              = (N + 1) * (N + 2) / 2
		L1, L2, L3, blend1, blend2, blend3, warp1, warp2, warp3, warpf1, warpf2, warpf3 utils.Vector
	)
	L1, L2, L3, blend1, blend2, blend3, warp1, warp2, warp3, warpf1, warpf2, warpf3, x, y =
		utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np),
		utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np),
		utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np), utils.NewVector(Np)

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

	_, _, _, _, _, _, _, _, _, _ = alpha, blend1, blend2, blend3, warp1, warp2, warp3, warpf1, warpf2, warpf3
	// Create equidistributed nodes on equilateral triangle
	l1d, l2d, l3d, xd, yd := L1.Data(), L2.Data(), L3.Data(), x.Data(), y.Data()
	fn := 1. / float64(N)
	var sk int
	for n := 0; n < N+1; n++ {
		for m := 0; m < (N + 2 - n); m++ {
			l1d[sk] = float64(n-1) * fn
			l3d[sk] = float64(m-1) * fn
			l2d[sk] = 1 - l1d[sk] - l3d[sk]
			xd[sk] = l3d[sk] - l2d[sk]
			yd[sk] = (-l3d[sk] - l2d[sk] + 2*l1d[sk]) / math.Sqrt(3)
			sk++
		}
	}
	_, _ = xd, yd
	/*
	     // Compute blending function at each node for each edge
	     blend1 = 4.0 * L2.dm(L3);
	     blend2 = 4.0 * L1.dm(L3);
	     blend3 = 4.0 * L1.dm(L2);

	     // Amount of warp for each node, for each edge
	     warpf1 = Warpfactor(N, L3-L2);
	     warpf2 = Warpfactor(N, L1-L3);
	     warpf3 = Warpfactor(N, L2-L1);

	     // Combine blend & warp
	     warp1 = blend1.dm(warpf1);  warp1 *= (1.0 + sqr(alpha*L1));
	     warp2 = blend2.dm(warpf2);  warp2 *= (1.0 + sqr(alpha*L2));
	     warp3 = blend3.dm(warpf3);  warp3 *= (1.0 + sqr(alpha*L3));

	     // Accumulate deformations associated with each edge
	     x += (1.0*warp1 + cos(2.0*PI/3.0)*warp2 + cos(4.0*PI/3.0)*warp3);
	     y += (0.0*warp1 + sin(2.0*PI/3.0)*warp2 + sin(4.0*PI/3.0)*warp3);
	   }
	*/
	return
}

func Warpfactor(N int, rout utils.Vector) (warp utils.Vector) {
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
	fmt.Println(Lmat.Print("Lmat"))
	_ = LGLr
	/*
	     Lmat = trans(Veq)|Pmat;

	     // Compute warp factor
	     warp = trans(Lmat)*(LGLr - req);

	     // Scale factor
	     zerof = rout.lt_abs(1.0 - 1e-10); sf = 1.0 - sqr(zerof.dm(rout));
	     warp = warp.dd(sf) + warp.dm(zerof-1.0);

	     warp.set_mode(OBJ_temp);  // adjust mode flag
	     return warp;
	   }
	*/
	return
}

func (el *Elements2D) Startup2D() {
	el.Nfp = el.N
	el.Np = (el.N + 1) * (el.N + 2) / 2
	el.NFaces = 3
	el.NODETOL = 1.e-12
	/*
	     // Compute nodal set
	     DVec x1,y1; Nodes2D(N, x1,y1);  xytors(x1,y1, r,s);

	     // Build reference element matrices
	     V = Vandermonde2D(N,r,s); invV = inv(V);
	     MassMatrix = trans(invV)*invV;
	     ::Dmatrices2D(N,r,s,V, Dr,Ds);

	     // build coordinates of all the nodes
	     IVec va = EToV(All,1), vb = EToV(All,2), vc = EToV(All,3);

	     // Note: outer products of (Vector,MappedRegion1D)
	     x = 0.5 * (-(r+s)*VX(va) + (1.0+r)*VX(vb) + (1.0+s)*VX(vc));
	     y = 0.5 * (-(r+s)*VY(va) + (1.0+r)*VY(vb) + (1.0+s)*VY(vc));

	     // find all the nodes that lie on each edge
	     IVec fmask1,fmask2,fmask3;
	     fmask1 = find( abs(s+1.0), '<', NODETOL);
	     fmask2 = find( abs(r+s  ), '<', NODETOL);
	     fmask3 = find( abs(r+1.0), '<', NODETOL);
	     Fmask.resize(Nfp,3);                    // set shape (M,N) before concat()
	     Fmask = concat(fmask1,fmask2,fmask3);   // load vector into shaped matrix

	     Fx = x(Fmask, All); Fy = y(Fmask, All);

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

	     return true;
	   }
	*/
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
