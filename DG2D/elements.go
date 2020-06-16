package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

type Elements2D struct {
	K, N, Np, NFaces                  int
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
func (el *Elements2D) Startup2D() {
	/*
	   //---------------------------------------------------------
	   bool NDG2D::StartUp2D()
	   //---------------------------------------------------------
	   {
	     // Purpose : Setup script, building operators, grid, metric,
	     //           and connectivity tables.

	     // Definition of constants
	     Nfp = N+1; Np = (N+1)*(N+2)/2; Nfaces=3; NODETOL = 1e-12;

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
