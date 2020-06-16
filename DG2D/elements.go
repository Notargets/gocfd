package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
)

type Elements2D struct {
	K, Np, NFaces                     int
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
		K:      K,
		Np:     (N + 1) * (N + 2) / 2,
		NFaces: 3,
	}
	el.ReadGambit2d(meshFile, plotMesh)
	el.NewCube2D(CubatureOrder)
	el.Startup2D()
	return
}

/*
	//---------------------------------------------------------
	Cub2D& NDG2D::CubatureVolumeMesh2D(int Corder)
	//---------------------------------------------------------
	{
	  // function cub = CubatureVolumeMesh2D(Corder)
	  // purpose: build cubature nodes, weights and geometric factors for all elements
	  //
	  // Note: m_cub is member of Globals2D

	  // set up cubature nodes
	  Cubature2D(Corder, m_cub);

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
