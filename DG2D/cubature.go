package DG2D

import "fmt"

func Cubature2D(COrder int) {
	// function [cubR,cubS,cubW, Ncub] = Cubature2D(COrder)
	// Purpose: provide multidimensional quadrature (i.e. cubature)
	//          rules to integrate up to COrder polynomials

	if COrder > 28 {
		COrder = 28
	}

	if COrder <= 28 {
		/*
		   DMat RSW
		   Cub2D_data(COrder, RSW)
		   cub.r = RSW(All,1)
		   cub.s = RSW(All,2)
		   cub.w = RSW(All,3)
		   cub.Ncub = RSW.num_rows()
		*/
	} else {
		fmt.Printf("Cubature2D(%d): COrder > 28 not yet tested\n", COrder)

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
}

/*

 */

/*
//---------------------------------------------------------
Cub2D& NDG2D::CubatureVolumeMesh2D(int COrder)
//---------------------------------------------------------
{
  // function cub = CubatureVolumeMesh2D(COrder)
  // purpose: build cubature nodes, weights and geometric factors for all elements
  //
  // Note: m_cub is member of Globals2D

  // set up cubature nodes
  Cubature2D(COrder, m_cub);

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

#include "NDGLib_headers.h"
#include "CubatureData2D.h"

void Cub2D_data(int Cn, DMat& Cdata);

//---------------------------------------------------------
void Cubature2D(int COrder, Cub2D& cub)
//---------------------------------------------------------
{
  // function [cubR,cubS,cubW, Ncub] = Cubature2D(COrder)
  // Purpose: provide multidimensional quadrature (i.e. cubature)
  //          rules to integrate up to COrder polynomials

  COrder = std::min(28, COrder);

  if (COrder<=28) {
    DMat RSW;  Cub2D_data(COrder, RSW);
    cub.r = RSW(All,1);
    cub.s = RSW(All,2);
    cub.w = RSW(All,3);
    cub.Ncub = RSW.num_rows();
  } else {

    umMSG(1, "Cubature2D(%d): COrder > 28 not yet tested\n", COrder);

    DVec cuba,cubwa, cubb,cubwb;
    DMat cubA, cubB, cubR, cubS, cubW, tA,tB;

    int cubNA=(int)ceil((COrder+1.0)/2.0);
    int cubNB=(int)ceil((COrder+1.0)/2.0);

    JacobiGQ(0.0, 0.0, cubNA-1,  cuba,cubwa);
    JacobiGQ(1.0, 0.0, cubNB-1,  cubb,cubwb);

    cubA = outer( ones(cubNB), cuba );
    cubB = outer( cubb, ones(cubNA) );

    tA = 1.0+cubA; tB = 1.0-cubB;
    cubR = 0.5 * tA.dm(tB) - 1.0;
    cubS = cubB;
    cubW = 0.5 * outer(cubwb, cubwa);

    cub.r = cubR;
    cub.s = cubS;
    cub.w = cubW;
    cub.Ncub = cub.r.size();
  }
}


//---------------------------------------------------------
void Cubature2D(int COrder, DVec& r, DVec& s, DVec& w, int& Ncub)
//---------------------------------------------------------
{
  Cub2D cub; Cubature2D(COrder, cub);
  r = cub.r; s = cub.s; w = cub.w;
  Ncub = cub.Ncub;
}

}
*/
