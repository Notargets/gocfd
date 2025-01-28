package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

type RTBasis2DSimplex struct {
	P             int            // Polynomial Order
	Np            int            // Number of terms and nodes in basis
	NpInt, NpEdge int            // Number of nodes in interior and on each edge, Np = 2*NpInt + 3*NpEdge
	R, S          utils.Vector   // Node locations within [-1,1] reference triangle
	Scalar2DBasis *JacobiBasis2D // Basis used for part of the RT basis
	// construction
	V           [2]utils.Matrix
	Div, DivInt utils.Matrix
}

/*
This is a second Raviart-Thomas element basis implementation within this project.

This implementation is based on the paper "Computational Bases for RTk and BDMk on Triangles" by V. J. Ervin

In this approach, edges are specifically addressed with 1D Lagrange polynomials multiplied by edge specific basis
functions, while the interior points are composed of a supplied 2D scalar polynomial multiplied by a barycentric basis
for the triangle.
*/
func NewRTBasis2DSimplex(P int) (rtb *RTBasis2DSimplex) {
	var Rint, Sint utils.Vector
	fmt.Printf("Order of RT Element %d\n", P)
	rtb = &RTBasis2DSimplex{
		P:      P,
		Np:     (P + 1) * (P + 3),
		NpInt:  P * (P + 1) / 2, // Number of interior points is same as the 2D scalar space one order lesser
		NpEdge: P + 1,           // Each edge is P+1 nodes
	}
	if P > 0 {
		if P < 9 {
			Rint, Sint = NodesEpsilon(P - 1)
		} else {
			Rint, Sint = XYtoRS(Nodes2D(P - 1))
		}
		rtb.Scalar2DBasis = NewJacobiBasis2D(P-1, Rint, Sint, 0, 0) // Basis used as part of non-Normal basis elements
	}
	rtb.R, rtb.S = rtb.ExtendGeomToRT(Rint, Sint)
	rtb.CalculateBasis()
	return
}

func (rtb *RTBasis2DSimplex) getLocationType(i int) (locationType RTFunctionNumber) {
	var (
		NpInt  = rtb.NpInt
		NpEdge = rtb.NpEdge
	)
	switch {
	case i < NpInt:
		// Unit vector is [1,0]
		locationType = E4
	case i >= NpInt && i < 2*NpInt:
		// Unit vector is [0,1]
		locationType = E5
	case i >= 2*NpInt && i < 2*NpInt+NpEdge:
		// E1: Unit vector is [0,-1]
		locationType = E1
	case i >= 2*NpInt+NpEdge && i < 2*NpInt+2*NpEdge:
		// E2: Unit vector is [1/sqrt(2), 1/sqrt(2)]
		locationType = E2
	case i >= 2*NpInt+2*NpEdge && i < 2*NpInt+3*NpEdge:
		// E3: Unit vector is [-1,0]
		locationType = E3
	default:
		fmt.Errorf("unable to match node location for node %d and Np = %d", i, rtb.Np)
	}
	return
}

func (rtb *RTBasis2DSimplex) EvaluateBasisAtLocation(r, s float64, derivO ...DerivativeDirection) (p0, p1 []float64) {
	var (
		deriv  DerivativeDirection
		RGauss = DG1D.LegendreZeros(rtb.P)
	)
	if len(derivO) != 0 {
		deriv = derivO[0]
	} else {
		deriv = None
	}
	Basis2DTerm := func(r, s float64, i, j int, dv DerivativeDirection) (p float64) {
		switch dv {
		case None:
			p = rtb.Scalar2DBasis.PolynomialTerm(r, s, i, j)
		case Dr:
			p = rtb.Scalar2DBasis.PolynomialTermDr(r, s, i, j)
		case Ds:
			p = rtb.Scalar2DBasis.PolynomialTermDs(r, s, i, j)
		}
		return
	}
	e1rf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p, _ = rtb.getCoreBasisTerm(E1, r, s, deriv)
		return
	}
	e1sf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		_, p = rtb.getCoreBasisTerm(E1, r, s, deriv)
		return
	}
	e2rf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p, _ = rtb.getCoreBasisTerm(E2, r, s, deriv)
		return
	}
	e2sf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		_, p = rtb.getCoreBasisTerm(E2, r, s, deriv)
		return
	}
	e3rf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p, _ = rtb.getCoreBasisTerm(E3, r, s, deriv)
		return
	}
	e3sf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		_, p = rtb.getCoreBasisTerm(E3, r, s, deriv)
		return
	}
	e4rf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p, _ = rtb.getCoreBasisTerm(E4, r, s, deriv)
		return
	}
	e4sf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		_, p = rtb.getCoreBasisTerm(E4, r, s, deriv)
		return
	}
	e5rf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p, _ = rtb.getCoreBasisTerm(E5, r, s, deriv)
		return
	}
	e5sf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		_, p = rtb.getCoreBasisTerm(E5, r, s, deriv)
		return
	}
	l1f := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.LinearPoly(r, s, 0, deriv)
		return
	}
	l2f := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.LinearPoly(r, s, 1, deriv)
		return
	}
	l3f := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.LinearPoly(r, s, 2, deriv)
		return
	}
	q1rf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.Lagrange1DPoly(r, RGauss, 0, RDir, deriv)
		return
	}
	q2rf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.Lagrange1DPoly(r, RGauss, 1, RDir, deriv)
		return
	}
	q3rf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.Lagrange1DPoly(r, RGauss, 2, RDir, deriv)
		return
	}
	q1sf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.Lagrange1DPoly(r, RGauss, 0, SDir, deriv)
		return
	}
	q2sf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.Lagrange1DPoly(r, RGauss, 1, SDir, deriv)
		return
	}
	q3sf := func(r, s float64, deriv DerivativeDirection) (p float64) {
		p = rtb.Lagrange1DPoly(r, RGauss, 2, SDir, deriv)
		return
	}
	CRP := func(pL, dpL, pR, dpR float64) (dp float64) {
		// Chain Rule Product d(P*F)/dr = F*(dP/dr)+ P*(dF/dr)
		dp = dpL*pR + pL*dpR
		return
	}
	CRP1 := func(r, s float64, deriv DerivativeDirection, left, right func(r, s float64, deriv DerivativeDirection) (p float64)) (p float64) {
		l := left(r, s, None)
		ld := left(r, s, deriv)
		rt := right(r, s, None)
		rtd := right(r, s, deriv)
		p = CRP(l, ld, rt, rtd)
		return
	}
	p0, p1 = make([]float64, rtb.Np), make([]float64, rtb.Np)
	var sk int
	/*
		The first three element degrees, P=[0,1,2] are covered by special case code
		Degrees equal and higher to P=3 use a slightly different scalar basis than the special case code.
		We use our Lagrange or Jacobi basis for P=[3...], which means there's a difference between quadratic and
		higher elements that needs some documentation...
	*/
	// Non-normal basis functions first
	switch rtb.P {
	case 0:
		// No "non-Normal" basis functions
	case 1:
		// Two basis functions, basics e4 and e5
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(E4, r, s, deriv)
		sk++
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(E5, r, s, deriv)
		sk++
	case 2:
		// Six "non-Normal" basis functions, use linear 2D polynomial for triangles multiplied by e4 and e5
		if deriv == None {
			e4r, e4s := rtb.getCoreBasisTerm(E4, r, s)
			e5r, e5s := rtb.getCoreBasisTerm(E5, r, s)
			l1 := rtb.LinearPoly(r, s, 0)
			l2 := rtb.LinearPoly(r, s, 1)
			l3 := rtb.LinearPoly(r, s, 2)
			p0[sk], p1[sk] = e4r*l1, e4s*l1
			sk++
			p0[sk], p1[sk] = e4r*l2, e4s*l2
			sk++
			p0[sk], p1[sk] = e4r*l3, e4s*l3
			sk++

			p0[sk], p1[sk] = e5r*l1, e5s*l1
			sk++
			p0[sk], p1[sk] = e5r*l2, e5s*l2
			sk++
			p0[sk], p1[sk] = e5r*l3, e5s*l3
			sk++
		} else {
			// This covers either derivative direction, Dr or Ds
			p0[sk], p1[sk] = CRP1(r, s, deriv, e4rf, l1f), CRP1(r, s, deriv, e4sf, l1f)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e4rf, l2f), CRP1(r, s, deriv, e4sf, l2f)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e4rf, l3f), CRP1(r, s, deriv, e4sf, l3f)
			sk++

			p0[sk], p1[sk] = CRP1(r, s, deriv, e5rf, l1f), CRP1(r, s, deriv, e5sf, l1f)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e5rf, l2f), CRP1(r, s, deriv, e5sf, l2f)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e5rf, l3f), CRP1(r, s, deriv, e5sf, l3f)
			sk++
		}
	default:
		e4r, e4s := rtb.getCoreBasisTerm(E4, r, s)
		e5r, e5s := rtb.getCoreBasisTerm(E5, r, s)
		var (
			e4Rderiv, e4Sderiv, e5Rderiv, e5Sderiv float64
		)
		if deriv != None {
			e4Rderiv, e4Sderiv = rtb.getCoreBasisTerm(E4, r, s, deriv)
			e5Rderiv, e5Sderiv = rtb.getCoreBasisTerm(E5, r, s, deriv)
		}
		Pint := rtb.P - 1
		for i := 0; i <= Pint; i++ {
			for j := 0; j <= Pint-i; j++ {
				p := Basis2DTerm(r, s, i, j, None)
				if deriv == None {
					p0[sk], p1[sk] = p*e4r, p*e4s
					sk++
					p0[sk], p1[sk] = p*e5r, p*e5s
					sk++
				} else {
					dpdr := Basis2DTerm(r, s, i, j, deriv)
					p0[sk], p1[sk] = CRP(e4r, e4Rderiv, p, dpdr), CRP(e4s, e4Sderiv, p, dpdr)
					sk++
					p0[sk], p1[sk] = CRP(e5r, e5Rderiv, p, dpdr), CRP(e5s, e5Sderiv, p, dpdr)
					sk++
				}
			}
		}
	}
	/*
		Now we do the "Normal" basis functions associated with the edges
		These use the 1D Lagrange polynomials multiplied by the core basis functions
	*/
	switch rtb.P {
	case 0:
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(E1, r, s, deriv)
		sk++
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(E2, r, s, deriv)
		sk++
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(E3, r, s, deriv)
		sk++
	case 1:
		if deriv == None {
			e1r, e1s := rtb.getCoreBasisTerm(E1, r, s)
			e2r, e2s := rtb.getCoreBasisTerm(E2, r, s)
			e3r, e3s := rtb.getCoreBasisTerm(E3, r, s)
			l1xi := rtb.Lagrange1DPoly(r, RGauss, 0, RDir)
			l2xi := rtb.Lagrange1DPoly(r, RGauss, 1, RDir)
			l1eta := rtb.Lagrange1DPoly(s, RGauss, 0, SDir)
			l2eta := rtb.Lagrange1DPoly(s, RGauss, 1, SDir)
			p0[sk], p1[sk] = l1eta*e1r, l1eta*e1s
			sk++
			p0[sk], p1[sk] = l2eta*e1r, l2eta*e1s
			sk++

			p0[sk], p1[sk] = l2eta*e2r, l2eta*e2s
			sk++
			p0[sk], p1[sk] = l1eta*e2r, l1eta*e2s
			sk++

			p0[sk], p1[sk] = l1xi*e3r, l1xi*e3s
			sk++
			p0[sk], p1[sk] = l2xi*e3r, l2xi*e3s
			sk++
		} else {
			p0[sk], p1[sk] = CRP1(r, s, deriv, e1rf, q1sf), CRP1(r, s, deriv, e1sf, q1sf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e1rf, q2sf), CRP1(r, s, deriv, e1sf, q2sf)
			sk++

			p0[sk], p1[sk] = CRP1(r, s, deriv, e2rf, q2sf), CRP1(r, s, deriv, e2sf, q2sf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e2rf, q1sf), CRP1(r, s, deriv, e2sf, q1sf)
			sk++

			p0[sk], p1[sk] = CRP1(r, s, deriv, e3rf, q1rf), CRP1(r, s, deriv, e3sf, q1rf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e3rf, q2rf), CRP1(r, s, deriv, e3sf, q2rf)
			sk++
		}
	case 2:
		if deriv == None {
			e1r, e1s := rtb.getCoreBasisTerm(E1, r, s)
			e2r, e2s := rtb.getCoreBasisTerm(E2, r, s)
			e3r, e3s := rtb.getCoreBasisTerm(E3, r, s)
			q1xi := rtb.Lagrange1DPoly(r, RGauss, 0, RDir)
			q2xi := rtb.Lagrange1DPoly(r, RGauss, 1, RDir)
			q3xi := rtb.Lagrange1DPoly(r, RGauss, 2, RDir)
			q1eta := rtb.Lagrange1DPoly(s, RGauss, 0, SDir)
			q2eta := rtb.Lagrange1DPoly(s, RGauss, 1, SDir)
			q3eta := rtb.Lagrange1DPoly(s, RGauss, 2, SDir)
			p0[sk], p1[sk] = q1eta*e1r, q1eta*e1s
			sk++
			p0[sk], p1[sk] = q2eta*e1r, q2eta*e1s
			sk++
			p0[sk], p1[sk] = q3eta*e1r, q3eta*e1s
			sk++

			p0[sk], p1[sk] = q3eta*e2r, q3eta*e2s
			sk++
			p0[sk], p1[sk] = q2eta*e2r, q2eta*e2s
			sk++
			p0[sk], p1[sk] = q1eta*e2r, q1eta*e2s
			sk++

			p0[sk], p1[sk] = q1xi*e3r, q1xi*e3s
			sk++
			p0[sk], p1[sk] = q2xi*e3r, q2xi*e3s
			sk++
			p0[sk], p1[sk] = q3xi*e3r, q3xi*e3s
			sk++
		} else {
			p0[sk], p1[sk] = CRP1(r, s, deriv, e1rf, q1sf), CRP1(r, s, deriv, e1sf, q1sf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e1rf, q2sf), CRP1(r, s, deriv, e1sf, q2sf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e1rf, q3sf), CRP1(r, s, deriv, e1sf, q3sf)
			sk++

			p0[sk], p1[sk] = CRP1(r, s, deriv, e2rf, q3sf), CRP1(r, s, deriv, e2sf, q3sf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e2rf, q2sf), CRP1(r, s, deriv, e2sf, q2sf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e2rf, q1sf), CRP1(r, s, deriv, e2sf, q1sf)
			sk++

			p0[sk], p1[sk] = CRP1(r, s, deriv, e3rf, q1rf), CRP1(r, s, deriv, e3sf, q1rf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e3rf, q2rf), CRP1(r, s, deriv, e3sf, q2rf)
			sk++
			p0[sk], p1[sk] = CRP1(r, s, deriv, e3rf, q3rf), CRP1(r, s, deriv, e3sf, q3rf)
			sk++
		}
	default:
		e1r, e1s := rtb.getCoreBasisTerm(E1, r, s)
		e2r, e2s := rtb.getCoreBasisTerm(E2, r, s)
		e3r, e3s := rtb.getCoreBasisTerm(E3, r, s)
		if deriv == None {
			for j := 0; j < rtb.P+1; j++ {
				leta := rtb.Lagrange1DPoly(s, RGauss, j, SDir)
				p0[sk], p1[sk] = leta*e1r, leta*e1s
				sk++
			}
			for j := 0; j < rtb.P+1; j++ {
				leta := rtb.Lagrange1DPoly(s, RGauss, rtb.P-j, SDir)
				p0[sk], p1[sk] = leta*e2r, leta*e2s
				sk++
			}
			for j := 0; j < rtb.P+1; j++ {
				lxi := rtb.Lagrange1DPoly(r, RGauss, j, RDir)
				p0[sk], p1[sk] = lxi*e3r, lxi*e3s
				sk++
			}
		} else {
			e1Rderiv, e1Sderiv := rtb.getCoreBasisTerm(E1, r, s, deriv)
			e2Rderiv, e2Sderiv := rtb.getCoreBasisTerm(E2, r, s, deriv)
			e3Rderiv, e3Sderiv := rtb.getCoreBasisTerm(E3, r, s, deriv)
			for j := 0; j < rtb.P+1; j++ {
				leta := rtb.Lagrange1DPoly(s, RGauss, j, SDir)
				letaDeriv := rtb.Lagrange1DPoly(s, RGauss, j, SDir, deriv)
				p0[sk], p1[sk] = CRP(e1r, e1Rderiv, leta, letaDeriv), CRP(e1s, e1Sderiv, leta, letaDeriv)
				sk++
			}
			for j := 0; j < rtb.P+1; j++ {
				leta := rtb.Lagrange1DPoly(s, RGauss, rtb.P-j, SDir)
				letaDeriv := rtb.Lagrange1DPoly(s, RGauss, rtb.P-j, SDir, deriv)
				p0[sk], p1[sk] = CRP(e2r, e2Rderiv, leta, letaDeriv), CRP(e2s, e2Sderiv, leta, letaDeriv)
				sk++
			}
			for j := 0; j < rtb.P+1; j++ {
				lxi := rtb.Lagrange1DPoly(r, RGauss, j, RDir)
				lxiDeriv := rtb.Lagrange1DPoly(r, RGauss, j, RDir, deriv)
				p0[sk], p1[sk] = CRP(e3r, e3Rderiv, lxi, lxiDeriv), CRP(e3s, e3Sderiv, lxi, lxiDeriv)
				sk++
			}
		}
	}
	return
}

func (rtb *RTBasis2DSimplex) CalculateBasis() {
	/*
		We follow the basis function construction of V.J. Ervin "Computational Bases for RTk and BDMk on Triangles"
	*/
	var (
		Np     = rtb.Np
		Rd, Sd = rtb.R.DataP, rtb.S.DataP
		p0, p1 []float64
	)
	P := utils.NewMatrix(Np, Np)
	if len(Rd) != Np {
		err := fmt.Errorf("calculated dimension %d doesn't match length of basis nodes %d", Np, len(Rd))
		panic(err)
	}
	// Evaluate at geometric locations
	rowEdge := make([]float64, Np)
	oosr2 := 1. / math.Sqrt(2)
	for ii, rr := range Rd {
		ss := Sd[ii]
		/*
			First, evaluate the polynomial at the (r,s) coordinates
			This is the same set that will be used for all dot products to form the basis matrix
		*/
		p0, p1 = rtb.EvaluateBasisAtLocation(rr, ss)
		// Implement dot product of (unit vector)_ii with each vector term in the polynomial evaluated at location ii
		switch rtb.getLocationType(ii) {
		case E4:
			// Unit vector is [1,0]
			P.M.SetRow(ii, p0)
		case E5:
			// Unit vector is [0,1]
			P.M.SetRow(ii, p1)
		case E1:
			for i := range rowEdge {
				// E3: // Unit vector is [0,-1]
				rowEdge[i] = -p1[i]
			}
			P.M.SetRow(ii, rowEdge)
		case E2:
			for i := range rowEdge {
				// E1: Unit vector is [1/sqrt(2), 1/sqrt(2)]
				rowEdge[i] = oosr2 * (p0[i] + p1[i])
			}
			P.M.SetRow(ii, rowEdge)
		case E3:
			for i := range rowEdge {
				// E2: Unit vector is [-1,0]
				rowEdge[i] = -p0[i]
			}
			P.M.SetRow(ii, rowEdge)
		}
	}
	// Invert [P] = [A] to obtain the coefficients (columns) of polynomials (rows), each row is a polynomial
	A := P.InverseWithCheck()
	// Evaluate 2D polynomial basis at geometric locations, also evaluate derivatives Dr and Ds for Rd and Sd
	P0, P1 := utils.NewMatrix(Np, Np), utils.NewMatrix(Np, Np)
	Pdr0, Pds1 := utils.NewMatrix(Np, Np), utils.NewMatrix(Np, Np)
	for ii, rr := range rtb.R.DataP {
		ss := rtb.S.DataP[ii]
		p0, p1 = rtb.EvaluateBasisAtLocation(rr, ss) // each of p1,p2 stores the polynomial terms for the Rd and Sd directions
		P0.M.SetRow(ii, p0)
		P1.M.SetRow(ii, p1)
		p0, _ = rtb.EvaluateBasisAtLocation(rr, ss, Dr) // each of p0,p1 stores the polynomial terms for the Rd and Sd directions
		_, p1 = rtb.EvaluateBasisAtLocation(rr, ss, Ds) // each of p0,p1 stores the polynomial terms for the Rd and Sd directions
		Pdr0.M.SetRow(ii, p0)                           // Only need the dP/dr(r) term from P(r,s)
		Pds1.M.SetRow(ii, p1)                           // Only need the dP/ds(s) term from P(r,s)
	}
	// Construct the Vandermonde matrices for each direction by multiplying coefficients of constrained basis
	rtb.V[0] = P0.Mul(A)
	rtb.V[1] = P1.Mul(A)
	rtb.Div = Pdr0.Mul(A).Add(Pds1.Mul(A))
	if rtb.P != 0 {
		rtb.DivInt = utils.NewMatrix(rtb.NpInt, Np)
		for i := 0; i < rtb.NpInt; i++ {
			rtb.DivInt.M.SetRow(i, rtb.Div.Row(i).DataP)
		}
	}
	return
}

func (rtb *RTBasis2DSimplex) ExtendGeomToRT(Rint, Sint utils.Vector) (R, S utils.Vector) {
	var (
		N            = rtb.P
		NpEdge       = N + 1
		rData, sData = Rint.DataP, Sint.DataP
		Rd, Sd       []float64
	)
	/*
		Determine geometric locations of edge points, located at Gauss locations in 1D, projected onto the edges
	*/
	GQR := utils.NewVector(N+1, DG1D.LegendreZeros(N))
	/*
		Double the number of interior points to match each direction of the basis
	*/
	if N == 0 { // Special case: when N=0, the interior of the RT element is empty
		Rd, Sd = []float64{}, []float64{}
	} else {
		for i := 0; i < 2; i++ {
			Rd = append(Rd, rData...)
			Sd = append(Sd, sData...)
		}
	}

	// Calculate the triangle edges
	GQRData := GQR.DataP
	rEdgeData := make([]float64, NpEdge*3)
	sEdgeData := make([]float64, NpEdge*3)
	for i := 0; i < NpEdge; i++ {
		gp := GQRData[i]
		// Edge 1
		rEdgeData[i] = gp
		sEdgeData[i] = -1
		// Edge 2 (hypotenuse)
		gpT := 0.5 * (gp + 1)
		rEdgeData[i+NpEdge] = 1 - 2*gpT
		sEdgeData[i+NpEdge] = -1 + 2*gpT
		// Edge 3
		rEdgeData[i+2*NpEdge] = -1
		sEdgeData[i+2*NpEdge] = -gp
	}
	Rd = append(Rd, rEdgeData...)
	Sd = append(Sd, sEdgeData...)
	nn := len(Rd)
	R, S = utils.NewVector(nn, Rd), utils.NewVector(nn, Sd)
	return
}

func (rtb *RTBasis2DSimplex) GetEdgeLocations(F []float64) (Fedge []float64) {
	var (
		Nint     = rtb.NpInt
		NedgeTot = rtb.NpEdge * 3
	)
	Fedge = make([]float64, NedgeTot)
	for i := 0; i < NedgeTot; i++ {
		Fedge[i] = F[i+2*Nint]
	}
	return
}

func (rtb *RTBasis2DSimplex) GetInternalLocations(F []float64) (Finternal []float64) {
	var (
		Nint = rtb.NpInt
	)
	Finternal = make([]float64, Nint)
	for i := 0; i < Nint; i++ {
		Finternal[i] = F[i]
	}
	return
}

/*
Set up the five core vector basis functions used for all order of terms
*/
func (rtb *RTBasis2DSimplex) getCoreBasisTerm(tt RTFunctionNumber, r, s float64,
	derivO ...DerivativeDirection) (p0, p1 float64) {
	var (
		sr2   = math.Sqrt(2)
		oosr2 = 1. / sr2
		deriv DerivativeDirection
	)
	if len(derivO) != 0 {
		deriv = derivO[0]
	} else {
		deriv = None
	}
	switch tt {
	case E1:
		switch deriv {
		case None:
			// p0 = sr2 * xi
			// p1 = sr2 * eta
			p0 = oosr2 * (r + 1)
			p1 = oosr2 * (s + 1)
		case Dr:
			p0 = oosr2
		case Ds:
			p1 = oosr2
		}
	case E2:
		switch deriv {
		case None:
			// p0 = xi - 1
			// p1 = eta
			p0 = 0.5*r - 0.5
			p1 = 0.5*s + 0.5
		case Dr:
			p0 = 0.5
		case Ds:
			p1 = 0.5
		}
	case E3:
		switch deriv {
		case None:
			// p0 = xi
			// p1 = eta - 1
			p0 = 0.5*r + 0.5
			p1 = 0.5*s - 0.5
		case Dr:
			p0 = 0.5
		case Ds:
			p1 = 0.5
		}
	case E4:
		switch deriv {
		case None:
			// p0 = eta * xi
			// p1 = eta * (eta - 1)
			p0 = 0.25 * (r*s + r + s + 1)
			p1 = 0.25 * (s*s - 1)
		case Dr:
			p0 = 0.25 * (s + 1)
		case Ds:
			p0 = 0.25 * (r + 1)
			p1 = 0.5 * s
		}
	case E5:
		switch deriv {
		case None:
			// p0 = xi * (xi - 1)
			// p1 = xi * eta
			p0 = 0.25 * (r*r - 1)
			p1 = 0.25 * (r*s + r + s + 1)
		case Dr:
			p0 = 0.5 * r
			p1 = 0.25 * (s + 1)
		case Ds:
			p1 = 0.25 * (r + 1)
		}
	}
	return
}

type Direction uint8

const (
	RDir Direction = iota
	SDir
)

func (rtb *RTBasis2DSimplex) Lagrange1DPoly(r float64, R []float64, j int,
	dir Direction, derivO ...DerivativeDirection) (p float64) {
	// This is a 1D polynomial, but the nodes are either in the RDir or SDir direction within 2D
	var (
		deriv DerivativeDirection
	)
	if len(derivO) != 0 {
		deriv = derivO[0]
	} else {
		deriv = None
	}
	switch deriv {
	case None:
		var (
			Np1D = len(R)
			XJ   = R[j]
			XI   = r
		)
		if j > Np1D-1 || j < 0 {
			panic("value of j larger than array or less than zero")
		}
		p = 1
		for m, XM := range R {
			if m != j {
				p *= (XI - XM) / (XJ - XM)
			}
		}
	case Dr, Ds:
		switch {
		case dir == SDir && deriv == Dr:
			return
		case dir == RDir && deriv == Ds:
			return
		}
		/*
			From https://en.wikipedia.org/wiki/Lagrange_polynomial#Derivatives
		*/
		var (
			XJ = R[j]
			X  = r
		)
		for i, XI := range R {
			if i != j {
				pp := 1.
				for m, XM := range R {
					if m != i && m != j {
						pp *= (X - XM) / (XJ - XM)
					}
				}
				p += pp / (XJ - XI)
			}
		}
	}
	return
}

func (rtb *RTBasis2DSimplex) LinearPoly(r, s float64, j int, derivO ...DerivativeDirection) (p float64) {
	var (
		deriv DerivativeDirection
	)
	if j > 2 || j < 0 {
		err := fmt.Errorf("linear 2D polynomial called with bad parameters [j]=[%d]", j)
		panic(err)
	}
	if len(derivO) != 0 {
		deriv = derivO[0]
	} else {
		deriv = None
	}
	switch deriv {
	case None:
		switch j {
		case 0:
			p = -0.5 * (r + s)
		case 1:
			p = 0.5*r + 0.5
		case 2:
			p = 0.5*s + 0.5
		}
	case Dr:
		switch j {
		case 0:
			p = -0.5
		case 1:
			p = 0.5
		case 2:
			p = 0
		}
	case Ds:
		switch j {
		case 0:
			p = -0.5
		case 1:
			p = 0
		case 2:
			p = 0.5
		}
	}
	return
}
