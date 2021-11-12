package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

type RTBasis2DSimplex struct {
	P             int       // Polynomial Order
	Np            int       // Number of terms and nodes in basis
	NpInt, NpEdge int       // Number of nodes in interior and on each edge, Np = 2*NpInt + 3*NpEdge
	R, S          []float64 // Node locations within [-1,1] reference triangle
	Scalar2DBasis Basis2D   // Basis used for part of the RT basis construction
}

func NewRTBasis2DSimplex(P int) (rtb *RTBasis2DSimplex) {
	var Rint, Sint utils.Vector
	rtb = &RTBasis2DSimplex{
		P:      P,
		Np:     (P + 1) * (P + 3),
		NpInt:  P * (P + 1) / 2, // Number of interior points is same as the 2D scalar space one order lesser
		NpEdge: P + 1,           // Each edge is P+1 nodes
	}
	if P > 0 {
		Rint, Sint = NodesEpsilon(P - 1)
		rtb.Scalar2DBasis = NewLagrangeBasis2D(P-1, Rint, Sint) // Basis used as part of non-Normal basis elements
	}
	rtb.R, rtb.S = rtb.ExtendGeomToRT(Rint, Sint)
	rtb.CalculateBasis()
	return
}

func convertToUnitTriangle(r, s float64) (xi, eta float64) {
	// Convert R,S coordinates from a [-1,1] triangle to a [0,1] triangle to match up with the basis definition
	xi = (r + 1) / 2
	eta = (s + 1) / 2
	return
}

/*
	Set up the five core vector basis functions used for all order of terms
*/
type TermType uint8

const (
	e1 TermType = iota
	e2
	e3
	e4
	e5
)

func (rtb *RTBasis2DSimplex) getCoreBasisTerm(tt TermType, r, s float64, derivO ...DerivativeDirection) (p0, p1 float64) {
	var (
		xi, eta = convertToUnitTriangle(r, s)
		sr2     = math.Sqrt(2)
		deriv   DerivativeDirection
	)
	if len(derivO) != 0 {
		deriv = derivO[0]
	} else {
		deriv = None
	}
	switch tt {
	case e1:
		switch deriv {
		case None:
			p0 = sr2 * xi
			p1 = sr2 * eta
		case Dr:
			p0 = 0.5 * sr2
		case Ds:
			p1 = 0.5 * sr2
		}
	case e2:
		switch deriv {
		case None:
			p0 = xi - 1
			p1 = eta
		case Dr:
			p0 = 0.5
		case Ds:
			p1 = 0.5
		}
	case e3:
		switch deriv {
		case None:
			p0 = xi
			p1 = eta - 1
		case Dr:
			p0 = 0.5
		case Ds:
			p1 = 0.5
		}
	case e4:
		switch deriv {
		case None:
			p0 = eta * xi
			p1 = eta * (eta - 1)
		case Dr:
			p0 = 0.5 * eta
		case Ds:
			p0 = 0.5 * xi
			p1 = eta - 0.5
		}
	case e5:
		switch deriv {
		case None:
			p0 = xi * (xi - 1)
			p1 = xi * eta
		case Dr:
			p0 = xi - 0.5
			p1 = 0.5 * eta
		case Ds:
			p1 = 0.5 * xi
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
		p = rtb.LagrangePolyAtJ(r, R, j)
	case Dr:
		if dir == RDir {
			return
		}
		p = rtb.LagrangePolyDerivAtJ(r, R, j)
	case Ds:
		if dir == SDir {
			return
		}
		p = rtb.LagrangePolyDerivAtJ(r, R, j)
	}
	return
}

func (rtb *RTBasis2DSimplex) LagrangePolyAtJ(r float64, R []float64, j int) (p float64) {
	var (
		Np1D = len(R)
	)
	if j > Np1D-1 || j < 0 {
		panic("value of j larger than array or less than zero")
	}
	s := 1.
	for i, rBasis := range R {
		if i == j {
			continue
		}
		metric := (r - rBasis) / (R[j] - rBasis)
		s *= metric
	}
	return
}

func (rtb *RTBasis2DSimplex) LagrangePolyDerivAtJ(r float64, R []float64, j int) (dp float64) {
	var (
		p = rtb.LagrangePolyAtJ(r, R, j)
	)
	s := 1.
	for i, rBasis := range R {
		if i == j {
			continue
		}
		s *= 1 / (r - rBasis)
	}
	dp = p * s
	return
}

func (rtb *RTBasis2DSimplex) LinearPoly(r, s float64, i, j int, derivO ...DerivativeDirection) (p float64) {
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
		p = rtb.LinearPolyTerm2D(r, s, i, j)
	case Dr:
		p = rtb.LinearPolyTerm2DDr(i, j)
	case Ds:
		p = rtb.LinearPolyTerm2DDs(i, j)
	}
	return
}

func (rtb *RTBasis2DSimplex) LinearPolyTerm2D(r, s float64, i, j int) (p float64) {
	var (
		xi, eta = convertToUnitTriangle(r, s)
	)
	switch {
	case i == 0 && j == 0:
		p = 1 - xi - eta
	case i == 1 && j == 0:
		p = xi
	case i == 0 && j == 1:
		p = eta
	default:
		err := fmt.Errorf("linear 2D polynomial called with bad parameters [i,j]=[%d,%d]", i, j)
		panic(err)
	}
	return
}

func (rtb *RTBasis2DSimplex) LinearPolyTerm2DDr(i, j int) (dpdr float64) {
	switch {
	case i == 0 && j == 0:
		dpdr = -0.5
	case i == 1 && j == 0:
		dpdr = 0.5
	case i == 0 && j == 1:
		dpdr = 0
	default:
		err := fmt.Errorf("linear 2D polynomial called with bad parameters [i,j]=[%d,%d]", i, j)
		panic(err)
	}
	return
}

func (rtb *RTBasis2DSimplex) LinearPolyTerm2DDs(i, j int) (dpds float64) {
	switch {
	case i == 0 && j == 0:
		dpds = -0.5
	case i == 1 && j == 0:
		dpds = 0
	case i == 0 && j == 1:
		dpds = 0.5
	default:
		err := fmt.Errorf("linear 2D polynomial called with bad parameters [i,j]=[%d,%d]", i, j)
		panic(err)
	}
	return
}

func (rtb *RTBasis2DSimplex) getLocationType(i int) (locationType RTPointType) {
	var (
		NpInt  = rtb.NpInt
		NpEdge = rtb.NpEdge
	)
	switch {
	case i < NpInt:
		// Unit vector is [1,0]
		locationType = InteriorR
	case i >= NpInt && i < 2*NpInt:
		// Unit vector is [0,1]
		locationType = InteriorS
	case i >= 2*NpInt && i < 2*NpInt+NpEdge:
		// Edge1: Unit vector is [0,-1]
		locationType = Edge1
	case i >= 2*NpInt+NpEdge && i < 2*NpInt+2*NpEdge:
		// Edge2: Unit vector is [1/sqrt(2), 1/sqrt(2)]
		locationType = Edge2
	case i >= 2*NpInt+2*NpEdge && i < 2*NpInt+3*NpEdge:
		// Edge3: Unit vector is [-1,0]
		locationType = Edge3
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
	CRP := func(pL, dpL, pR, dpR float64) (dp float64) {
		dp = dpL*pR + pL*dpR
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
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(e4, r, s, deriv)
		sk++
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(e5, r, s, deriv)
		sk++
	case 2:
		// Six "non-Normal" basis functions, use linear 2D polynomial for triangles multiplied by e4 and e5
		e4r, e4s := rtb.getCoreBasisTerm(e4, r, s)
		e5r, e5s := rtb.getCoreBasisTerm(e5, r, s)
		l1 := rtb.LinearPoly(r, s, 0, 0)
		l2 := rtb.LinearPoly(r, s, 1, 0)
		l3 := rtb.LinearPoly(r, s, 0, 1)
		if deriv == None {
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
			l1deriv := rtb.LinearPoly(r, s, 0, 0, deriv)
			l2deriv := rtb.LinearPoly(r, s, 1, 0, deriv)
			l3deriv := rtb.LinearPoly(r, s, 0, 1, deriv)
			e4Rderiv, e4Sderiv := rtb.getCoreBasisTerm(e4, r, s, deriv)
			e5Rderiv, e5Sderiv := rtb.getCoreBasisTerm(e5, r, s, deriv)
			p0[sk], p1[sk] = CRP(e4r, e4Rderiv, l1, l1deriv), CRP(e4s, e4Sderiv, l1, l1deriv)
			sk++
			p0[sk], p1[sk] = CRP(e4r, e4Rderiv, l2, l2deriv), CRP(e4s, e4Sderiv, l2, l2deriv)
			sk++
			p0[sk], p1[sk] = CRP(e4r, e4Rderiv, l3, l3deriv), CRP(e4s, e4Sderiv, l3, l3deriv)
			sk++

			p0[sk], p1[sk] = CRP(e5r, e5Rderiv, l1, l1deriv), CRP(e5s, e5Sderiv, l1, l1deriv)
			sk++
			p0[sk], p1[sk] = CRP(e5r, e5Rderiv, l2, l2deriv), CRP(e5s, e5Sderiv, l2, l2deriv)
			sk++
			p0[sk], p1[sk] = CRP(e5r, e5Rderiv, l3, l3deriv), CRP(e5s, e5Sderiv, l3, l3deriv)
			sk++
		}
	default:
		e4r, e4s := rtb.getCoreBasisTerm(e4, r, s)
		e5r, e5s := rtb.getCoreBasisTerm(e5, r, s)
		var (
			e4Rderiv, e4Sderiv, e5Rderiv, e5Sderiv float64
		)
		if deriv != None {
			e4Rderiv, e4Sderiv = rtb.getCoreBasisTerm(e4, r, s, deriv)
			e5Rderiv, e5Sderiv = rtb.getCoreBasisTerm(e5, r, s, deriv)
		}
		Pint := rtb.P - 1
		for i := 0; i <= Pint; i++ {
			for j := 0; j <= Pint-i; j++ {
				p := rtb.Scalar2DBasis.PolynomialTerm(r, s, i, j)
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
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(e1, r, s, deriv)
		sk++
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(e2, r, s, deriv)
		sk++
		p0[sk], p1[sk] = rtb.getCoreBasisTerm(e3, r, s, deriv)
		sk++
	case 1:
		e1r, e1s := rtb.getCoreBasisTerm(e1, r, s)
		e2r, e2s := rtb.getCoreBasisTerm(e2, r, s)
		e3r, e3s := rtb.getCoreBasisTerm(e3, r, s)
		l1xi := rtb.Lagrange1DPoly(r, RGauss, 0, RDir)
		l1eta := rtb.Lagrange1DPoly(s, RGauss, 0, SDir)
		l2xi := rtb.Lagrange1DPoly(r, RGauss, 1, RDir)
		l2eta := rtb.Lagrange1DPoly(s, RGauss, 1, SDir)
		if deriv == None {
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
			e1Rderiv, e1Sderiv := rtb.getCoreBasisTerm(e1, r, s, deriv)
			e2Rderiv, e2Sderiv := rtb.getCoreBasisTerm(e2, r, s, deriv)
			e3Rderiv, e3Sderiv := rtb.getCoreBasisTerm(e3, r, s, deriv)
			l1xiDeriv := rtb.Lagrange1DPoly(r, RGauss, 0, RDir, deriv)
			l1etaDeriv := rtb.Lagrange1DPoly(s, RGauss, 0, SDir, deriv)
			l2xiDeriv := rtb.Lagrange1DPoly(r, RGauss, 1, RDir, deriv)
			l2etaDeriv := rtb.Lagrange1DPoly(s, RGauss, 1, SDir, deriv)
			p0[sk], p1[sk] = CRP(e1r, e1Rderiv, l1eta, l1etaDeriv), CRP(e1s, e1Sderiv, l1eta, l1etaDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e1r, e1Rderiv, l2eta, l2etaDeriv), CRP(e1s, e1Sderiv, l2eta, l2etaDeriv)
			sk++

			p0[sk], p1[sk] = CRP(e2r, e2Rderiv, l2eta, l2etaDeriv), CRP(e2s, e2Sderiv, l2eta, l2etaDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e2r, e2Rderiv, l1eta, l1etaDeriv), CRP(e2s, e2Sderiv, l1eta, l1etaDeriv)
			sk++

			p0[sk], p1[sk] = CRP(e3r, e3Rderiv, l1xi, l1xiDeriv), CRP(e3s, e3Sderiv, l1xi, l1xiDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e3r, e3Rderiv, l2xi, l2xiDeriv), CRP(e3s, e3Sderiv, l2xi, l2xiDeriv)
			sk++
		}
	case 2:
		e1r, e1s := rtb.getCoreBasisTerm(e1, r, s)
		e2r, e2s := rtb.getCoreBasisTerm(e2, r, s)
		e3r, e3s := rtb.getCoreBasisTerm(e3, r, s)
		l1xi := rtb.Lagrange1DPoly(r, RGauss, 0, RDir)
		l1eta := rtb.Lagrange1DPoly(s, RGauss, 0, SDir)
		l2xi := rtb.Lagrange1DPoly(r, RGauss, 1, RDir)
		l2eta := rtb.Lagrange1DPoly(s, RGauss, 1, SDir)
		l3xi := rtb.Lagrange1DPoly(r, RGauss, 2, RDir)
		l3eta := rtb.Lagrange1DPoly(s, RGauss, 2, SDir)
		if deriv == None {
			p0[sk], p1[sk] = l1eta*e1r, l1eta*e1s
			sk++
			p0[sk], p1[sk] = l2eta*e1r, l2eta*e1s
			sk++
			p0[sk], p1[sk] = l3eta*e1r, l3eta*e1s
			sk++
			p0[sk], p1[sk] = l3eta*e2r, l3eta*e2s
			sk++
			p0[sk], p1[sk] = l2eta*e2r, l2eta*e2s
			sk++
			p0[sk], p1[sk] = l1eta*e2r, l1eta*e2s
			sk++
			p0[sk], p1[sk] = l1xi*e3r, l1xi*e3s
			sk++
			p0[sk], p1[sk] = l2xi*e3r, l2xi*e3s
			sk++
			p0[sk], p1[sk] = l3xi*e3r, l3xi*e3s
			sk++
		} else {
			e1Rderiv, e1Sderiv := rtb.getCoreBasisTerm(e1, r, s, deriv)
			e2Rderiv, e2Sderiv := rtb.getCoreBasisTerm(e2, r, s, deriv)
			e3Rderiv, e3Sderiv := rtb.getCoreBasisTerm(e3, r, s, deriv)
			l1xiDeriv := rtb.Lagrange1DPoly(r, RGauss, 0, RDir, deriv)
			l1etaDeriv := rtb.Lagrange1DPoly(s, RGauss, 0, SDir, deriv)
			l2xiDeriv := rtb.Lagrange1DPoly(r, RGauss, 1, RDir, deriv)
			l2etaDeriv := rtb.Lagrange1DPoly(s, RGauss, 1, SDir, deriv)
			l3xiDeriv := rtb.Lagrange1DPoly(r, RGauss, 2, RDir, deriv)
			l3etaDeriv := rtb.Lagrange1DPoly(s, RGauss, 2, SDir, deriv)
			p0[sk], p1[sk] = CRP(e1r, e1Rderiv, l1eta, l1etaDeriv), CRP(e1s, e1Sderiv, l1eta, l1etaDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e1r, e1Rderiv, l2eta, l2etaDeriv), CRP(e1s, e1Sderiv, l2eta, l2etaDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e1r, e1Rderiv, l3eta, l3etaDeriv), CRP(e1s, e1Sderiv, l3eta, l3etaDeriv)
			sk++

			p0[sk], p1[sk] = CRP(e2r, e2Rderiv, l3eta, l3etaDeriv), CRP(e2s, e2Sderiv, l3eta, l3etaDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e2r, e2Rderiv, l2eta, l2etaDeriv), CRP(e2s, e2Sderiv, l2eta, l2etaDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e2r, e2Rderiv, l1eta, l1etaDeriv), CRP(e2s, e2Sderiv, l1eta, l1etaDeriv)
			sk++

			p0[sk], p1[sk] = CRP(e3r, e3Rderiv, l1xi, l1xiDeriv), CRP(e3s, e3Sderiv, l1xi, l1xiDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e3r, e3Rderiv, l2xi, l2xiDeriv), CRP(e3s, e3Sderiv, l2xi, l2xiDeriv)
			sk++
			p0[sk], p1[sk] = CRP(e3r, e3Rderiv, l3xi, l3xiDeriv), CRP(e3s, e3Sderiv, l3xi, l3xiDeriv)
			sk++
		}
	default:
	}
	return
}

func (rtb *RTBasis2DSimplex) CalculateBasis() (P utils.Matrix) {
	/*
		We follow the basis function construction of V.J. Ervin "Computational Bases for RTk and BDMk on Triangles"
	*/
	var (
		Np = rtb.Np
		b1 = rtb.Scalar2DBasis
		//N2DBasis = rtb.NpInt
		R, S   = rtb.R, rtb.S
		deriv  = None
		tFunc  func(r, s float64, i, j int) (val float64)
		p0, p1 []float64
	)
	if len(R) != Np {
		err := fmt.Errorf("calculated dimension %d doesn't match length of basis nodes %d", Np, len(R))
		panic(err)
	}
	P = utils.NewMatrix(Np, Np)
	if rtb.P > 0 {
		switch deriv {
		case None:
			tFunc = b1.PolynomialTerm
		case Dr:
			tFunc = b1.PolynomialTermDr
		case Ds:
			tFunc = b1.PolynomialTermDs
		}
	}
	_ = tFunc
	// Evaluate at geometric locations
	P = utils.NewMatrix(Np, Np)
	rowEdge := make([]float64, Np)
	oosr2 := math.Sqrt(2)
	for ii, rr := range R {
		ss := S[ii]
		/*
			First, evaluate the polynomial at the (r,s) coordinates
			This is the same set that will be used for all dot products to form the basis matrix
		*/
		p0, p1 = rtb.EvaluateBasisAtLocation(rr, ss)
		// Implement dot product of (unit vector)_ii with each vector term in the polynomial evaluated at location ii
		switch rtb.getLocationType(ii) {
		case InteriorR:
			// Unit vector is [1,0]
			P.M.SetRow(ii, p0)
		case InteriorS:
			// Unit vector is [0,1]
			P.M.SetRow(ii, p1)
		case Edge1:
			for i := range rowEdge {
				// Edge3: // Unit vector is [0,-1]
				rowEdge[i] = -p1[i]
			}
			P.M.SetRow(ii, rowEdge)
		case Edge2:
			for i := range rowEdge {
				// Edge1: Unit vector is [1/sqrt(2), 1/sqrt(2)]
				rowEdge[i] = oosr2 * (p0[i] + p1[i])
			}
			P.M.SetRow(ii, rowEdge)
		case Edge3:
			for i := range rowEdge {
				// Edge2: Unit vector is [-1,0]
				rowEdge[i] = -p0[i]
			}
			P.M.SetRow(ii, rowEdge)
		}
	}
	return
}

func (rtb *RTBasis2DSimplex) ExtendGeomToRT(Rint, Sint utils.Vector) (R, S []float64) {
	var (
		N            = rtb.P
		NpEdge       = N + 1
		rData, sData = Rint.DataP, Sint.DataP
	)
	/*
		Determine geometric locations of edge points, located at Gauss locations in 1D, projected onto the edges
	*/
	GQR := utils.NewVector(N+1, DG1D.LegendreZeros(N))
	/*
		Double the number of interior points to match each direction of the basis
	*/
	if N == 0 { // Special case: when N=0, the interior of the RT element is empty
		R, S = []float64{}, []float64{}
	} else {
		for i := 0; i < 2; i++ {
			R = append(R, rData...)
			S = append(S, sData...)
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
	R = append(R, rEdgeData...)
	S = append(S, sEdgeData...)
	return
}
