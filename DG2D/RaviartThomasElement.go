package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

/*
	Definition:
		Raviart Thomas element
		Total of (N+1)(N+3) degrees of freedom. (N)(N+1) interior degrees of
		freedom split between R and S directions. 3*(N+1) edge degrees of
		freedom, (N+1) per each of 3 edges on the triangle.

	Inputs:
		(N)(N+2)/2 [r,s] points from the interior of the [-1,1] triangle

	Outputs:
		[R,S]: Coordinates of points within element

		First (N)(N+1)/2 points: Interior points, excluding edges in
			[-1,1] triangle coordinates corresponding to the R
			direction DOFs
		Next (N)(N+1)/2 points: Same as above for the S direction
		Next (N+1) points: Edge 1 locations
		Next (N+1) points: Edge 2 (Hypotenuse) locations
		Next (N+1) points: Edge 3 locations


	Methods:
		The evaluation matrix rt.A relates the coefficients of the RT element
		basis to the solution vectors. For the RT element at degree 1 (RT1),
		this equation is:

		⎡ φ₁(P₁)  φ₂(P₁)  φ₃(P₁)  φ₄(P₁)  φ₅(P₁)  φ₆(P₁)  φ₇(P₁)  φ₈(P₁) ⎤   ⎡ c₁ ⎤   ⎡ f₁ ⎤
		⎢ φ₁(P₁)  φ₂(P₁)  φ₃(P₁)  φ₄(P₁)  φ₅(P₁)  φ₆(P₁)  φ₇(P₁)  φ₈(P₁) ⎥   ⎢ c₂ ⎥   ⎢ f₂ ⎥
		⎢ φ₁(P₂)  φ₂(P₂)  φ₃(P₂)  φ₄(P₂)  φ₅(P₂)  φ₆(P₂)  φ₇(P₂)  φ₈(P₂) ⎥   ⎢ c₃ ⎥   ⎢ f₃ ⎥
		⎢ φ₁(P₃)  φ₂(P₃)  φ₃(P₃)  φ₄(P₃)  φ₅(P₃)  φ₆(P₃)  φ₇(P₃)  φ₈(P₃) ⎥   ⎢ c₄ ⎥   ⎢ f₄ ⎥
		⎢ φ₁(P₄)  φ₂(P₄)  φ₃(P₄)  φ₄(P₄)  φ₅(P₄)  φ₆(P₄)  φ₇(P₄)  φ₈(P₄) ⎥ * ⎢ c₅ ⎥ = ⎢ f₅ ⎥
		⎢ φ₁(P₅)  φ₂(P₅)  φ₃(P₅)  φ₄(P₅)  φ₅(P₅)  φ₆(P₅)  φ₇(P₅)  φ₈(P₅) ⎥   ⎢ c₆ ⎥   ⎢ f₆ ⎥
		⎢ φ₁(P₆)  φ₂(P₆)  φ₃(P₆)  φ₄(P₆)  φ₅(P₆)  φ₆(P₆)  φ₇(P₆)  φ₈(P₆) ⎥   ⎢ c₇ ⎥   ⎢ f₇ ⎥
		⎣ φ₁(P₇)  φ₂(P₇)  φ₃(P₇)  φ₄(P₇)  φ₅(P₇)  φ₆(P₇)  φ₇(P₇)  φ₈(P₇) ⎦   ⎢ c₈ ⎥   ⎢ f₈ ⎥

		Notes for RT1:
		1) there is one internal position in [R,S] shared by the two interior DOFs,
		one for each of R and S directions.
		2) there are 6 positions, one for each of the DOFs on the edges, 2
		per edge.
		3) each basis function φ is evaluated fully at each position. This
		amounts to summation of the full polynomial range of each function
		at the corresponding position.

		Note that (3) above is in contrast to the Vendermonde matrix, which
		is very different in that each column in a Vandermonde matrix is an
		evaluation of a single polynomial term. The Vandermonde matrix
		relates the coefficients of the polynomial terms to the function
		representation in a polynomial, which is NOT what we are doing here.

		For an RT element, the function is distinct on the interior and the
		edges. There is no value of the edge basis functions in the
		interior, and there is no value of the interior functions on the
		edge. We correct the above evaluation matrix for RT1 to reflect the
		nature of the RT element:
*/
/*
   ⎡ φ₁(P₁)  φ₂(P₁)    0       0       0       0       0       0    ⎤   ⎡ c₁ ⎤   ⎡ f₁ ⎤
   ⎢ φ₁(P₁)  φ₂(P₁)    0       0       0       0       0       0    ⎥   ⎢ c₂ ⎥   ⎢ f₂ ⎥
   ⎢   0       0     φ₃(P₂)  φ₄(P₂)    0       0       0       0    ⎥   ⎢ c₃ ⎥   ⎢ f₃ ⎥
   ⎢   0       0     φ₃(P₃)  φ₄(P₃)    0       0       0       0    ⎥   ⎢ c₄ ⎥   ⎢ f₄ ⎥
   ⎢   0       0       0       0     φ₅(P₄)  φ₆(P₄)    0       0    ⎥   ⎢ c₅ ⎥ = ⎢ f₅ ⎥
   ⎢   0       0       0       0     φ₅(P₅)  φ₆(P₅)    0       0    ⎥   ⎢ c₆ ⎥   ⎢ f₆ ⎥
   ⎢   0       0       0       0       0       0     φ₇(P₆)  φ₈(P₆) ⎥   ⎢ c₇ ⎥   ⎢ f₇ ⎥
   ⎢   0       0       0       0       0       0     φ₇(P₇)  φ₈(P₇) ⎥   ⎢ c₈ ⎥   ⎢ f₈ ⎥
*/
/* Here is the RT2 evaluation matrix to make the pattern clear.
Note for RT2:
1) There are 3 internal positions in [R,S] (N=2, (N)*(N+1)/2 = 3
2) There are 3 edge DOFs for each edge, total of 9 positions
3) The edge basis functions are only evaluated on each edge and the interior
basis functions are only evaluated at each interior position.
*/
/*
   ⎡ φ₁(P₁)  φ₂(P₁)  φ₃(P₁)  φ₄(P₁)  φ₅(P₁)  φ₆(P₁)    0       0       0       0       0       0    0        0        0     ⎤   ⎡ c₁  ⎤   ⎡ f₁  ⎤
   ⎢ φ₁(P₂)  φ₂(P₂)  φ₃(P₂)  φ₄(P₂)  φ₅(P₂)  φ₆(P₂)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₂  ⎥   ⎢ f₂  ⎥
   ⎢ φ₁(P₃)  φ₂(P₃)  φ₃(P₃)  φ₄(P₃)  φ₅(P₃)  φ₆(P₃)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₃  ⎥   ⎢ f₃  ⎥
   ⎢ φ₁(P₁)  φ₂(P₁)  φ₃(P₁)  φ₄(P₁)  φ₅(P₁)  φ₆(P₁)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₄  ⎥   ⎢ f₄  ⎥
   ⎢ φ₁(P₂)  φ₂(P₂)  φ₃(P₂)  φ₄(P₂)  φ₅(P₂)  φ₆(P₂)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₅  ⎥   ⎢ f₅  ⎥
   ⎢ φ₁(P₃)  φ₂(P₃)  φ₃(P₃)  φ₄(P₃)  φ₅(P₃)  φ₆(P₃)    0       0       0       0       0       0    0        0        0     ⎥   ⎢ c₆  ⎥   ⎢ f₆  ⎥
   ⎢   0       0       0       0       0       0     φ₇(P₄)  φ₈(P₄)  φ₉(P₄)    0       0       0    0        0        0     ⎥   ⎢ c₇  ⎥   ⎢ f₇  ⎥
   ⎢   0       0       0       0       0       0     φ₇(P₅)  φ₈(P₅)  φ₉(P₅)    0       0       0    0        0        0     ⎥   ⎢ c₈  ⎥ = ⎢ f₈  ⎥
   ⎢   0       0       0       0       0       0     φ₇(P₆)  φ₈(P₆)  φ₉(P₆)    0       0       0    0        0        0     ⎥   ⎢ c₉  ⎥   ⎢ f₉  ⎥
   ⎢   0       0       0       0       0       0       0       0       0 φ₁₀(P₇) φ₁₁(P₇) φ₁₂(P₇)    0        0        0     ⎥   ⎢ c₁₀ ⎥   ⎢ f₁₀ ⎥
   ⎢   0       0       0       0       0       0       0       0       0 φ₁₀(P₈) φ₁₁(P₈) φ₁₂(P₈)    0        0        0     ⎥   ⎢ c₁₁ ⎥   ⎢ f₁₁ ⎥
   ⎢   0       0       0       0       0       0       0       0       0 φ₁₀(P₉) φ₁₁(P₉) φ₁₂(P₉)    0        0        0     ⎥   ⎢ c₁₂ ⎥   ⎢ f₁₂ ⎥
   ⎢   0       0       0       0       0       0       0       0       0       0       0       0  φ₁₃(P₁₀) φ₁₄(P₁₀) φ₁₅(P₁₀)⎥   ⎢ c₁₃ ⎥   ⎢ f₁₃ ⎥
   ⎢   0       0       0       0       0       0       0       0       0       0       0       0  φ₁₃(P₁₁) φ₁₄(P₁₁) φ₁₅(P₁₁)⎥   ⎢ c₁₄ ⎥   ⎢ f₁₄ ⎥
   ⎢   0       0       0       0       0       0       0       0       0       0       0       0  φ₁₃(P₁₂) φ₁₄(P₁₂) φ₁₅(P₁₂)⎥   ⎢ c₁₅ ⎥   ⎢ f₁₅ ⎥
*/

type DerivativeDirection uint8

const (
	None DerivativeDirection = iota
	Dr
	Ds
)

type RTElement struct {
	P                                int             // Order of element
	Np                               int             // Number of points in element
	NpEdge, NpInt                    int             // Number of Edge and Interior points
	A                                utils.Matrix    // Polynomial coefficient matrix, NpxNp
	V                                [2]utils.Matrix // Vandermonde matrix for each direction r and s, [2]xNpxNp
	Div, DivInt                      utils.Matrix    // Divergence matrix, NpxNp for all, NintxNp Interior Points
	R, S                             utils.Vector    // Point locations defining element in [-1,1] Triangle, NpxNp
	RTPolyBasis2D_A, RTPolyBasis2D_B *JacobiBasis2D
	RTPolyBasis1D_Edge1              []float64 // Edge1 polynomial
	RTPolyBasis1D_Edge2              []float64 // Edge2 polynomial
	RTPolyBasis1D_Edge3              []float64 // Edge3 polynomial
}

func NewRTElement(P int) (rt *RTElement) {
	// We expect that there are points in R and S to match the dimension of dim(P(NFlux-1))
	/*
		<---- NpInt ----><---- NpInt ----><---NpEdge----><---NpEdge----><---NpEdge---->
		         Solution Points          Edge 1 pts	Edge 2 pts	  Edge 3 pts
		<---- NpInt ----><---- NpInt ----><---NpEdge----><---NpEdge----><---NpEdge---->
	*/
	var (
		RInt, SInt utils.Vector
	)
	rt = &RTElement{
		P:      P,
		Np:     (P + 1) * (P + 3),
		NpInt:  P * (P + 1) / 2, // Number of interior points is same as the 2D scalar space one order lesser
		NpEdge: P + 1,           // Each edge is P+1 nodes
	}
	if P > 0 {
		if P < 9 {
			RInt, SInt = NodesEpsilon(P - 1)
		} else {
			RInt, SInt = XYtoRS(Nodes2D(P - 1))
		}
	}
	rt.R, rt.S = rt.ExtendGeomToRT(RInt, SInt)
	rt.RTPolyBasis2D_A = NewJacobiBasis2D(rt.P-1, RInt, SInt, 1, 0)
	rt.RTPolyBasis2D_B = NewJacobiBasis2D(rt.P-1, RInt, SInt, 0, 1)
	// For edges 1 and 3, we can use the same definition on R, as S is the same
	rt.RTPolyBasis1D_Edge1 = DG1D.JacobiP(rt.getEdgeCoordinates(1), 0, 0, rt.P)
	rt.RTPolyBasis1D_Edge2 = DG1D.JacobiP(rt.getEdgeCoordinates(2), 0, 0, rt.P)
	rt.RTPolyBasis1D_Edge3 = DG1D.JacobiP(rt.getEdgeCoordinates(3), 0, 0, rt.P)
	rt.CalculateBasis()
	return
}

func (rt *RTElement) getEdgeCoordinates(edgeNum int) (SS utils.Vector) {
	// Edge 1 is S=-1 (bottom of tri)
	// Edge 2 is Hypotenuse
	// Edge 3 is R=-1 (left side of tri)
	switch edgeNum {
	case 1:
		SS = rt.R.Subset(2*rt.NpInt, 2*rt.NpInt+rt.NpEdge-1)
	case 2:
		SS = rt.R.Subset(2*rt.NpInt+rt.NpEdge, 2*rt.NpInt+2*rt.NpEdge-1)
	case 3:
		SS = rt.S.Subset(2*rt.NpInt+2*rt.NpEdge, 2*rt.NpInt+3*rt.NpEdge-1)
	default:
		panic("invalid edgeNum")
	}
	return
}

func (rt *RTElement) ExtendGeomToRT(Rint, Sint utils.Vector) (R, S utils.Vector) {
	var (
		N            = rt.P
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

func (rt *RTElement) EvaluateRTBasisFunction(functionNumber int, r, s float64,
	derivO ...DerivativeDirection) (psi float64) {
	// This function evaluates all of the terms of a poly for one basis location
	switch rt.getLocationType(functionNumber) {
	case InteriorR:
		psi = rt.RTPolyBasis2D_A.GetPolynomialEvaluation(r, s, derivO...)
	case InteriorS:
		psi = rt.RTPolyBasis2D_B.GetPolynomialEvaluation(r, s, derivO...)
	case Edge1:
		// TODO:
	}
	return
}

func (rt *RTElement) CalculateBasis() {
	// Invert [P] = [A] to obtain the coefficients (columns) of polynomials (rows), each row is a polynomial
	// rt.A = P.InverseWithCheck()
	return
}

func NodesEpsilon(N int) (R, S utils.Vector) {
	/*
		From the 2017 paper "A Direct Flux Reconstruction Scheme for Advection Diffusion Problems on Triangular Grids"

		This is a node set that is compatible with DFR in that it implements colocated solution and flux points for the
		interior nodes, while enabling a set of face nodes for the N+1 degree flux polynomial

		There are two node sets, one for N=3 and one for N=4. They were computed via an optimization, and are only
		available for N=3 and N=4. Also, the convergence of N=3 is degraded for diffusion problems.

		Therefore, only the N=4 points should be used for Viscous solutions, while the N=3 nodes are fine for inviscid
	*/
	var (
		Np       = (N + 1) * (N + 2) / 2
		epsD     []float64
		types    = []string{"Linear", "Williams and Shun", "Romero and Jameson(2017)"}
		nameType string
	)
	switch N {
	// Cases 3,4 from Romero and Jameson, Others from Williams and Shun
	case 0:
		epsD = []float64{
			0.3333333333333333,
			0.3333333333333333,
			0.3333333333333333,
		}
		nameType = types[0]
	case 1:
		epsD = []float64{
			0.666666666666667, 0.166666666666667, 0.166666666666667,
			0.166666666666667, 0.666666666666667, 0.166666666666667,
			0.166666666666667, 0.166666666666667, 0.666666666666667,
		}
		nameType = types[0]
	case 2:
		epsD = []float64{
			0.816847572980440, 0.091576213509780, 0.091576213509780, 0.445948490915964, 0.445948490915964, 0.108103018168071,
			0.091576213509780, 0.816847572980440, 0.091576213509780, 0.445948490915964, 0.108103018168071, 0.445948490915964,
			0.091576213509780, 0.091576213509780, 0.816847572980440, 0.108103018168071, 0.445948490915964, 0.445948490915964,
		}
		nameType = types[1]
	case 3:
		epsD = []float64{
			0.3333333333333333, 0.055758983558155, 0.88848203288369, 0.055758983558155, 0.290285227512689, 0.6388573870878149, 0.290285227512689, 0.6388573870878149, 0.070857385399496, 0.070857385399496,
			0.3333333333333333, 0.055758983558155, 0.055758983558155, 0.88848203288369, 0.070857385399496, 0.290285227512689, 0.6388573870878149, 0.070857385399496, 0.290285227512689, 0.6388573870878149,
			0.3333333333333333, 0.88848203288369, 0.055758983558155, 0.055758983558155, 0.6388573870878149, 0.070857385399496, 0.070857385399496, 0.290285227512689, 0.6388573870878149, 0.290285227512689,
		}
		nameType = types[2]
	case 4:
		epsD = []float64{
			0.034681580220044, 0.9306368395599121, 0.034681580220044, 0.243071555674492, 0.513856888651016, 0.243071555674492, 0.473372556704605, 0.05325488659079003, 0.473372556704605, 0.200039998995093, 0.752666332493468, 0.200039998995093, 0.752666332493468, 0.047293668511439, 0.047293668511439,
			0.034681580220044, 0.034681580220044, 0.9306368395599121, 0.243071555674492, 0.243071555674492, 0.513856888651016, 0.473372556704605, 0.473372556704605, 0.05325488659079003, 0.047293668511439, 0.200039998995093, 0.752666332493468, 0.047293668511439, 0.200039998995093, 0.752666332493468,
			0.9306368395599121, 0.034681580220044, 0.034681580220044, 0.513856888651016, 0.243071555674492, 0.243071555674492, 0.05325488659079003, 0.473372556704605, 0.473372556704605, 0.752666332493468, 0.047293668511439, 0.047293668511439, 0.200039998995093, 0.752666332493468, 0.200039998995093,
		}
		nameType = types[2]
	case 5:
		epsD = []float64{
			0.943774095634672, 0.028112952182664, 0.028112952182664, 0.645721803061365, 0.177139098469317, 0.177139098469317, 0.405508595867433, 0.405508595867433, 0.188982808265134, 0.148565812270887, 0.148565812270887, 0.033533207700614, 0.817900980028499, 0.817900980028499, 0.033533207700614, 0.357196298615681, 0.357196298615681, 0.037824789609186, 0.604978911775132, 0.604978911775132, 0.037824789609186,
			0.028112952182664, 0.943774095634672, 0.028112952182664, 0.177139098469317, 0.645721803061365, 0.177139098469317, 0.405508595867433, 0.188982808265134, 0.405508595867433, 0.817900980028499, 0.033533207700614, 0.148565812270887, 0.148565812270887, 0.033533207700614, 0.817900980028499, 0.604978911775132, 0.037824789609186, 0.357196298615681, 0.357196298615681, 0.037824789609186, 0.604978911775132,
			0.028112952182664, 0.028112952182664, 0.943774095634672, 0.177139098469317, 0.177139098469317, 0.645721803061365, 0.188982808265134, 0.405508595867433, 0.405508595867433, 0.033533207700614, 0.817900980028499, 0.817900980028499, 0.033533207700614, 0.148565812270887, 0.148565812270887, 0.037824789609186, 0.604978911775132, 0.604978911775132, 0.037824789609186, 0.357196298615681, 0.357196298615681,
		}
		nameType = types[1]
	case 6:
		epsD = []float64{
			0.960045625755613, 0.019977187122193, 0.019977187122193, 0.736556464940005, 0.131721767529998, 0.131721767529998, 0.333333333333333, 0.485135346793461, 0.485135346793461, 0.029729306413079, 0.107951981846011, 0.107951981846011, 0.024136808036039, 0.867911210117951, 0.867911210117951, 0.024136808036039, 0.270840772921567, 0.270840772921567, 0.028286656697710, 0.700872570380723, 0.700872570380723, 0.028286656697710, 0.316549598844617, 0.316549598844617, 0.146795716949245, 0.536654684206138, 0.536654684206138, 0.146795716949245,
			0.019977187122193, 0.960045625755613, 0.019977187122193, 0.131721767529998, 0.736556464940005, 0.131721767529998, 0.333333333333333, 0.485135346793461, 0.029729306413079, 0.485135346793461, 0.867911210117951, 0.024136808036039, 0.107951981846011, 0.107951981846011, 0.024136808036039, 0.867911210117951, 0.700872570380723, 0.028286656697710, 0.270840772921567, 0.270840772921567, 0.028286656697710, 0.700872570380723, 0.536654684206138, 0.146795716949245, 0.316549598844617, 0.316549598844617, 0.146795716949245, 0.536654684206138,
			0.019977187122193, 0.019977187122193, 0.960045625755613, 0.131721767529998, 0.131721767529998, 0.736556464940005, 0.333333333333333, 0.029729306413079, 0.485135346793461, 0.485135346793461, 0.024136808036039, 0.867911210117951, 0.867911210117951, 0.024136808036039, 0.107951981846011, 0.107951981846011, 0.028286656697710, 0.700872570380723, 0.700872570380723, 0.028286656697710, 0.270840772921567, 0.270840772921567, 0.146795716949245, 0.536654684206138, 0.536654684206138, 0.146795716949245, 0.316549598844617, 0.316549598844617,
		}
		nameType = types[1]
	case 7:
		epsD = []float64{
			0.957657154441070, 0.021171422779465, 0.021171422779465, 0.798831205208225, 0.100584397395888, 0.100584397395888, 0.457923384576135, 0.271038307711932, 0.271038307711932, 0.440191258403832, 0.440191258403832, 0.119617483192335, 0.101763679498021, 0.101763679498021, 0.018256679074748, 0.879979641427232, 0.879979641427232, 0.018256679074748, 0.394033271669987, 0.394033271669987, 0.023404705466341, 0.582562022863673, 0.582562022863673, 0.023404705466341, 0.226245530909229, 0.226245530909229, 0.022223854547989, 0.751530614542782, 0.751530614542782, 0.022223854547989, 0.635737183263105, 0.635737183263105, 0.115183589115563, 0.249079227621332, 0.249079227621332, 0.115183589115563,
			0.021171422779465, 0.957657154441070, 0.021171422779465, 0.100584397395888, 0.798831205208225, 0.100584397395888, 0.271038307711932, 0.457923384576135, 0.271038307711932, 0.440191258403832, 0.119617483192335, 0.440191258403832, 0.879979641427232, 0.018256679074748, 0.101763679498021, 0.101763679498021, 0.018256679074748, 0.879979641427232, 0.582562022863673, 0.023404705466341, 0.394033271669987, 0.394033271669987, 0.023404705466341, 0.582562022863673, 0.751530614542782, 0.022223854547989, 0.226245530909229, 0.226245530909229, 0.022223854547989, 0.751530614542782, 0.249079227621332, 0.115183589115563, 0.635737183263105, 0.635737183263105, 0.115183589115563, 0.249079227621332,
			0.021171422779465, 0.021171422779465, 0.957657154441070, 0.100584397395888, 0.100584397395888, 0.798831205208225, 0.271038307711932, 0.271038307711932, 0.457923384576135, 0.119617483192335, 0.440191258403832, 0.440191258403832, 0.018256679074748, 0.879979641427232, 0.879979641427232, 0.018256679074748, 0.101763679498021, 0.101763679498021, 0.023404705466341, 0.582562022863673, 0.582562022863673, 0.023404705466341, 0.394033271669987, 0.394033271669987, 0.022223854547989, 0.751530614542782, 0.751530614542782, 0.022223854547989, 0.226245530909229, 0.226245530909229, 0.115183589115563, 0.249079227621332, 0.249079227621332, 0.115183589115563, 0.635737183263105, 0.635737183263105,
		}
		nameType = types[1]
	default:
		panic(fmt.Errorf("Epsilon nodes not defined for N = %v\n", N))
	}
	fmt.Printf("Node distribution type is: [%s]\n", nameType)
	eps := utils.NewMatrix(3, Np, epsD)
	T := utils.NewMatrix(2, 3, []float64{
		-1, 1, -1,
		-1, -1, 1,
	})
	RS := T.Mul(eps)
	R = RS.Row(0)
	S = RS.Row(1)
	return
}

func (rt *RTElement) GetInternalLocations(F []float64) (Finternal []float64) {
	var (
		Nint = rt.NpInt
	)
	Finternal = make([]float64, Nint)
	for i := 0; i < Nint; i++ {
		Finternal[i] = F[i]
	}
	return
}

func (rt *RTElement) GetEdgeLocations(F []float64) (Fedge []float64) {
	var (
		Nint     = rt.NpInt
		NedgeTot = rt.NpEdge * 3
	)
	Fedge = make([]float64, NedgeTot)
	for i := 0; i < NedgeTot; i++ {
		Fedge[i] = F[i+2*Nint]
	}
	return
}

func (rt *RTElement) ProjectFunctionOntoDOF(s1, s2 []float64) (sp []float64) {
	// For each location in {R,S}, project the input vector function [s1,s2]
	// on to the degrees of freedom of the element (not the basis)
	var (
		Np = len(s1)
	)
	sp = make([]float64, Np)
	oosr2 := 1 / math.Sqrt(2)
	for i := range s1 {
		switch rt.getLocationType(i) {
		case InteriorR:
			// Unit vector is [1,0]
			sp[i] = s1[i]
		case InteriorS:
			// Unit vector is [0,1]
			sp[i] = s2[i]
		case Edge1:
			// Edge1: // Unit vector is [0,-1]
			sp[i] = -s2[i]
		case Edge2:
			// Edge2: Unit vector is [1/sqrt(2), 1/sqrt(2)]
			sp[i] = (s1[i] + s2[i]) * oosr2
		case Edge3:
			// Edge3: Unit vector is [-1,0]
			sp[i] = -s1[i]
		}
	}
	return
}

type RTPointType uint8

const (
	All       RTPointType = iota
	InteriorR             // R component of vector field
	InteriorS             // S component of vector field
	Edge1                 // Edge from vertex 0-1, Normal is [0,-1]
	Edge2                 // Edge from vertex 1-2, Normal is [1./sqrt(2),1./sqrt(2)]
	Edge3                 // Edge from vertex 2-0, Normal is [-1,0]
)

func (rt *RTElement) getLocationType(i int) (rtt RTPointType) {
	var (
		NpInt  = rt.NpInt
		NpEdge = rt.NpEdge
	)
	switch {
	case i < NpInt:
		// Unit vector is [1,0]
		rtt = InteriorR
	case i >= NpInt && i < 2*NpInt:
		// Unit vector is [0,1]
		rtt = InteriorS
	case i >= 2*NpInt && i < 2*NpInt+NpEdge:
		// Edge1: Unit vector is [0,-1]
		rtt = Edge1
	case i >= 2*NpInt+NpEdge && i < 2*NpInt+2*NpEdge:
		// Edge2: Unit vector is [1/sqrt(2), 1/sqrt(2)]
		rtt = Edge2
	case i >= 2*NpInt+2*NpEdge && i < 2*NpInt+3*NpEdge:
		// Edge3: Unit vector is [-1,0]
		rtt = Edge3
	}
	return
}
