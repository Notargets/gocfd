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

		Each point in this basis corresponds to a basis vector with 2 components
							e.g. [v]j = [v1,v2]j

		Interior points in [R,S] each have two basis vectors and each basis
		vector has two components (it's a vector!). That means the interior
		point indices are duplicated for the interior points, because we need
		to refer to each basis vector independently.

		Edge points in [R,S] each have one basis vector, and each basis vector
		has two components (it's a vector!).

		When we evaluate a basis function at a point in [R,S], we will generally
		take the dot product of that basis function [v]j with the basis vector
		at that location [v]i. When we project the physical flux [F] at a point
		i within the element, we'll take the dot product between [F]j and the
		basis vector [v]j to do the projection: [v]j ⋅ [F]j

		The indexing of the RT element is: first come the "left" basis vectors
		for each of the interior points, then come the "right" basis vectors for
		each of the interior points, then come each of the edge basis vectors.

		First (N)(N+1)/2 points:
            Interior points, excluding edges in [-1,1] triangle coordinates
            corresponding to the "left" component DOFs
		Next (N)(N+1)/2 points: Same as above for the "right" component DOFs

		Edges run in counterclockwise order for right hand rule reasons,
		the indices start at bottom left, run right along edge 1, then up+left
		along hypotenuse edge 2, then down along edge 3

		Next (N+1) points: Edge 1 (Bottom) locations
		Next (N+1) points: Edge 2 (Hypotenuse) locations
		Next (N+1) points: Edge 3 (Left) locations

	Methods:
		The evaluation matrix rt.A relates the coefficients of the RT element
		basis to the solution vectors. For the RT element at degree 1 (RT1),
		this equation is:

⎡ ϕ₁ ⋅ ϕ₁   ϕ₂ ⋅ ϕ₁   ϕ₃ ⋅ ϕ₁   ϕ₄ ⋅ ϕ₁   ϕ₅ ⋅ ϕ₁   ϕ₆ ⋅ ϕ₁   ϕ₇ ⋅ ϕ₁   ϕ₈ ⋅ ϕ₁ ⎤   ⎡ c₁ ⎤   ⎡ ϕ₁ ⋅ f₁ ⎤
⎢ ϕ₁ ⋅ ϕ₂   ϕ₂ ⋅ ϕ₂   ϕ₃ ⋅ ϕ₂   ϕ₄ ⋅ ϕ₂   ϕ₅ ⋅ ϕ₂   ϕ₆ ⋅ ϕ₂   ϕ₇ ⋅ ϕ₂   ϕ₈ ⋅ ϕ₂ ⎥   ⎢ c₂ ⎥   ⎢ ϕ₂ ⋅ f₂ ⎥
⎢ ϕ₁ ⋅ ϕ₃   ϕ₂ ⋅ ϕ₃   ϕ₃ ⋅ ϕ₃   ϕ₄ ⋅ ϕ₃   ϕ₅ ⋅ ϕ₃   ϕ₆ ⋅ ϕ₃   ϕ₇ ⋅ ϕ₃   ϕ₈ ⋅ ϕ₃ ⎥   ⎢ c₃ ⎥   ⎢ ϕ₃ ⋅ f₃ ⎥
⎢ ϕ₁ ⋅ ϕ₄   ϕ₂ ⋅ ϕ₄   ϕ₃ ⋅ ϕ₄   ϕ₄ ⋅ ϕ₄   ϕ₅ ⋅ ϕ₄   ϕ₆ ⋅ ϕ₄   ϕ₇ ⋅ ϕ₄   ϕ₈ ⋅ ϕ₄ ⎥   ⎢ c₄ ⎥   ⎢ ϕ₄ ⋅ f₄ ⎥
⎢ ϕ₁ ⋅ ϕ₅   ϕ₂ ⋅ ϕ₅   ϕ₃ ⋅ ϕ₅   ϕ₄ ⋅ ϕ₅   ϕ₅ ⋅ ϕ₅   ϕ₆ ⋅ ϕ₅   ϕ₇ ⋅ ϕ₅   ϕ₈ ⋅ ϕ₅ ⎥ * ⎢ c₅ ⎥ = ⎢ ϕ₅ ⋅ f₅ ⎥
⎢ ϕ₁ ⋅ ϕ₆   ϕ₂ ⋅ ϕ₆   ϕ₃ ⋅ ϕ₆   ϕ₄ ⋅ ϕ₆   ϕ₅ ⋅ ϕ₆   ϕ₆ ⋅ ϕ₆   ϕ₇ ⋅ ϕ₆   ϕ₈ ⋅ ϕ₆ ⎥   ⎢ c₆ ⎥   ⎢ ϕ₆ ⋅ f₆ ⎥
⎢ ϕ₁ ⋅ ϕ₇   ϕ₂ ⋅ ϕ₇   ϕ₃ ⋅ ϕ₇   ϕ₄ ⋅ ϕ₇   ϕ₅ ⋅ ϕ₇   ϕ₆ ⋅ ϕ₇   ϕ₇ ⋅ ϕ₇   ϕ₈ ⋅ ϕ₇ ⎥   ⎢ c₇ ⎥   ⎢ ϕ₇ ⋅ f₇ ⎥
⎣ ϕ₁ ⋅ ϕ₈   ϕ₂ ⋅ ϕ₈   ϕ₃ ⋅ ϕ₈   ϕ₄ ⋅ ϕ₈   ϕ₅ ⋅ ϕ₈   ϕ₆ ⋅ ϕ₈   ϕ₇ ⋅ ϕ₈   ϕ₈ ⋅ ϕ₈ ⎦   ⎢ c₈ ⎥   ⎢ ϕ₈ ⋅ f₈ ⎥


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
		edge. The definition of the interior basis functions establish that the
		dot product of the interior basis functions with the edge functions
		automatically is zero. Additionally, each edge function is implemented
		using Lagrange basis functions that ensure each edge basis dot product
		with other edge functions is zero. So the resulting evaluation matrix
		will feature zeros that enforce orthogonality of the proper parts of
		the basis.
*/
/*
⎡ ϕ₁ ⋅ ϕ₁   ϕ₂ ⋅ ϕ₁      0         0         0         0         0         0    ⎤   ⎡ c₁ ⎤   ⎡ ϕ₁ ⋅ f₁ ⎤
⎢ ϕ₁ ⋅ ϕ₂   ϕ₂ ⋅ ϕ₂      0         0         0         0         0         0    ⎥   ⎢ c₂ ⎥   ⎢ ϕ₂ ⋅ f₂ ⎥
⎢    0         0      ϕ₃ ⋅ ϕ₃      0         0         0         0         0    ⎥   ⎢ c₃ ⎥   ⎢ ϕ₃ ⋅ f₃ ⎥
⎢    0         0         0      ϕ₄ ⋅ ϕ₄      0         0         0         0    ⎥   ⎢ c₄ ⎥   ⎢ ϕ₄ ⋅ f₄ ⎥
⎢    0         0         0         0      ϕ₅ ⋅ ϕ₅      0         0         0    ⎥ * ⎢ c₅ ⎥ = ⎢ ϕ₅ ⋅ f₅ ⎥
⎢    0         0         0         0         0      ϕ₆ ⋅ ϕ₆      0         0    ⎥   ⎢ c₆ ⎥   ⎢ ϕ₆ ⋅ f₆ ⎥
⎢    0         0         0         0         0         0      ϕ₇ ⋅ ϕ₇      0    ⎥   ⎢ c₇ ⎥   ⎢ ϕ₇ ⋅ f₇ ⎥
⎣    0         0         0         0         0         0         0      ϕ₈ ⋅ ϕ₈ ⎦   ⎢ c₈ ⎥   ⎢ ϕ₈ ⋅ f₈ ⎥
*/
type DerivativeDirection uint8

const (
	None DerivativeDirection = iota
	Dr
	Ds
)

type RTElement struct {
	P             int // Order of element
	Np            int // Number of points in element
	NpEdge, NpInt int // Number of Edge and Interior points
	// At every point, the sum of basis functions equals the flux vector
	// Each basis function is multiplied by it's constant, Ci
	// [A] is the matrix of basis functions, evaluated at each point (row)
	// [A] relates the constants [C] to the flux vector [F]
	// [A] x [C] = [F]
	// ===> [C] = [AInv] x [F]
	A, AInv utils.Matrix // Basis evaluation matrix, NpxNp
	// The divergence of [F] at every point is the sum of the basis derivatives
	// Div[F] = Dr([F])+Ds([F]) = Dr([A]x[C])+Ds([A]x[C]) = [Dr]x[C]+[Ds]x[C]
	// Div[F] = ([Dr]+[Ds]) x [C] = ([Dr]+[Ds]) x [AInv] x [F]
	Dr, Ds utils.Matrix // Derivative of basis functions
	// Commutation gives a useful matrix we can use to calculate Flux Divergence
	// [Div] == ([Dr] + [Ds]) x [AInv]
	// Div[F] = [Div] x [F]
	Div        utils.Matrix // Divergence matrix, NpxNp for all, NintxNp Interior Points
	R, S       utils.Vector // Point locations defining element in [-1,1] Triangle, NpxNp
	RInt, SInt utils.Vector
}

func NewRTElement(P int) (rt *RTElement) {
	// We expect that there are points in R and S to match the dimension of dim(P(NFlux-1))
	/*
		<---- NpInt ----><---- NpInt ----><---NpEdge----><---NpEdge----><---NpEdge---->
		         Solution Points          Edge 1 pts	Edge 2 pts	  Edge 3 pts
		<---- NpInt ----><---- NpInt ----><---NpEdge----><---NpEdge----><---NpEdge---->
	*/
	var (
		Np = (P + 1) * (P + 3)
	)
	rt = &RTElement{
		P:      P,
		Np:     Np,
		NpInt:  P * (P + 1) / 2, // Number of interior points is same as the 2D scalar space one order lesser
		NpEdge: P + 1,           // Each edge is P+1 nodes
		A:      utils.NewMatrix(Np, Np),
		AInv:   utils.NewMatrix(Np, Np),
		Dr:     utils.NewMatrix(Np, Np),
		Ds:     utils.NewMatrix(Np, Np),
		Div:    utils.NewMatrix(Np, Np),
	}
	if P > 0 {
		if P < 9 {
			rt.RInt, rt.SInt = NodesEpsilon(P - 1)
		} else {
			rt.RInt, rt.SInt = XYtoRS(Nodes2D(P - 1))
		}
	}
	rt.R, rt.S = rt.ExtendGeomToRT(rt.RInt, rt.SInt)
	// For edges 1 and 3, we can use the same definition on R, as S is the same
	rt.CalculateBasis()
	return
}

/*
This follows the development in the V.J. Ervin paper: Computational Bases for
RTk and BDMk on Triangles
All are defined in terms of Eta and Xi which go from 0 to 1 compared to here
where R and S are defined from -1 to 1

	xi = (s + 1)/2.
	eta = (r + 1)/2.

There are three "Normal" or "edge" functions and on each of the 3 triangle edges
there are some number of points with a single vector basis function.

For the RT2 element, there are 3 points on each of the 3 edges. The vectors are
defined below such that each column is one edge, having three rows, each
defining the basis function for that point on the edge:
Φ1 (ξ, η) = q1 (η) ê1 (ξ, η) , Φ2 (ξ, η) = q2 (η) ê1 (ξ, η) , Φ3 (ξ, η) = q3 (η) ê1 (ξ, η) ,
Φ1 (ξ, η) = q3 (η) ê2 (ξ, η) , Φ2 (ξ, η) = q2 (η) ê2 (ξ, η) , Φ3 (ξ, η) = q1 (η) ê2 (ξ, η) ,
Φ1 (ξ, η) = q1 (ξ) ê3 (ξ, η) , Φ2 (ξ, η) = q2 (ξ) ê3 (ξ, η) , Φ3 (ξ, η) = q3 (ξ) ê3 (ξ, η) .

For the RT2 element, there are 3 interior points, and in the below example, we
see there are two basis functions (vectors) for each of the 3 points. Each
point is a column, and each column has two rows defining the basis vectors for
that point:
Φ1 (ξ, η) = (1 − ξ − η) ê4 (ξ, η) , Φ2 (ξ, η) = ξ ê4 (ξ, η) , Φ3 (ξ, η) = η ê4 (ξ, η) ,
Φ1 (ξ, η) = (1 − ξ − η) ê5 (ξ, η) , Φ2 (ξ, η) = ξ ê5 (ξ, η) , Φ3 (ξ, η) = η ê5 (ξ, η) .

Given that the edge functions are only evaluated on the edges, and that the
lagrange functions multiplying them cause them to vanish on all but their
defining point, it seems unnecessary to evaluate anything but the e1-e3
functions in practice.

For the general case of the interior basis functions, they take the form:
Φj (ξ, η) = bj (ξ, η) ê4 (ξ, η) , j = 1, 2, . . . , k(k + 1)/2 ,
Φj (ξ, η) = bj (ξ, η) ê5 (ξ, η) , j = 1, 2, . . . , k(k + 1)/2 .

where bj is a 2D polynomial of order P-1
*/
func baseBasisFunctions(r, s float64, edgeNum int) (ef [2]float64) {
	// Mapping:
	// Edge 1 (Left edge)   => Ervin Edge 2
	// Edge 2 (Hypotenuse)  => Ervin Edge 1
	// Edge 3 (Bottom edge) => Ervin Edge 3
	// Note: Here we fix an error in Ervin - the midpoint of RT0 edges are now
	// correctly unit normal
	var (
		sr2 = math.Sqrt(2.)
		// Transform our unit triangle coords to Ervin:
		xi  = 0.5 * (r + 1)
		eta = 0.5 * (s + 1)
	)
	switch edgeNum {
	case 2:
		ef = [2]float64{sr2 * xi, sr2 * eta}
	case 1:
		ef = [2]float64{xi - 1, eta - 0.5}
	case 3:
		ef = [2]float64{xi - 0.5, eta - 1}
	case 4:
		ef = [2]float64{eta * xi, eta * (eta - 1)}
	case 5:
		ef = [2]float64{xi * (xi - 1), xi * eta}
	default:
		panic("wrong basis function number (1-5)")
	}
	return
}
func divergenceOfEdgeFunctions(r, s float64, edgeNum, j int) (div float64) {
	// ---------------------------------------------------
	// the edge basis vector function [v] varies along the edge.
	// It is the product of a 1D edge function f(xi) and [v], so we have:
	// div(edgeFunction) = div(f(xi)*[v]) =
	//       [df(xi)/dr,df(xi)/ds] ⋅ [v] + f(xi) * ([div] ⋅ [v])
	//
	// div = df(xi)/dxi * (v1*dxi/dr + v2*dxi/ds) + f(xi) * (dv1/dr + dv2/ds)
	//
	// Conversion from Ervin coordinates:
	// xi  = 0.5 * (r + 1)
	// eta = 0.5 * (s + 1)
	// r = 2 * xi  - 1
	// s = 2 * eta - 1
	//
	// Left Edge (1) divergence:
	// 			[v] = [(r - 1)/2, s/2]
	// v1 = (r-1)/2, v2 = s/2, dv1/dr = 1/2, dv2/ds = 1/2
	// The left edge  in a counter-clockwise direction is parameterized:
	// r = -1 (constant), xi = -s => dxi/dr = 0, dxi/ds = -1,
	// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
	//     = df/dxi*(v1*(  0   )+v2*(  -1  )) + f(xi) * (  1/2  +   1/2 )
	//     = df/dxi*(          v2           ) + f(xi)
	//         div(edge1) = (df/dxi) * v2 + f(xi)
	//
	// Hypotenuse (2) divergence:
	// 			[v] = [Sqrt2/2 * (r+1), Sqrt2/2 * (s+1)]
	// v1 = Sqrt2/2 * (r+1), v2 = Sqrt2/2 * (s+1), dv1/dr = Sqrt2/2 = dv2/ds
	//
	// The hypotenuse in a counter-clockwise direction is parameterized:
	// xi = -r = s, => dxi/dr = -1, dxi/ds = 1
	//
	// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
	//     = df/dxi*(v1*( -1   )+v2*(   1  )) + f(xi) * (Sqrt2/2+Sqrt2/2)
	//     = df/dxi*(         v2-v1         ) + f(xi) * Sqrt2
	//         div(edge2) = (df/dxi) * (v2-v1) + Sqrt2 * f(xi)
	//
	// Bottom Edge (3) divergence:
	// 			[v] = [r/2, (s - 1)/2]
	// v1 = r/2, v2 = (s-1)/2, dv1/dr = 1/2, dv2/ds = 1/2
	// The bottom edge  in a counter-clockwise direction is parameterized:
	// xi = r, s = -1 (const) => dxi/dr = 1, dxi/ds = 0
	// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
	//     = df/dxi*(v1*(  1   )+v2*(   0  )) + f(xi) * (  1/2  +   1/2 )
	//     = df/dxi*(          v1           ) + f(xi)
	//        div(edge3) = (df/dxi) * v1 + f(xi)
	return
}

func (rt *RTElement) basisPolynomialValue(r, s float64, j int,
	derivO ...DerivativeDirection) (val float64) {
	// This evaluates the j-th polynomial or derivative at r,
	// s used to multiply each of the 5 base basis vectors in Ervin's RT basis.
	// The input parameter "j" is the function index within the element
	var (
		NpInt   = rt.NpInt
		NpEdge  = rt.NpEdge
		funcNum int
	)
	switch {
	case j >= 0 && j < NpInt:
		funcNum = 4 // Ervin's polynomial 4 which multiplies the 4th vector
	case j >= NpInt && j < 2*NpInt:
		funcNum = 5 // Ervin's polynomial 5 which multiplies the 5th vector
	case j >= 2*NpInt && j < 2*NpInt+NpEdge:
		funcNum = 1
	case j >= 2*NpInt+NpEdge && j < 2*NpInt+2*NpEdge:
		funcNum = 2
	case j >= 2*NpInt+2*NpEdge && j < 2*NpInt+3*NpEdge:
		funcNum = 3
	default:
		panic("j polynomial index out of range")
	}
	var Xi []float64
	var deriv = 0
	var jb2d *JacobiBasis2D
	if funcNum < 4 {
		// Parameterized edge coordinate Xi
		Xi = rt.getEdgeXiParameter(funcNum).DataP
		if len(derivO) > 0 {
			deriv = 1
		}
	} else if funcNum <= 5 {
		// Calculate alpha based on the value of j within NpInt, scaled to run
		// Beta from -0.5 to 0.5.
		// j runs from 0 to NpInt-1 for this function
		var jj int
		if funcNum == 5 {
			jj = j - NpInt
		} else {
			jj = j
		}
		alpha := 0.
		beta := float64(jj)/(float64(NpInt)-1.) - 0.5
		jb2d = NewJacobiBasis2D(rt.P-1, rt.RInt, rt.SInt, alpha, beta)
	}
	switch funcNum {
	case 1: // Edge 1
		if r != -1 {
			panic("edge 1 polynomial r must be -1")
		}
		// jj is the edge local index
		jj := j - 2*rt.NpInt
		xi := -s
		val = DG1D.Lagrange1DPoly(xi, Xi, jj, deriv)
	case 2: // Edge 2
		if r+s != 0 {
			panic("edge 2 polynomial must have s+r = 0")
		}
		// jj is the edge local index
		jj := j - 2*rt.NpInt - rt.NpEdge
		xi := -r
		val = DG1D.Lagrange1DPoly(xi, Xi, jj, deriv)
	case 3: // Edge 3
		if s != -1 {
			panic("edge 3 polynomial must have s = -1")
		}
		// jj is the edge local index
		jj := j - 2*rt.NpInt - 2*rt.NpEdge
		xi := r
		val = DG1D.Lagrange1DPoly(xi, Xi, jj, deriv)
	case 4:
		fallthrough
	case 5:
		// We represent the interior polynomial for the special cases of the
		// RT1 and RT2 elements. For RT3 and above, we use a Jacobi2D
		// polynomial, but for RT1 it's a constant and for RT2 it's a basic poly set
		var derivA DerivativeDirection
		if len(derivO) > 0 {
			derivA = derivO[0]
		}
		switch rt.P {
		case 0:
			panic("RT element isn't defined here for RT0")
		case 1:
			if derivA == None {
				val = 1
			} else {
				val = 0
			}
		case 2:
			var jj int
			if funcNum == 5 {
				jj = j - NpInt
			} else {
				jj = j
			}
			switch derivA {
			case None:
				xi := 0.5 * (r + 1)
				eta := 0.5 * (s + 1)
				switch jj {
				case 0:
					val = 1. - xi - eta
				case 1:
					val = xi
				case 2:
					val = eta
				}
			case Dr:
				switch jj {
				case 0:
					val = -0.5
				case 1:
					val = 0.5
				case 2:
					val = 0
				}
			case Ds:
				switch jj {
				case 0:
					val = -0.5
				case 1:
					val = 0
				case 2:
					val = 0.5
				}
			}
		default:
			val = jb2d.GetPolynomialEvaluation(r, s, derivO...)
		}
	default:
		panic("wrong edge number")
	}
	return
}

func (rt *RTElement) getEdgeXiParameter(edgeNum int) (Xi utils.Vector) {
	// Edge 1 is S=-1 (bottom of tri)
	// Edge 2 is Hypotenuse
	// Edge 3 is R=-1 (left side of tri)
	switch edgeNum {
	case 1:
		Xi = rt.R.Subset(2*rt.NpInt, 2*rt.NpInt+rt.NpEdge-1)
	case 2:
		Xi = rt.R.Subset(2*rt.NpInt+rt.NpEdge, 2*rt.NpInt+2*rt.NpEdge-1)
	case 3:
		Xi = rt.S.Subset(2*rt.NpInt+2*rt.NpEdge, 2*rt.NpInt+3*rt.NpEdge-1)
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
		// psi = rt.RTPolyBasis2D_A.GetPolynomialEvaluation(r, s, derivO...)
	case InteriorS:
		// psi = rt.RTPolyBasis2D_B.GetPolynomialEvaluation(r, s, derivO...)
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
