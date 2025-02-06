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

type RTFunctionNumber uint8

const (
	All RTFunctionNumber = iota
	E4                   // Ervin E4 Interior basis
	E5                   // Ervin E5 Interior basis
	E1                   // Bottom Edge, Normal is [0,-1], Ervin E3
	E2                   // Hypotenuse, Normal is [1/Sqrt2,1/Sqrt2], Ervin E1
	E3                   // Left Edge, Normal is [-1,0], Ervin E2
)

func (fn RTFunctionNumber) String() string {
	switch fn {
	case All:
		return "All"
	case E4:
		return "E4"
	case E5:
		return "E5"
	case E1:
		return "E1"
	case E2:
		return "E2"
	case E3:
		return "E3"
	default:
		return "Unknown"
	}
}

type RTElement struct {
	P             int // Order of element
	Np            int // Number of points in element
	NpEdge, NpInt int // Number of Edge and Interior points
	// At every point, the sum of basis functions equals the flux vector
	// Each basis function is multiplied by its constant, Ci
	// [A] is the matrix of basis functions, evaluated at each point (row)
	// [A] relates the constants [C] to the flux vector [F]
	// [A] x [C] = [F]
	// ===> [C] = [AInv] x [F]
	V, VInv utils.Matrix // Vandermonde Matrix
	A, AInv utils.Matrix // Basis evaluation matrix, NpxNp
	// The divergence of [F] at every point is the sum of the basis derivatives
	// Div[F] = Dr([F])+Ds([F]) = Dr([A]x[C])+Ds([A]x[C]) = [Dr]x[C]+[Ds]x[C]
	// Div[F] = ([Dr]+[Ds]) x [C] = ([Dr]+[Ds]) x [AInv] x [F]
	Dr, Ds utils.Matrix // Derivative of basis functions
	// Commutation gives a useful matrix we can use to calculate Flux Divergence
	// [Div] == ([Dr] + [Ds]) x [AInv]
	// Div[F] = [Div] x [F]
	Div         utils.Matrix // Divergence matrix, NpxNp for all, NintxNp Interior Points
	R, S        utils.Vector // Point locations defining element in [-1,1] Triangle, NpxNp
	RInt, SInt  utils.Vector
	BasisVector [2]utils.Matrix
	Projection  utils.Matrix
	// Technically, each edge is
	lp1d         *LagrangePolynomial1D
	lp2d         *LagrangePolynomial2D
	jb2d         *JacobiBasis2D
	PolyTerms    utils.Matrix // Each column is a polynomial term evaluated at rows
	BasisVectors []BasisVectorStruct
	Phi          []BasisPolynomialTerm
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
		V:      utils.NewMatrix(Np, Np),
		A:      utils.NewMatrix(Np, Np),
		AInv:   utils.NewMatrix(Np, Np),
		Dr:     utils.NewMatrix(Np, Np),
		Ds:     utils.NewMatrix(Np, Np),
		Div:    utils.NewMatrix(Np, Np),
		BasisVector: [2]utils.Matrix{utils.NewMatrix(Np, 1),
			utils.NewMatrix(Np, 1)},
		Projection: utils.NewMatrix(Np, 1),
		PolyTerms:  utils.NewMatrix(Np, Np),
	}

	if P > 0 {
		if P < 9 {
			rt.RInt, rt.SInt = NodesEpsilon(P - 1)
		} else {
			rt.RInt, rt.SInt = XYtoRS(Nodes2D(P - 1))
		}
	}
	if P > 2 {
		rt.lp2d = NewLagrangePolynomialBasis2D(rt.P-1, rt.RInt, rt.SInt)
	}

	rt.R, rt.S = rt.ExtendGeomToRT(rt.RInt, rt.SInt)
	// Each edge is identical WRT point distribution along the edge,
	// so we just need one lagrange polynomial for each of them
	edgePoints := rt.R.Subset(2*rt.NpInt, 2*rt.NpInt+rt.NpEdge-1)
	rt.lp1d = NewLagrangePolynomial1D(edgePoints, rt.P, 0, 0)

	// fmt.Printf("RT%d - Calculating Basis...", P)
	rt.CalculateBasis()
	rt.V = rt.ComposeV()
	rt.VInv = rt.V.InverseWithCheck()
	// fmt.Printf("\nRT%d - Calculating Divergence Matrix...", P)
	rt.Div = rt.ComputeDivergenceMatrix()
	// fmt.Printf("Done.\n")
	return
}

type BasisVectorStruct struct {
	Eval       func(r, s float64) (v [2]float64)
	Dot        func(r, s float64, f [2]float64) (dot float64)
	Project    func(r, s float64, psi float64) (v [2]float64) // scalar mult psi
	Divergence func(r, s float64) (div float64)               // Div of vector
	Sum        func(r, s float64) (sum float64)               // Sum of vector
}

type BasisPolynomialMultiplier struct {
	// This multiplies the BasisVector to produce a term
	Eval        func(r, s float64) (val float64)
	Divergence  func(r, s float64) (div float64)
	OrderOfTerm int
}
type BasisPolynomialTerm struct {
	IsScaled       bool // After scaling, is part of the element polynomial
	PolyMultiplier BasisPolynomialMultiplier
	BasisVector    BasisVectorStruct
	Coefficient    float64 // Initially set to 1. to allow scaling calculations
}

func (pt BasisPolynomialTerm) Eval(r, s float64) (v [2]float64) {
	scale := pt.PolyMultiplier.Eval(r, s) * pt.Coefficient
	v = pt.BasisVector.Project(r, s, scale)
	return
}
func (pt BasisPolynomialTerm) Dot(r, s float64, b [2]float64) (dot float64) {
	v := pt.Eval(r, s)
	dot = b[0]*v[0] + b[1]*v[1]
	return
}

func (pt BasisPolynomialTerm) Divergence(r, s float64) (div float64) {
	var (
		divBasis = pt.BasisVector.Divergence(r, s)
		divPoly  = pt.PolyMultiplier.Divergence(r, s)
		polyEval = pt.PolyMultiplier.Eval(r, s)
		sumBasis = pt.BasisVector.Sum(r, s)
	)
	div = polyEval*divBasis + divPoly*sumBasis
	return
}

func (rt *RTElement) CalculateBasis() {
	var (
		v [2]float64
	)
	switch rt.P {
	case 1:
		rt.Phi = rt.ComposePhiRT1()
	}
	// fmt.Printf("RT%d Element Basis\n", P)
	// fmt.Printf("Np:%d; NpInt:%d; NpEdge:%d\n", rt.Np, rt.NpInt, rt.NpEdge)
	for j := 0; j < rt.Np; j++ {
		r, s := rt.R.AtVec(j), rt.S.AtVec(j)
		v, _ = rt.basisEvaluation(r, s, j)
		rt.BasisVector[0].Set(j, 0, v[0])
		rt.BasisVector[1].Set(j, 0, v[1])
	}
	return
}

var (
	sr2   = math.Sqrt2
	conv  = func(r float64) float64 { return (r + 1.) / 2. }
	Dot   = func(v1, v2 [2]float64) float64 { return v1[0]*v2[0] + v1[1]*v2[1] }
	Scale = func(v [2]float64, scale float64) [2]float64 {
		return [2]float64{v[0] * scale, v[1] * scale}
	}
	Sum       = func(v [2]float64) float64 { return v[0] + v[1] }
	normalize = func(v [2]float64) (normed [2]float64) {
		norm := 1. / math.Sqrt(v[0]*v[0]+v[1]*v[1])
		normed = [2]float64{v[0] * norm, v[1] * norm}
		return
	}
	e1v = func(r, s float64) (v [2]float64) {
		xi, eta := conv(r), conv(s)
		v = [2]float64{xi, eta - 1.}
		return
	}
	e2v = func(r, s float64) (v [2]float64) {
		sr2 := math.Sqrt2
		xi, eta := conv(r), conv(s)
		v = [2]float64{sr2 * xi, sr2 * eta}
		return
	}
	e3v = func(r, s float64) (v [2]float64) {
		xi, eta := conv(r), conv(s)
		v = [2]float64{xi - 1., eta}
		return
	}
	e4v = func(r, s float64) (v [2]float64) {
		xi, eta := conv(r), conv(s)
		v = [2]float64{eta * xi, eta * (eta - 1.)}
		// v = normalize(v)
		return
	}
	e5v = func(r, s float64) (v [2]float64) {
		xi, eta := conv(r), conv(s)
		v = [2]float64{xi * (xi - 1.), eta * xi}
		// v = normalize(v)
		return
	}
	E4Vector = BasisVectorStruct{
		// Interior vector E4
		Eval: func(r, s float64) [2]float64 { return e4v(r, s) },
		Dot: func(r, s float64, f [2]float64) float64 {
			return Dot(e4v(r, s), f)
		},
		Project: func(r, s float64, scale float64) [2]float64 {
			return Scale(e4v(r, s), scale)
		},
		Divergence: func(r, s float64) (div float64) {
			div = (3.*s + 1.) / 4.
			return
		},
		Sum: func(r, s float64) float64 { return Sum(e4v(r, s)) },
	}
	E5Vector = BasisVectorStruct{
		// Interior vector E5
		Eval: func(r, s float64) [2]float64 { return e5v(r, s) },
		Dot: func(r, s float64, f [2]float64) (dot float64) {
			return Dot(e5v(r, s), f)
		},
		Project: func(r, s float64, scale float64) (v [2]float64) {
			return Scale(e5v(r, s), scale)
		},
		Divergence: func(r, s float64) (div float64) {
			div = (3.*r + 1.) / 4.
			return
		},
		Sum: func(r, s float64) float64 { return Sum(e5v(r, s)) },
	}
	Edge1NormalVector = BasisVectorStruct{
		// Bottom edge
		Eval: func(r, s float64) (v [2]float64) { return e1v(r, s) },
		Dot: func(r, s float64, f [2]float64) (dot float64) {
			return Dot(e1v(r, s), f)
		},
		Project: func(r, s float64, scale float64) (v [2]float64) {
			return Scale(e1v(r, s), scale)
		},
		Divergence: func(r, s float64) (div float64) { return },
		Sum:        func(r, s float64) float64 { return Sum(e1v(r, s)) },
	}
	Edge2NormalVector = BasisVectorStruct{
		// Hypotenuse
		Eval: func(r, s float64) (v [2]float64) { return e2v(r, s) },
		Dot: func(r, s float64, f [2]float64) (dot float64) {
			return Dot(e2v(r, s), f)
		},
		Project: func(r, s float64, scale float64) (v [2]float64) {
			return Scale(e2v(r, s), scale)
		},
		Divergence: func(r, s float64) (div float64) { return },
		Sum:        func(r, s float64) float64 { return Sum(e2v(r, s)) },
	}
	Edge3NormalVector = BasisVectorStruct{
		// Left edge
		Eval: func(r, s float64) (v [2]float64) { return e3v(r, s) },
		Dot: func(r, s float64, f [2]float64) (dot float64) {
			return Dot(e3v(r, s), f)
		},
		Project: func(r, s float64, scale float64) (v [2]float64) {
			return Scale(e3v(r, s), scale)
		},
		Divergence: func(r, s float64) (div float64) { return },
		Sum:        func(r, s float64) float64 { return Sum(e3v(r, s)) },
	}
)

func (rt *RTElement) ComposeV() (V utils.Matrix) {
	var (
		Np   = rt.Np
		R, S = rt.R, rt.S
	)
	V = utils.NewMatrix(Np, Np)
	for i := 0; i < Np; i++ {
		r_i, s_i := R.DataP[i], S.DataP[i]
		b_i := rt.Phi[i].BasisVector
		for j := 0; j < Np; j++ {
			// Don't evaluate edge basis on other edges
			switch rt.getFunctionNumber(j) {
			case E1, E2, E3:
				switch rt.getFunctionNumber(i) {
				case E1, E2, E3:
					if rt.getFunctionNumber(i) != rt.getFunctionNumber(j) {
						continue
					}
				}
			}
			// Don't evaluate Interior basis on edges (likely wrong)
			// switch rt.getFunctionNumber(i) {
			// case E4, E5:
			// 	switch rt.getFunctionNumber(j) {
			// 	case E1, E2, E3:
			// 		continue
			// 	}
			// }
			v_j := rt.Phi[j].Eval(r_i, s_i)
			V.Set(i, j, b_i.Dot(r_i, s_i, v_j))
		}
	}
	return
}

func (rt *RTElement) ComposePhiRT1() (phi []BasisPolynomialTerm) {
	// Compose the basis for the RT1 element based on Ervin with normal edges
	// Get the two edge points for use in the lagrange terms
	var (
		g1   = rt.R.AtVec(2 * rt.NpInt)
		g2   = rt.R.AtVec(2*rt.NpInt + 1)
		l1xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return (r - g2) / (g1 - g2)
			},
			Divergence: func(r, s float64) (div float64) {
				div = 1 / (g1 - g2)
				return
			},
			OrderOfTerm: 1,
		}
		l1eta = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return (s - g2) / (g1 - g2)
			},
			Divergence: func(r, s float64) (div float64) {
				div = 1 / (g1 - g2)
				return
			},
			OrderOfTerm: 1,
		}
		l2xi = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return (r - g1) / (g2 - g1)
			},
			Divergence: func(r, s float64) (div float64) {
				div = 1 / (g2 - g1)
				return
			},
			OrderOfTerm: 1,
		}
		l2eta = BasisPolynomialMultiplier{
			Eval: func(r, s float64) float64 {
				return (s - g1) / (g2 - g1)
			},
			Divergence: func(r, s float64) (div float64) {
				div = 1 / (g2 - g1)
				return
			},
			OrderOfTerm: 1,
		}
		constant = BasisPolynomialMultiplier{
			Eval:        func(r, s float64) float64 { return 1. },
			Divergence:  func(r, s float64) (div float64) { return 0. },
			OrderOfTerm: 0,
		}
	)
	phi = []BasisPolynomialTerm{
		{
			Coefficient:    1,
			PolyMultiplier: constant,
			BasisVector:    E4Vector,
		},
		{
			Coefficient:    1,
			PolyMultiplier: constant,
			BasisVector:    E5Vector,
		},
		{
			Coefficient:    1,
			PolyMultiplier: l1xi,
			BasisVector:    Edge1NormalVector,
		},
		{
			Coefficient:    1,
			PolyMultiplier: l2xi,
			BasisVector:    Edge1NormalVector,
		},
		{
			Coefficient:    1,
			PolyMultiplier: l1eta,
			BasisVector:    Edge2NormalVector,
		},
		{
			Coefficient:    1,
			PolyMultiplier: l2eta,
			BasisVector:    Edge2NormalVector,
		},
		{
			Coefficient:    1,
			PolyMultiplier: l2eta,
			BasisVector:    Edge3NormalVector,
		},
		{
			Coefficient:    1,
			PolyMultiplier: l1eta,
			BasisVector:    Edge3NormalVector,
		},
	}
	return
}

func (rt *RTElement) ComputeDivergenceMatrix() (Div utils.Matrix) {
	// This must be called before the element is "finished" - we break this out
	// separately to allow for testing
	// BasisDotInverse := rt.ComputeBasisDotInverse()
	BasisDotInverse := rt.V.InverseWithCheck()

	// Build basis divergence matrix
	DivBasis := utils.NewMatrix(rt.Np, rt.Np)
	for i := 0; i < rt.Np; i++ {
		r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		for j := 0; j < rt.Np; j++ {
			// fmt.Printf("NpInt, j = %d, %d\n", rt.NpInt, j)
			DivBasis.Set(i, j, rt.Phi[j].Divergence(r, s))
		}
	}
	// DivBasis.Print("Basis Divergence Matrix")
	Div = DivBasis.Mul(BasisDotInverse)
	return
}

func (rt *RTElement) ComputeBasisDotInverse() (BasisDotInverse utils.Matrix) {
	// This is the Basis Vector dotted with other basis functions per location
	// The inverse enables the computation of the constants that project a
	// target vector field onto the element.
	// If this matrix does not invert, it reflects issues with the basis
	BasisMatrix := [2]utils.Matrix{utils.NewMatrix(rt.Np, rt.Np),
		utils.NewMatrix(rt.Np, rt.Np)}
	for i := 0; i < rt.Np; i++ {
		r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		v := [2]float64{rt.BasisVector[0].At(i, 0), rt.BasisVector[1].At(i, 0)}
		// fmt.Printf("Base[%d] = [%f,%f]\n", i, v[0], v[1])
		for j := 0; j < rt.Np; j++ {
			v2, _ := rt.basisEvaluation(r, s, j)
			// fmt.Printf("Psi[%d,%d] = [%f,%f]\n", i, j, v2[0], v2[1])
			BasisMatrix[0].Set(i, j, v[0]*v2[0])
			BasisMatrix[1].Set(i, j, v[1]*v2[1])
		}
	}
	BasisDot := BasisMatrix[0].Add(BasisMatrix[1])
	// BasisDot.Print("BasisDot")
	BasisDotInverse = BasisDot.InverseWithCheck()
	// BasisDotInverse.Print("BasisDotInverse")
	return
}

func (rt *RTElement) ProjectFunctionOntoDOF(s1, s2 []float64) {
	// For each location in {R,S}, project the input vector function [s1,s2]
	// on to the degrees of freedom of the element
	for j := range s1 {
		r, s := rt.R.AtVec(j), rt.S.AtVec(j)
		b_j := rt.Phi[j].BasisVector
		f := [2]float64{s1[j], s2[j]}
		rt.Projection.Set(j, 0, b_j.Dot(r, s, f))
	}
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
E4:	Φj (r, s) = bj (r, s) * ê4 (r, s) , j = 1, 2, . . . , k(k + 1)/2; {j = 1, NpInt}
E5:	Φj (r, s) = bj (r, s) * ê5 (r, s) , j = 1, 2, . . . , k(k + 1)/2 .

where bj is a 2D polynomial of order P-1
*/
func (rt *RTElement) baseBasisVectors(r, s float64, j int) (ef [2]float64) {
	// Mapping:
	// Edge 1 (Bottom edge)   => Ervin Edge 2
	// Edge 2 (Hypotenuse)  => Ervin Edge 1
	// Edge 3 (Left edge) => Ervin Edge 3
	// Note: Here we fix an error in Ervin - the midpoint of RT0 edges are now
	// correctly unit normal
	// Note2: The Ervin basis is wrong, edge basis vectors should be unit normal
	var (
		sr2 = math.Sqrt2
		// Transform our unit triangle coords to Ervin:
		xi  = 0.5 * (r + 1)
		eta = 0.5 * (s + 1)
	)
	_, _ = xi, eta
	switch rt.getFunctionNumber(j) {
	case E1:
		// ef = [2]float64{xi, eta - 1}
		// ef = [2]float64{xi - 0.5, eta - 1} // Fixed Ervin basis
		ef = [2]float64{0, -1}
	case E2:
		// ef = [2]float64{sr2 * xi, sr2 * eta} // Ervin basis
		ef = [2]float64{0.5 * sr2, 0.5 * sr2}
	case E3:
		// ef = [2]float64{xi - 1, eta}
		// ef = [2]float64{xi - 1, eta - 0.5} // Ervin basis
		ef = [2]float64{-1, 0}
	// case E4:
	// 	ef = [2]float64{1, 0}
	// case E5:
	// 	ef = [2]float64{0, 1}
	case E4:
		ef = [2]float64{eta * xi, eta * (eta - 1)}
	case E5:
		ef = [2]float64{xi * (xi - 1), xi * eta}
	default:
		panic("wrong basis function number (1-5)")
	}
	return
}

func (rt *RTElement) basisEvaluation_new(r, s float64, j int) (psi float64) {
	funcNum := rt.getFunctionNumber(j)
	switch funcNum {
	case E4, E5:
	}

	return
}
func (rt *RTElement) basisPolynomialValue_new(r, s float64, j int,
	derivO ...DerivativeDirection) (val float64) {
	var (
		funcNum = rt.getFunctionNumber(j)
		xi      float64
		jj      int
	)
	if rt.P != 2 {
		panic("wrong polynomial order - testing")
	}

	switch funcNum {
	case E1, E2, E3:
		edgeNum := rt.locationOnEdge(r, s)
		// Do not calculate edge basis functions when not on own edge
		if funcNum != edgeNum {
			val = 0
			return
		}
		jj = j - 2*rt.NpInt
		if edgeNum == E2 {
			jj -= rt.NpEdge
		} else if edgeNum == E3 {
			jj -= 2 * rt.NpEdge
		}
		// val = DG1D.Lagrange1DPoly(xi, Xi, jj, int(deriv))
		val = rt.lp1d.getPolynomial(rt.getEdgeXiParameter(r, s, funcNum), jj,
			derivO...)
	case E4, E5:
		xi = 0.5 * (r + 1)
		eta := 0.5 * (s + 1)
		// Here JJ runs from 0 to 2, and each of the 3 internal points
		// corresponds to a distinct polynomial for each point in Ervin
		// that multiplies either ℯ̂₄ or ℯ̂₅.
		jj = j
		if funcNum == E5 {
			jj -= rt.NpInt
		}
		switch jj {
		case 0:
			val = 1. - xi - eta
		case 1:
			val = xi
		case 2:
			val = eta
		default:
			panic("j out of range")
		}
	}

	return
}

func (rt *RTElement) basisEvaluation(r, s float64, j int) (v [2]float64,
	div float64) {
	// This routine produces the j-th basis vector and divergence at location
	// [r,s]. Note that the Alpha constants are set after using the basis to
	// solve the divergence conformance equation,
	// which sets the constants for each basis vector. For simplicity,
	// we use the Alpha constants set to 1 initially to minimize code
	// duplication.
	e := rt.baseBasisVectors(r, s, j)
	P := rt.basisPolynomialValue(r, s, j)
	// fmt.Printf("[%f,%f]%d P=%f, e=[%f,%f]\n", r, s, j, P, e[0], e[1])
	v = [2]float64{P * e[0], P * e[1]}
	dPdr := rt.basisPolynomialValue(r, s, j, Dr)
	dPds := rt.basisPolynomialValue(r, s, j, Ds)
	dPdXi := dPdr // for edge functions, either of the derivatives will do
	funcNum := rt.getFunctionNumber(j)
	switch funcNum {
	case E4, E5:
		//		// Interior Basis functions
		// The first NpInt index positions in the interior correspond to the E4
		// Ervin basis vectors. The second NpInt index positions correspond to the
		// E5 Ervin basis vectors.
		// The RT1 and RT2 elements are treated differently from the RT3+
		// elements in that they are not multiplied by the same 2D polynomial as
		// RT3+. This is already factored into the support function
		// basisPolynomialValue, so that the divergence can be calculated the same
		// way for all basis functions using the chain rule:
		//           [v]j = P(r,s)j * [E4]; {0     <= j < NpInt}
		//           [v]j = P(r,s)j * [E5]; {NpInt <= j < 2*NpInt}
		//           [v] = [P(r,s)*E4_1, P(r,s)*E4_2]; {0     <= j < NpInt}
		//           [v] = [P(r,s)*E5_1, P(r,s)*E5_2]; {NpInt <= j < 2*NpInt}
		// Conversion from Ervin coordinates:
		// 					xi  = (r + 1)/2
		// 					eta = (s + 1)/2
		// 					r = 2 * xi  - 1
		// 					s = 2 * eta - 1
		//     Div([v]) = [Div]⋅[v] = [d/dr,d/ds]⋅[v1,v2] = dv1/dr + dv2/ds
		if funcNum == E4 {
			// E4 Div:
			//         div = d/dr(P(r,s)*E4_1) + d/ds(P(r,s)*E4_2)
			// div = (dP/dr)*(E4_1) + P*(dE4_1/dr) + (dP/ds)*(E4_2) + P*(dE4_2/ds)
			//	  div = P*(dE4_1/dr+dE4_2/ds) + E4_1*(dP/dr) + E4_2*(dP/ds)
			// 			 [E4] = [eta * xi, eta * (eta - 1)]
			//       = [(s+1)/2 * (r+1)/2   , (s+1)/2 * ((s+1)/2 -1)]
			//       = [(1/4)*(s+1)*(r+1)   , (s+1)/2 * (s-1)/2
			//       = [(1/4)*(s*r+r+s+1)   , (1/4) * (s+1) * (s-1)
			//       = [(1/4)*(s*r +r+s +1) , (1/4) * (s*s - 1)
			//   E4_1 = (1/4)*(s*r +r+s +1) , E4_2 = (1/4)*(s*s -1)
			//          dE4_1/dr = (s+1)/4  , dE4_2/ds = s/2
			//	  div = P*(dE4_1/dr+dE4_2/ds) + E4_1*(dP/dr) + E4_2*(dP/ds)
			//    div = P*((s+1)/4 + s/2) + (dP/dr)*(1/4)*(s*r +r+s +1) +
			//                              (dP/ds)*(1/4)*(s*s -1)
			//    div = (1/4)*P*(3*s+1) + (1/4)*(dP/dr)*(s*r +r+s +1) +
			//                             (1/4)*(dP/ds)*(s*s -1)
			//    div = (1/4)*(P*(3*s+1) + (dP/dr)*(s*r +r+s +1) + (dP/ds)*(s*s-1))
			// div = (1. / 4.) * (P*(3*s+1) + dPdr*(s*r+r+s+1) + dPds*(s*s-1))
			div = dPdr
			return
		} else {
			// E5 Div:
			//         div = d/dr(P(r,s)*E5_1) + d/ds(P(r,s)*E5_2)
			// div = (dP/dr)*(E5_1) + P*(dE5_1/dr) + (dP/ds)*(E5_2) + P*(dE5_2/ds)
			//	  div = P*(dE5_1/dr+dE5_2/ds) + E5_1*(dP/dr) + E5_2*(dP/ds)
			// 			 [E5] = [xi * (xi - 1), eta * xi)]
			//                (1/4)*(r*r -1)  , (s+1)/2 * (r+1)/2]
			//                                , (1/4)*(s+1)*(r+1)]
			//                                , (1/4)*(s*r+r+s+1)]
			//                                , (1/4)*(s*r +r+s +1)]
			//          E5_1 = (1/4)*(r*r -1) , E5_2 = (1/4)*(s*r +r+s +1)]
			//                 dE5_1/dr = r/2 , dE5_2/ds = (r+1)/4
			//	  div = P*(dE5_1/dr+dE5_2/ds) + E5_1*(dP/dr) + E5_2*(dP/ds)
			//    div = P*((r+1)/4 + r/2) + (dP/dr)*(1/4)*(r*r -1) +
			//                              (dP/ds)*(1/4)*(s*r +r+s +1)
			//    div = (1/4)*P*(3*r+1) + (1/4)*(dP/dr)*(r*r -1) +
			//                             (1/4)*(dP/ds)*(s*r +r+s +1)
			// div = (1/4)*(P*(3*r+1) + (dP/dr)*(r*r-1)) + (dP/ds)*(s*r +r+s +1)
			// div = (1. / 4.) * (P*(3*r+1) + dPdr*(r*r-1) + dPds*(s*r+r+s+1))
			div = dPds
			return
		}
	//
	// Edge Divergence and basis vectors
	// the edge basis vector function [v] varies along the edge.
	// It is the product of a 1D edge function f(xi) and [v], so we have:
	// ************************************************************************
	// div(edgeFunction) = div(f(xi)*[v]) =
	//     [df(xi)/dr,df(xi)/ds] ⋅ [v]  +  f(xi) * ([div] ⋅ [v])
	// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
	// ************************************************************************
	//
	case E1:
		// Bottom Edge (1):
		//          [E1] = [xi, eta-1] Ervin Edge 3
		//          [E1] = [xi-1/2, eta-1] Ervin Edge 3 (fixed errata)
		// 			[E1] = [(r+1)/2-1/2, (s+1)/2-1]
		// 			[E1] = [r/2, (s+1)/2-1]
		// 			[E1] = [r/2, (s-1)/2] (fixed errata)
		// E1_1 = r/2, E1_2 = (s-1)/2, dE1_1/dr = 1/2, dE1_2/ds = 1/2
		// The bottom edge in a counter-clockwise direction is parameterized:
		// xi = r, s = -1 (const) => dxi/dr = 1, dxi/ds = 0
		// div = df/dxi*(E1_1*(dxi/dr)+E1_2*(dxi/ds)) + f(xi)*(dv1/dr + dv2/ds)
		//     = df/dxi*(E1_1*(  1   )+E1_2*(   0  )) + f(xi)*( 1/2  +   1/2 )
		//     = df/dxi*(         E1_1          ) + f(xi)
		//        div(edge1) = f(xi) + v1 * (df/dxi)
		// div = P + e[0]*dPdXi // Ervin basis
		//          [E1] = [0, -1]  Bottom edge
		//           div = df/dxi * (0 -1) = -dfdxi
		div = -0.5 * dPdXi
		return
	case E2:
		// Hypotenuse (2) divergence:
		//          [E2] = [Sqrt2 * xi, Sqrt2 * eta] Ervin Edge 1
		// 			[E2] = [Sqrt2/2 * (r+1), Sqrt2/2 * (s+1)]
		// v1 = Sqrt2/2 * (r+1), v2 = Sqrt2/2 * (s+1), dv1/dr = Sqrt2/2 = dv2/ds
		// The hypotenuse in a counter-clockwise direction is parameterized:
		// xi = s, => dxi/dr = 0, dxi/ds = 1
		// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
		//     = df/dxi*(v1*(   0  )+v2*(   1  )) + f(xi) * (Sqrt2/2+Sqrt2/2)
		//     = df/dxi*(          v2           ) + f(xi) * Sqrt2
		//         div(edge2) = Sqrt2 * f(xi) + (v2-v1) * (df/dxi)
		// div = 0.5*math.Sqrt2*P + e[1]*dPdXi // Ervin basis
		// 			[E2] = [Sqrt2/2, Sqrt2/2]
		//           div = df/dxi * (Sqrt2/2 + Sqrt2/2) = Sqrt2 * df/dxi
		div = math.Sqrt2 * dPdXi
		return
	case E3:
		// Left Edge (3) divergence:
		//          [E3] = [xi-1, eta] Ervin Edge 2
		//          [E3] = [xi-1, eta-1/2] Ervin Edge 2 (fixed errata)
		// 			[E3] = [(r+1)/2 - 1, (s+1)/2-1/2]
		// 			[E3] = [(r-1)/2, s/2] (fixed errata)
		// v1 = (r-1)/2, v2 = s/2, dv1/dr = 1/2, dv2/ds = 1/2
		// The left edge  in a counter-clockwise direction is parameterized:
		// r = -1 (constant), xi = -s => dxi/dr = 0, dxi/ds = -1,
		// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
		//     = df/dxi*(v1*(  0   )+v2*(  -1  )) + f(xi) * (  1/2  +   1/2 )
		//     = df/dxi*(         -v2           ) + f(xi)
		//         div(edge3) = f(xi) - v2 * (df/dxi)
		// div = P - e[1]*dPdXi // Ervin Basis
		//			[E3] = [-1,0]
		//           div = df/dxi * (-1 + 0) = -df/dxi
		div = -0.5 * dPdXi
		return
	}
	return
}

func (rt *RTElement) getEdgeXiParameter(r, s float64,
	funcNum RTFunctionNumber) (xi float64) {
	// Edge 1 is S=-1 (bottom of tri)
	// Edge 2 is Hypotenuse
	// Edge 3 is R=-1 (left side of tri)
	switch funcNum {
	case E1:
		xi = r
	case E2:
		xi = s
	case E3:
		xi = -s
	default:
		panic("invalid edgeNum")
	}
	return
}

func (rt *RTElement) getFunctionNumber(j int) (funcNum RTFunctionNumber) {
	var (
		NpInt  = rt.NpInt
		NpEdge = rt.NpEdge
	)
	switch {
	case j >= 0 && j < NpInt:
		// Unit vector is [1,0]
		funcNum = E4
	case j >= NpInt && j < 2*NpInt:
		// Unit vector is [0,1]
		funcNum = E5
	case j >= 2*NpInt && j < 2*NpInt+NpEdge:
		// E1: Unit vector is [0,-1]
		funcNum = E1
	case j >= 2*NpInt+NpEdge && j < 2*NpInt+2*NpEdge:
		// E2: Unit vector is [1/sqrt(2), 1/sqrt(2)]
		funcNum = E2
	case j >= 2*NpInt+2*NpEdge && j < 2*NpInt+3*NpEdge:
		// E3: Unit vector is [-1,0]
		funcNum = E3
	default:
		panic("j out of range")
	}
	return
}

func (rt *RTElement) locationOnEdge(r, s float64) RTFunctionNumber {
	var (
		eps = 0.000001
	)
	switch {
	case math.Abs(s+1) < eps: // s == -1, on bottom edge
		return E1
	case math.Abs(r+s) < eps: // Hypotenuse
		return E2
	case math.Abs(r+1) < eps: // r == -1, on Left edge
		return E3
	}
	return All
}

func (rt *RTElement) basisPolynomialValue(r, s float64, j int,
	derivO ...DerivativeDirection) (val float64) {
	// This evaluates the j-th polynomial or derivative at r,s used to multiply
	// each of the 5 base basis vectors in Ervin's RT basis, named
	//               ℯ̂₁, ℯ̂₂, ℯ̂₃, ℯ̂₄, ℯ̂₅
	// The first three are the edge basis vectors, but in Ervin's paper
	// e1 is the hypotenuse, e2 is left, e3 is bottom, where here
	// edge 1 is bottom, edge 2 is hypotenuse, edge 3 is left
	//
	// This function provides the value of the polynomial that multiplies
	// each basis vector, not including the basis vector values or derivatives.
	// Note that to compute the derivatives of the full basis including the
	// base vector and polynomial, you have to use the chain rule. This function
	// is a helper for that, in that it provides the value and derivative of the
	// multiplying polynomial for use in calculating the basis or derivatives.
	//
	// The input parameter "j" is the function index within the element
	var (
		funcNum = rt.getFunctionNumber(j)
		xi      float64
		jj      int
	)

	switch funcNum {
	case E1, E2, E3:
		edgeNum := rt.locationOnEdge(r, s)
		// Do not calculate edge basis functions when not on own edge
		if funcNum != edgeNum {
			val = 0
			return
		}
		jj = j - 2*rt.NpInt
		if edgeNum == E2 {
			jj -= rt.NpEdge
		} else if edgeNum == E3 {
			jj -= 2 * rt.NpEdge
		}
		// val = DG1D.Lagrange1DPoly(xi, Xi, jj, int(deriv))
		val = rt.lp1d.getPolynomial(rt.getEdgeXiParameter(r, s, funcNum), jj,
			derivO...)
		// fmt.Printf("Edge[%s], Val[%f]jj=%d =%f, Deriv = %v\n",
		// 	funcNum.String(), xi, jj, val, deriv)
	case E4, E5:
		// Don't calculate interior basis on edges
		switch rt.locationOnEdge(r, s) {
		case E1, E2, E3:
			val = 0
			return
		}
		// We represent the interior polynomial for the special cases of the
		// RT1 and RT2 elements. For RT3 and above, we use a Jacobi2D
		// polynomial, but for RT1 it's a constant and for RT2 it's a basic poly set
		var derivA DerivativeDirection
		if len(derivO) > 0 {
			derivA = derivO[0]
		}
		switch {
		case rt.P == 0:
			panic("RT element isn't defined here for RT0")
		case rt.P == 1:
			// No polynomial multiplies the basis vectors - constant
			if derivA == None {
				val = 1
			} else {
				val = 0
			}
			return
		case rt.P == 2:
			switch derivA {
			case None:
				xi = 0.5 * (r + 1)
				eta := 0.5 * (s + 1)
				// Here JJ runs from 0 to 2, and each of the 3 internal points
				// corresponds to a distinct polynomial for each point in Ervin
				// that multiplies either ℯ̂₄ or ℯ̂₅.
				jj = j
				if funcNum == E5 {
					jj -= rt.NpInt
				}
				switch jj {
				case 0:
					val = 1. - xi - eta
				case 1:
					val = xi
				case 2:
					val = eta
				default:
					panic("j out of range")
				}
			case Dr:
				switch jj {
				case 0:
					val = -0.5
				case 1:
					val = 0.5
				case 2:
					val = 0
				default:
					panic("j out of range")
				}
			case Ds:
				switch jj {
				case 0:
					val = -0.5
				case 1:
					val = 0
				case 2:
					val = 0.5
				default:
					panic("j out of range")
				}
			}
		case rt.P >= 3: // RT Element Polynomial Order is RT3 or greater
			jj = j
			if funcNum == E5 {
				jj -= rt.NpInt
			}
			val = rt.lp2d.GetPolynomialEvaluation(r, s, jj, derivO...)
		default:
			panic("RT polynomial degree out of range")
		}
	default:
		panic("j basis function number (1-5)")
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
