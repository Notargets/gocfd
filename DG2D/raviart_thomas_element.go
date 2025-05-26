package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
)

/*
	Definition:
		Raviart Thomas element
		Total of (N+1)(N+3) degrees of freedom. (N)(N+1) interior degrees of
		freedom split between R and S directions. 3*(N+1) edge degrees of
		freedom, (N+1) per each of 3 edges on the triangle.

	Inputs:
		(N)(N+2)/2 [R,S] points from the interior of the [-1,1] triangle

	Outputs:
		[R,S]: Coordinates of points within element

		Each point in this basis corresponds to a basis vector with 2 components
							e.g. [v]j = [v1,v2]j

		Interior points in [R,S] each have two basis vectors and each basis
		vector has two components (it'S a vector!). That means the interior
		point indices are duplicated for the interior points, because we need
		to refer to each basis vector independently.

		Edge points in [R,S] each have one basis vector, and each basis vector
		has two components (it'S a vector!).

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

type RTBasisFunctionType uint8

const (
	All RTBasisFunctionType = iota
	E4                      // Ervin E4 Interior basis
	E5                      // Ervin E5 Interior basis
	E1                      // Bottom Edge, Normal is [0,-1], Ervin E3
	E2                      // Hypotenuse, Normal is [1/Sqrt2,1/Sqrt2], Ervin E1
	E3                      // Left Edge, Normal is [-1,0], Ervin E2
)

func (fn RTBasisFunctionType) String() string {
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

type RTBasisType uint8

const (
	ErvinBasisRT = RTBasisType(iota)
	SimplexRTBasis
)

func (fn RTBasisType) String() string {
	switch fn {
	case ErvinBasisRT:
		return "ErvinBasisRT"
	case SimplexRTBasis:
		return "SimplexRTBasis"
	}
	return "Unknown"
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
	V, VInv     utils.Matrix // Vandermonde Matrix
	Div, DivInt utils.Matrix // divergence matrix, NpxNp for all,
	// NintxNp Interior Points
	R, S       utils.Vector // Point locations defining element in [-1,1] Triangle, NpxNp
	RInt, SInt utils.Vector
	DOFVectors []*ConstantVector // The unit vectors for each DOF
	RTBasis    RTBasisType
	Phi        []VectorI // Each term of the basis is a vector
}

func NewRTElement(P int, basisType RTBasisType) (rt *RTElement) {
	// We expect that there are points in R and S to match the dimension of dim(P(NFlux-1))
	/*
		<---- NpInt ----><---- NpInt ----><---NpEdge----><---NpEdge----><---NpEdge---->
		         Solution Points          Edge 1 pts	Edge 2 pts	  Edge 3 pts
		<---- NpInt ----><---- NpInt ----><---NpEdge----><---NpEdge----><---NpEdge---->
	*/
	var (
		Np     = (P + 1) * (P + 3)
		NpInt  = P * (P + 1) / 2 // Number of interior points is same as the 2D scalar space one order lesser
		NpEdge = P + 1           // Each edge is P+1 nodes
		oosr2  = 0.5 * math.Sqrt2
	)
	rt = &RTElement{
		P:          P,
		Np:         Np,
		NpInt:      NpInt,
		NpEdge:     NpEdge,
		V:          utils.NewMatrix(Np, Np),
		Div:        utils.NewMatrix(Np, Np),
		DOFVectors: make([]*ConstantVector, Np),
		RTBasis:    basisType,
		DivInt:     utils.NewMatrix(NpInt, Np),
	}

	if P < 1 {
		panic("P must be greater than or equal to 1")
	}
	rt.RInt, rt.SInt = MakeRSFromPoints(WilliamsShunnJameson(P - 1))

	// Construct the unit vectors for the DOFs

	for i := 0; i < NpInt; i++ {
		rt.DOFVectors[i] = NewConstantVector(1, 0)
		rt.DOFVectors[i+NpInt] = NewConstantVector(0, 1)
	}

	offset := 2 * NpInt
	for i := 0; i < NpEdge; i++ {
		rt.DOFVectors[offset+i] = NewConstantVector(0, -1)
		rt.DOFVectors[offset+i+NpEdge] = NewConstantVector(oosr2, oosr2)
		rt.DOFVectors[offset+i+2*NpEdge] = NewConstantVector(-1, 0)
	}

	rt.R, rt.S = rt.ExtendGeomToRT(rt.RInt, rt.SInt)
	rt.CalculateBasis()

	// Compose basis Vandermonde matrix
	rt.V = rt.ComposeV(rt.Phi)
	rt.VInv = rt.V.InverseWithCheck()
	rt.Div = rt.ComputeDivergenceMatrix()
	for i := 0; i < NpInt; i++ {
		for j := 0; j < Np; j++ {
			rt.DivInt.Set(i, j, rt.Div.At(i, j))
		}
	}
	return
}

func (rt *RTElement) CalculateBasis() {
	switch rt.RTBasis {
	case ErvinBasisRT:
		e := NewErvinRTBasis(rt.P, rt.R, rt.S)
		rt.Phi = e.Phi
	case SimplexRTBasis:
		rjb := NewRTBasisSimplex(rt.P, rt.R, rt.S)
		rt.Phi = rjb.Phi
	default:
		panic("No basis chosen")
	}
	return
}

func (rt *RTElement) ComposeV(Phi []VectorI) (V utils.Matrix) {
	// Computes a Vandermonde matrix from a set of VectorFunction and element
	// location vectors
	var (
		Np   = rt.Np
		R, S = rt.R, rt.S
	)
	// fmt.Printf("Length of Phi:%d\n", len(rt.Phi))
	V = utils.NewMatrix(Np, Np)
	for i := 0; i < Np; i++ {
		r_i, s_i := R.DataP[i], S.DataP[i]
		b_i := rt.DOFVectors[i]
		for j := 0; j < Np; j++ {
			V.Set(i, j, Phi[j].Dot(r_i, s_i, b_i.Eval()))
		}
	}
	return
}

func (rt *RTElement) ComputeDivergenceMatrix() (Div utils.Matrix) {
	// Build basis divergence matrix
	DivBasis := utils.NewMatrix(rt.Np, rt.Np)
	for i := 0; i < rt.Np; i++ {
		r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		for j := 0; j < rt.Np; j++ {
			// fmt.Printf("NpInt, j = %d, %d\n", rt.NpInt, j)
			DivBasis.Set(i, j, rt.Phi[j].Divergence(r, s))
		}
	}
	// DivBasis.String("Basis divergence Matrix")
	Div = DivBasis.Mul(rt.VInv)
	return
}

func (rt *RTElement) ProjectFunctionOntoDOF(s1, s2, FProj []float64) {
	// For each location in {R,S}, project the input vector function [s1,s2]
	// on to the degrees of freedom of the element
	// FProj should already be allocated to size [Np,1]
	for j := range s1 {
		b_j := rt.DOFVectors[j]
		f := [2]float64{s1[j], s2[j]}
		FProj[j] = b_j.Dot(f)
	}
	return
}

func (rt *RTElement) getLocationType(j int) (funcNum RTBasisFunctionType) {
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

func (rt *RTElement) locationOnEdge(r, s float64) RTBasisFunctionType {
	var (
		eps = 0.000001
	)
	switch {
	case math.Abs(s+1) < eps: // S == -1, on bottom edge
		return E1
	case math.Abs(r+s) < eps: // Hypotenuse
		return E2
	case math.Abs(r+1) < eps: // R == -1, on Left edge
		return E3
	}
	return All
}

func (rt *RTElement) GetEdgeLocations(F []float64) (
	Fedge []float64) {
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

func (rt *RTElement) GetInternalLocations(F []float64) (
	Finternal []float64) {
	var (
		Nint = rt.NpInt
	)
	Finternal = make([]float64, Nint)
	for i := 0; i < Nint; i++ {
		Finternal[i] = F[i]
	}
	return
}
func GetOptimizedEdgePointsEpsilon(NRT int) (Rdist []float64) {
	// Use optimized edge points from edge_point_distribution optimization
	switch NRT {
	case 1:
		Rdist = []float64{-0.38490018, 0.38490018}
	case 2:
		Rdist = []float64{-0.58378055, 0, 0.58378055}
	case 3:
		Rdist = []float64{-0.70460764, -0.25640858, 0.25640858, 0.70460764}
	case 4:
		Rdist = []float64{-0.78093647, -0.43143958, 0, 0.43143071, 0.780937}
	case 5:
		Rdist = []float64{-0.83155051, -0.55394627, -0.19423772, 0.19420846,
			0.55372084, 0.83157178}
	case 6:
		Rdist = []float64{-0.86736381, -0.64201719, -0.34118457, 0,
			0.34103526, 0.64187617, 0.86666209}
	case 7:
		Rdist = []float64{-0.89215475, -0.7065447, -0.45294463, -0.15588555,
			0.15590652, 0.45318727, 0.70679979, 0.89196653}
	case 8:
		Rdist = []float64{-0.91090107, -0.755625, -0.5402, -0.28137625,
			0, 0.28119241, 0.54000096, 0.75548225, 0.91083652}
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
	var Rdist []float64
	var GQR utils.Vector
	// switch nodeType {
	// case Hesthaven, Uniform, WSJ:
	// case Hesthaven, Uniform:
	// 	GQR = utils.NewVector(N+1, DG1D.LegendreZeros(N))
	// case WSJ, Epsilon:
	Rdist = GetOptimizedEdgePointsEpsilon(N)
	GQR = utils.NewVector(N+1, Rdist)
	// case Epsilon:
	// 	Use optimized edge points from edge_point_distribution optimization
	// Rdist = GetOptimizedEdgePointsEpsilon(N)
	// GQR = utils.NewVector(N+1, Rdist)
	// }
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
		Np    = (N + 1) * (N + 2) / 2
		epsD  []float64
		types = []string{"Equidistant", "Williams and Shun",
			"Romero and Jameson(2017)"}
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
	_ = nameType
	// fmt.Printf("Node distribution type is: [%S]\n", nameType)
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
