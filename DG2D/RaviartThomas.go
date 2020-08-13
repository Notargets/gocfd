package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type RTPointType uint

const (
	InteriorR RTPointType = iota
	InteriorS
	Edge1
	Edge2
	Edge3
)

func GetPointType(N int, i int) (rtt RTPointType) {
	var (
		NInterior = N * (N + 1) / 2 // one order less than RT element in (P_k)2
	)
	switch {
	case i < NInterior:
		// Unit vector is [1,0]
		rtt = InteriorR
	case i >= NInterior && i < 2*NInterior:
		// Unit vector is [0,1]
		rtt = InteriorS
	case i >= 2*NInterior && i < 2*NInterior+(N+1):
		// Edge1
		rtt = Edge1
	case i >= 2*NInterior+(N+1) && i < 2*NInterior+2*(N+1):
		// Edge2: Unit vector is [-1,0]
		rtt = Edge2
	case i >= 2*NInterior+2*(N+1) && i < 2*NInterior+3*(N+1):
		// Edge3: // Unit vector is [0,-1]
		rtt = Edge3
	}
	return
}

func RTBasis(N int, R, S utils.Vector) (V, Vinv utils.Matrix) {
	/*
				This is constructed from the defining space of the RT element:
								 2
					RT_k = [(P_k) ]   + [ X ] P_k
						 = [ b1(r,s)_i + r * b3(r,s)_j ]
						   [ b2(r,s)_i + s * b3(r,s)_j ]
					i := 1, (K+1)(K+2)/2
					j := 1, (K+1) (highest order terms in polynomial)

				The dimension of RT_k is (K+1)(K+3) and we can see from the above that the total
				number of terms in the polynomial will be:
					  (K+1)(K+2) + K+1 = (K+1)(K+3)

				The explanation for why the b3 polynomial sub-basis is partially consumed:
				When multiplied by [ X ], the b3 polynomial produces terms redundant with
				the b1 and b2 sub-bases. The redundancy is removed from the b3 sub-basis to
				compensate, producing the correct overall dimension.

				Another more physical way to think about it: The [ X ] * P_k term is tracing a
				1D shell within the 2D space - the [ X ] "pointer vector" is permuted through
				a 1D polynomial at high order to represent the outer surface of the shape.

				Two groups of bases, two ways, a pair of polynomial types and a pair of geometric types:
					1) The RT basis consists of two parts, the (P_k)2 basis and a P_k basis with (N+1) terms.
					The dimension of the first part is 2 * (K+1)(K+2)/2 and the second is (K+1). The total number of terms
					in the polynomial space is:
						2*(K+1)(K+2)/2 + (K+1) = (K+3)(K+1)

					2) The RT basis is also composed of two types of geometric bases, interior and exterior
					points. The number of interior points is (K)(K+1)/2, and the number of edge points is 3*(K+1).
					For each interior point, we have two basis vectors, [1,0] and [0,1]. The total degrees of freedom are:
						2*(K)(K+1)/2 + 3(K+1) = (K+3)*(K+1)

				The number of interior points matches a 2D Lagrangian element basis at order (K-1):
		    		There are (K)(K+1)/2 interior points in this RT element, which matches the point count of a Lagrangian
					element at order (K-1). This is very convenient and enabling for the DFR method, as it allows us to
					represent the flux vector function of a solution at degree (K-1) on an RT element of order (K) by
					simply transferring the values from the (K-1) solution element to the interior of the RT(K) element.
					We then provide the flux values along the triangle edges of the RT(K) element, after which we can
					calculate gradient, divergence, and curl using a polynomial of degree (K), yielding a gradient,
					divergence, curl of order (K-1), which is exactly what we need for the solution at (K-1).
	*/
	/*	   			Inputs:
							(N)(N+2)/2 [r,s] points from the interior of the [-1,1] triangle

		    		Outputs:
						[R,S]: Coordinates of points within element
				    	First (N)(N+1)/2 points: Interior points, excluding edges in [-1,1] triangle coordinates
				    	Next (N)(N+1)/2 points: Duplicate of above
						Next (N+1) points: Edge 1 locations
						Next (N+1) points: Edge 2 locations
						Next (N+1) points: Edge 3 locations

						[V1, V2]: Vandermonde matrix for each of R and S directions:
	    				V_ij = Psi_j([r,s]_i)
						Element V_ij of the Vandermonde matrix is the basis function Psi_j evaluated at [r,s]_i
						Since this is a vector valued element/basis, we have a interpolation matrix for each direction.
	*/
	var (
		err    error
		Np     = (N + 1) * (N + 3)
		p1, p2 = make([]float64, Np), make([]float64, Np)
	)
	/*
		Transform the incoming (r,s) locations into the unit triangle (0,1) from the (-1,1) basis
	*/
	R.AddScalar(1).Scale(0.5)
	S.AddScalar(1).Scale(0.5)
	// Add the edge points from the RT basis
	R, S = AddEdgePoints(N, R, S)
	/*
		First, evaluate the polynomial at the (r,s) coordinates
		This is the same set that will be used for all dot products to form the basis matrix
	*/
	Term2D := func(r, s float64, i, j int) (val float64) {
		// Note this only outputs the value of the (i,j)th term of the 2D polynomial
		val = utils.POW(r, j) * utils.POW(s, i)
		return
	}

	/*
		Form the basis matrix by forming a dot product with unit vectors, matching the coordinate locations in R,S
	*/
	V = utils.NewMatrix(Np, Np)
	rowEdge := make([]float64, Np)
	oosr2 := 1 / math.Sqrt(2)

	// Evaluate at geometric locations
	N2DBasis := (N + 1) * (N + 2) / 2 // Number of polynomial terms for each of R and S directions
	for ii, rr := range R.Data() {
		ss := S.Data()[ii]
		// Evaluate the full 2D polynomial basis first, once for each of two components
		var sk int
		for i := 0; i <= N; i++ {
			for j := 0; j <= (N - i); j++ {
				val := Term2D(rr, ss, i, j)
				p1[sk] = val
				p2[sk+N2DBasis] = val
				sk++
			}
		}
		sk += N2DBasis // Skip to the beginning of the second polynomial group
		// Evaluate the term ([ X ]*(Pk)) at only the top N+1 terms (highest order) of the 2D polynomial
		for i := 0; i <= N; i++ {
			j := N - i
			val := Term2D(rr, ss, i, j)
			p1[sk] = val * rr
			p2[sk] = val * ss
			sk++
		}
		switch GetPointType(N, ii) {
		case InteriorR:
			// Unit vector is [1,0]
			V.M.SetRow(ii, p1)
		case InteriorS:
			// Unit vector is [0,1]
			V.M.SetRow(ii, p2)
		case Edge1:
			for i := range rowEdge {
				// Edge1: Unit vector is [1/sqrt(2), 1/sqrt(2)]
				rowEdge[i] = oosr2 * (p1[i] + p2[i])
			}
			V.M.SetRow(ii, rowEdge)
		case Edge2:
			for i := range rowEdge {
				// Edge2: Unit vector is [-1,0]
				rowEdge[i] = -p1[i]
			}
			V.M.SetRow(ii, rowEdge)
		case Edge3:
			for i := range rowEdge {
				// Edge3: // Unit vector is [0,-1]
				rowEdge[i] = -p2[i]
			}
			V.M.SetRow(ii, rowEdge)
		}
	}
	if Vinv, err = V.Inverse(); err != nil {
		panic(err)
	}
	return V, Vinv
}

func AddEdgePoints(N int, rInt, sInt utils.Vector) (r, s utils.Vector) {
	var (
		NpEdge       = N + 1
		rData, sData = rInt.Data(), sInt.Data()
	)
	/*
		Determine geometric locations of edge points, located at Gauss locations in 1D, projected onto the edges
	*/
	GQR, _ := DG1D.JacobiGQ(1, 1, N)
	// Transform into the space (0,1) from (-1,1)
	GQR.AddScalar(1).Scale(0.5)

	/*
		Double the number of interior points to match each direction of the basis
	*/
	rData = append(rData, rData...)
	sData = append(sData, sData...)

	// Calculate the triangle edges
	GQRData := GQR.Data()
	rEdgeData := make([]float64, NpEdge*3)
	sEdgeData := make([]float64, NpEdge*3)
	for i := 0; i < NpEdge; i++ {
		// Edge 1 (hypotenuse)
		rEdgeData[i] = 1 - GQRData[i]
		sEdgeData[i] = GQRData[i]
		// Edge 2
		rEdgeData[i+NpEdge] = 0
		sEdgeData[i+NpEdge] = 1 - GQRData[i]
		// Edge 3
		rEdgeData[i+2*NpEdge] = GQRData[i]
		sEdgeData[i+2*NpEdge] = 0
	}
	rData = append(rData, rEdgeData...)
	sData = append(sData, sEdgeData...)
	r = utils.NewVector(len(rData), rData)
	s = utils.NewVector(len(sData), sData)
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
		Np   = (N + 1) * (N + 2) / 2
		epsD []float64
	)
	switch N {
	case 0:
		R = utils.NewVector(1, []float64{-.5})
		S = utils.NewVector(1, []float64{-.5})
		return
	case 3:
		epsD = []float64{
			0.3333333333333333, 0.055758983558155, 0.88848203288369, 0.055758983558155, 0.290285227512689, 0.6388573870878149, 0.290285227512689, 0.6388573870878149, 0.070857385399496, 0.070857385399496,
			0.3333333333333333, 0.055758983558155, 0.055758983558155, 0.88848203288369, 0.070857385399496, 0.290285227512689, 0.6388573870878149, 0.070857385399496, 0.290285227512689, 0.6388573870878149,
			0.3333333333333333, 0.88848203288369, 0.055758983558155, 0.055758983558155, 0.6388573870878149, 0.070857385399496, 0.070857385399496, 0.290285227512689, 0.6388573870878149, 0.290285227512689,
		}
	case 4:
		epsD = []float64{
			0.034681580220044, 0.9306368395599121, 0.034681580220044, 0.243071555674492, 0.513856888651016, 0.243071555674492, 0.473372556704605, 0.05325488659079003, 0.473372556704605, 0.200039998995093, 0.752666332493468, 0.200039998995093, 0.752666332493468, 0.047293668511439, 0.047293668511439,
			0.034681580220044, 0.034681580220044, 0.9306368395599121, 0.243071555674492, 0.243071555674492, 0.513856888651016, 0.473372556704605, 0.473372556704605, 0.05325488659079003, 0.047293668511439, 0.200039998995093, 0.752666332493468, 0.047293668511439, 0.200039998995093, 0.752666332493468,
			0.9306368395599121, 0.034681580220044, 0.034681580220044, 0.513856888651016, 0.243071555674492, 0.243071555674492, 0.05325488659079003, 0.473372556704605, 0.473372556704605, 0.752666332493468, 0.047293668511439, 0.047293668511439, 0.200039998995093, 0.752666332493468, 0.200039998995093,
		}
	default:
		panic(fmt.Errorf("Epsilon nodes not defined for N = %v, only defined for N=3 or N=4\n", N))
	}
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
