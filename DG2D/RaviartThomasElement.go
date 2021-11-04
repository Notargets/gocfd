package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"

	"github.com/notargets/gocfd/utils"
)

type RTElement struct {
	N           int             // Order of element
	Np          int             // Number of points in element
	Nedge, Nint int             // Number of Edge and Interior points
	A           utils.Matrix    // Polynomial coefficient matrix, NpxNp
	V           [2]utils.Matrix // Vandermonde matrix for each direction r and s, [2]xNpxNp
	Div, DivInt utils.Matrix    // Divergence matrix, NpxNp for all, NintxNp Interior Points
	R, S        utils.Vector    // Point locations defining element in [-1,1] Triangle, NpxNp
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

func NewRTElement(R, S utils.Vector, N int, useLagrangeBasis bool) (rt *RTElement) {
	// We expect that there are points in R and S to match the dimension of dim(P(N-1))
	/*
		<---- Nint ----><---- Nint ----><---Nedge----><---Nedge----><---Nedge---->
		         Solution Points          Edge 1 pts	Edge 2 pts	  Edge 3 pts
		<---- Nint ----><---- Nint ----><---Nedge----><---Nedge----><---Nedge---->
	*/
	var (
		NN         = N - 1
		NpInterior = (NN + 1) * (NN + 2) / 2
	)
	if R.Len() != NpInterior || S.Len() != NpInterior {
		panic("incorrect number of interior points supplied")
	}
	rt = &RTElement{
		N:     N,
		R:     R,
		S:     S,
		Nint:  N * (N + 1) / 2,
		Nedge: N + 1,
	}
	rt.CalculateBasis(useLagrangeBasis)
	return
}

func (rt *RTElement) ProjectFunctionOntoBasis(s1, s2 []float64) (sp []float64) {
	var (
		Np = len(s1)
	)
	sp = make([]float64, Np)
	oosr2 := 1 / math.Sqrt(2)
	for i := range s1 {
		switch rt.GetTermType(i) {
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

func (rt *RTElement) GetTermType(i int) (rtt RTPointType) {
	var (
		N     = rt.N
		Nint  = N * (N + 1) / 2 // one order less than RT element in (P_k)2
		Nedge = (N + 1)
	)
	switch {
	case i < Nint:
		// Unit vector is [1,0]
		rtt = InteriorR
	case i >= Nint && i < 2*Nint:
		// Unit vector is [0,1]
		rtt = InteriorS
	case i >= 2*Nint && i < 2*Nint+Nedge:
		// Edge1: Unit vector is [0,-1]
		rtt = Edge1
	case i >= 2*Nint+Nedge && i < 2*Nint+2*Nedge:
		// Edge2: Unit vector is [1/sqrt(2), 1/sqrt(2)]
		rtt = Edge2
	case i >= 2*Nint+2*Nedge && i < 2*Nint+3*Nedge:
		// Edge3: Unit vector is [-1,0]
		rtt = Edge3
	}
	return
}

func (rt *RTElement) CalculateBasis(useLagrangeBasis bool) {
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
		                    b1(r,s)        b2(r,s)     b3(r,s)
		                  (K+1)(K+2)/2 + (K+1)(K+2)/2 + (K+1)
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
		N    = rt.N
		R, S = rt.R, rt.S
		Np   = (N + 1) * (N + 3)
		Nint = N * (N + 1) / 2
		P    utils.Matrix
		L2Dp *LagrangeBasis2D // 2D Lagrange Polynomial Basis
	)
	// Add the edge and additional interior (duplicated) points to complete the RT geometry2D
	rt.R, rt.S = ExtendGeomToRT(N, R, S)
	if useLagrangeBasis {
		var RR, SS utils.Vector
		if N < 8 {
			RR, SS = NodesEpsilon(N)
		} else {
			RR, SS = XYtoRS(Nodes2D(N))
		}
		L2Dp = NewLagrangeBasis2D(N, RR, SS)
	}
	/*
		Form the basis matrix by forming a dot product with unit vectors, matching the coordinate locations in R,S
	*/
	P = utils.NewMatrix(Np, Np)
	rowEdge := make([]float64, Np)
	oosr2 := 1 / math.Sqrt(2)

	// Evaluate at geometric locations
	var p0, p1 []float64
	for ii, rr := range rt.R.DataP {
		ss := rt.S.DataP[ii]
		/*
			First, evaluate the polynomial at the (r,s) coordinates
			This is the same set that will be used for all dot products to form the basis matrix
		*/
		if useLagrangeBasis {
			p0, p1 = rt.EvaluateLagrangeRTBasis(L2Dp, rr, ss) // each of p1,p2 stores the polynomial terms for the R and S directions
		} else {
			p0, p1 = rt.EvaluateONRTBasis(rr, ss) // each of p1,p2 stores the polynomial terms for the R and S directions
		}
		// Implement dot product of (unit vector)_ii with each vector term in the polynomial evaluated at location ii
		switch rt.GetTermType(ii) {
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
	// Invert [P] = [A] to obtain the coefficients (columns) of polynomials (rows), each row is a polynomial
	rt.A = P.InverseWithCheck()
	// Evaluate 2D polynomial basis at geometric locations, also evaluate derivatives Dr and Ds for R and S
	P0, P1 := utils.NewMatrix(Np, Np), utils.NewMatrix(Np, Np)
	Pdr0, Pds1 := utils.NewMatrix(Np, Np), utils.NewMatrix(Np, Np)
	for ii, rr := range rt.R.DataP {
		ss := rt.S.DataP[ii]
		if useLagrangeBasis {
			p0, p1 = rt.EvaluateLagrangeRTBasis(L2Dp, rr, ss) // each of p1,p2 stores the polynomial terms for the R and S directions
		} else {
			p0, p1 = rt.EvaluateONRTBasis(rr, ss) // each of p0,p1 stores the polynomial terms for the R and S directions
		}
		P0.M.SetRow(ii, p0)
		P1.M.SetRow(ii, p1)
		if useLagrangeBasis {
			p0, _ = rt.EvaluateLagrangeRTBasis(L2Dp, rr, ss, Dr) // each of p0,p1 stores the polynomial terms for the R and S directions
			_, p1 = rt.EvaluateLagrangeRTBasis(L2Dp, rr, ss, Ds) // each of p0,p1 stores the polynomial terms for the R and S directions
		} else {
			p0, _ = rt.EvaluateONRTBasis(rr, ss, Dr) // each of p0,p1 stores the polynomial terms for the R and S directions
			_, p1 = rt.EvaluateONRTBasis(rr, ss, Ds) // each of p0,p1 stores the polynomial terms for the R and S directions
		}
		Pdr0.M.SetRow(ii, p0)
		Pds1.M.SetRow(ii, p1)
	}
	// Construct the Vandermonde matrices for each direction by multiplying coefficients of constrained basis
	rt.V[0] = P0.Mul(rt.A)
	rt.V[1] = P1.Mul(rt.A)
	rt.Div = Pdr0.Mul(rt.A).Add(Pds1.Mul(rt.A))
	rt.DivInt = utils.NewMatrix(Nint, Np)
	for i := 0; i < Nint; i++ {
		rt.DivInt.M.SetRow(i, rt.Div.Row(i).DataP)
	}
	rt.Np = Np
	return
}

type DerivativeDirection uint8

const (
	None DerivativeDirection = iota
	Dr
	Ds
)

func (rt *RTElement) EvaluateONRTBasis(r, s float64, derivO ...DerivativeDirection) (p1, p2 []float64) {
	var (
		sk       int
		N        = rt.N
		Np       = (N + 1) * (N + 3)
		N2DBasis = (N + 1) * (N + 2) / 2 // Number of polynomial terms for each of R and S directions
		deriv    = None
		tFunc    func(r, s float64, i, j int) (val float64)
	)
	if len(derivO) != 0 {
		deriv = derivO[0]
	}
	ONTerm2D := func(r, s float64, i, j int) (val float64) {
		return Simplex2DPTerm(r, s, i, j)
	}
	DrONTerm2D := func(r, s float64, i, j int) (val float64) {
		val, _ = GradSimplex2DPTerm(r, s, i, j)
		//fmt.Printf("Dr r,s,i,j,val = %8.5f,%8.5f,%d,%d,%8.5f,", r, s, i, j, val)
		return
	}
	DsONTerm2D := func(r, s float64, i, j int) (val float64) {
		_, val = GradSimplex2DPTerm(r, s, i, j)
		//fmt.Printf("Ds r,s,i,j,val = %8.5f,%8.5f,%d,%d,%8.5f,", r, s, i, j, val)
		return
	}
	switch deriv {
	case None:
		tFunc = ONTerm2D
	case Dr:
		tFunc = DrONTerm2D
	case Ds:
		tFunc = DsONTerm2D
	}
	p1, p2 = make([]float64, Np), make([]float64, Np)
	// Evaluate the full 2D polynomial basis first, once for each of two components
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			val := tFunc(r, s, i, j)
			p1[sk] = val
			p2[sk+N2DBasis] = val
			sk++
		}
	}
	// Evaluate the term ([ X ]*(Pk)) at only the top N+1 terms (highest order) of the 2D polynomial
	sk += N2DBasis // Skip to the beginning of the second polynomial group
	for i := 0; i <= N; i++ {
		j := N - i
		val := tFunc(r, s, i, j)
		switch deriv {
		case None:
			p1[sk] = val * r
			p2[sk] = val * s
		case Dr:
			val2 := Simplex2DPTerm(r, s, i, j)
			p1[sk] = val2 + val*r
			p2[sk] = val * s
		case Ds:
			val2 := Simplex2DPTerm(r, s, i, j)
			p1[sk] = val * r
			p2[sk] = val2 + val*s
		}
		sk++
	}
	return
}

func (rt *RTElement) EvaluateLagrangeRTBasis(l2dp *LagrangeBasis2D, r, s float64, derivO ...DerivativeDirection) (p1, p2 []float64) {
	var (
		sk       int
		N        = rt.N
		Np       = (N + 1) * (N + 3)
		N2DBasis = (N + 1) * (N + 2) / 2 // Number of polynomial terms for each of R and S directions
		deriv    = None
		tFunc    func(r, s float64, i, j int) (val float64)
	)
	if len(derivO) != 0 {
		deriv = derivO[0]
	}
	Term2D := func(r, s float64, i, j int) (val float64) {
		var (
			R, S = utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s})
		)
		return l2dp.BasisPolynomialTerm(R, S, i, j)[0]
	}
	DrTerm2D := func(r, s float64, i, j int) (val float64) {
		var (
			R, S = utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s})
		)
		valSlice, _ := l2dp.GradBasisPolynomialTerm(R, S, i, j)
		val = valSlice[0]
		//fmt.Printf("Dr r,s,i,j,val = %8.5f,%8.5f,%d,%d,%8.5f,", r, s, i, j, val)
		return
	}
	DsTerm2D := func(r, s float64, i, j int) (val float64) {
		var (
			R, S = utils.NewVector(1, []float64{r}), utils.NewVector(1, []float64{s})
		)
		_, valSlice := l2dp.GradBasisPolynomialTerm(R, S, i, j)
		val = valSlice[0]
		//fmt.Printf("Ds r,s,i,j,val = %8.5f,%8.5f,%d,%d,%8.5f,", r, s, i, j, val)
		return
	}
	switch deriv {
	case None:
		tFunc = Term2D
	case Dr:
		tFunc = DrTerm2D
	case Ds:
		tFunc = DsTerm2D
	}
	p1, p2 = make([]float64, Np), make([]float64, Np)
	// Evaluate the full 2D polynomial basis first, once for each of two components
	//fmt.Printf("N = %d\n", N)
	for i := 0; i <= N; i++ {
		for j := 0; j <= (N - i); j++ {
			val := tFunc(r, s, i, j)
			p1[sk] = val
			p2[sk+N2DBasis] = val
			sk++
		}
	}
	// Evaluate the term ([ X ]*(Pk)) at only the top N+1 terms (highest order) of the 2D polynomial
	sk += N2DBasis // Skip to the beginning of the second polynomial group
	for i := 0; i <= N; i++ {
		j := N - i
		val := tFunc(r, s, i, j)
		switch deriv {
		case None:
			p1[sk] = val * r
			p2[sk] = val * s
		case Dr:
			val2 := Term2D(r, s, i, j)
			p1[sk] = val2 + val*r
			p2[sk] = val * s
		case Ds:
			val2 := Term2D(r, s, i, j)
			p1[sk] = val * r
			p2[sk] = val2 + val*s
		}
		sk++
	}
	return
}

func ExtendGeomToRT(N int, rInt, sInt utils.Vector) (r, s utils.Vector) {
	var (
		NpEdge       = N + 1
		rData, sData = rInt.DataP, sInt.DataP
	)
	/*
		Determine geometric locations of edge points, located at Gauss locations in 1D, projected onto the edges
	*/
	GQR, _ := DG1D.JacobiGQ(1, 1, N)
	/*
		// Equi-spaced edge
		gqr := make([]float64, NpEdge)
		space := 2. / float64(NpEdge+1)
		start := -1.
		for i := 0; i < NpEdge; i++ {
			gqr[i] = start + space
			start = gqr[i]
		}
		GQR := utils.NewVector(NpEdge, gqr)
	*/
	/*
		Double the number of interior points to match each direction of the basis
	*/
	if N == 0 { // Special case: when N=0, the interior of the RT element is empty
		rData, sData = []float64{}, []float64{}
	} else {
		rData = append(rData, rData...)
		sData = append(sData, sData...)
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

func (rt *RTElement) GetInternalLocations(F utils.Vector) (Finternal []float64) {
	var (
		Nint = rt.Nint
	)
	Finternal = make([]float64, Nint)
	for i := 0; i < Nint; i++ {
		Finternal[i] = F.DataP[i]
	}
	return
}

func (rt *RTElement) GetEdgeLocations(F utils.Vector) (Fedge []float64) {
	var (
		Nint     = rt.Nint
		NedgeTot = rt.Nedge * 3
	)
	Fedge = make([]float64, NedgeTot)
	for i := 0; i < NedgeTot; i++ {
		Fedge[i] = F.DataP[i+2*Nint]
	}
	return
}
