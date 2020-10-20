package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"
)

type RTElement struct {
	N              int             // Order of element
	Np             int             // Number of points in element
	Nedge, Nint    int             // Number of Edge and Interior points
	Ainv           utils.Matrix    // Polynomial coefficient matrix
	V              [2]utils.Matrix // Vandermonde matrix for each direction r and s
	Vr, Vs         [2]utils.Matrix // Derivative matrices in r and s directions
	Vr0Int, Vs1Int utils.Matrix    // Divergence of basis, restricted to internal points only
	R, S           utils.Vector    // Point locations defining element in [-1,1] Triangle
}

type RTPointType uint

const (
	All RTPointType = iota
	InteriorR
	InteriorS
	Edge1
	Edge2
	Edge3
)

func NewRTElement(N int, R, S utils.Vector) (rt *RTElement) {
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
	rt.CalculateBasis()
	return
}

func (rt *RTElement) ProjectFunctionOntoBasis2(s1, s2 []float64) (s1p, s2p []float64) {
	var (
		Np    = len(s1)
		Nint  = rt.Nint
		Nedge = rt.Nedge
	)
	s1p, s2p = make([]float64, Np), make([]float64, Np)

	oosr2 := 1 / math.Sqrt(2)
	for i := 0; i < Nint; i++ {
		s1p[i] = s1[i]
		s2p[i] = s2[i]
	}
	for i := 0; i < Nedge; i++ {
		// Edge1: Unit vector is [1/sqrt(2), 1/sqrt(2)]
		dp := oosr2 * (s1[i+Nint] + s2[i+Nint])
		s1p[i+Nint] = oosr2 * dp
		s2p[i+Nint] = oosr2 * dp
		// Edge2: Unit vector is [-1,0]
		s1p[i+Nint+Nedge] = -s1[i+Nint+Nedge]
		s2p[i+Nint+Nedge] = 0
		// Edge3: // Unit vector is [0,-1]
		s1p[i+Nint+2*Nedge] = 0
		s2p[i+Nint+2*Nedge] = -s2[i+Nint+2*Nedge]
	}
	return
}

func (rt *RTElement) ProjectFunctionOntoBasis(s1, s2 []float64) (s1p, s2p []float64) {
	var (
		Np = len(s1)
	)
	s1p, s2p = make([]float64, Np), make([]float64, Np)

	oosr2 := 1 / math.Sqrt(2)
	for i := range s1 {
		switch rt.GetTermType(i) {
		case InteriorR:
			// Unit vector is [1,0]
			s1p[i] = s1[i]
		case InteriorS:
			// Unit vector is [0,1]
			s2p[i] = s2[i]
		case Edge1:
			// Edge1: Unit vector is [1/sqrt(2), 1/sqrt(2)]
			s1p[i] = s1[i] * oosr2
			s2p[i] = s2[i] * oosr2
		case Edge2:
			// Edge2: Unit vector is [-1,0]
			s1p[i] = -s1[i]
		case Edge3:
			// Edge3: // Unit vector is [0,-1]
			s2p[i] = -s2[i]
		}
	}
	return
}

func (rt *RTElement) ExtractFunctionFromBasis(s1p, s2p []float64) (s1, s2 []float64) {
	var (
		Np   = len(s1p)
		Nint = rt.Nint
	)
	s1, s2 = make([]float64, Np), make([]float64, Np)

	for i := range s1 {
		switch rt.GetTermType(i) {
		case InteriorR:
			// Unit vector is [1,0]
			s1[i] = s1p[i]
			s2[i] = s2p[i+Nint]
		case InteriorS:
			// Unit vector is [0,1]
			s1[i] = s1p[i-Nint]
			s2[i] = s2p[i]
		case Edge1, Edge2, Edge3:
			s1[i] = s1p[i]
			s2[i] = s2p[i]
		}
	}
	return
}

func (rt *RTElement) GetTermType(i int) (rtt RTPointType) {
	var (
		N         = rt.N
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

func (rt *RTElement) CalculateBasis() {
	/*
				This is constructed from the defining space of the RT element:
								 2
					RT_k = [(P_k) ]   + [ X ] P_k
						 = [ b1(r,s)_i + r * b3(r,s)_j ]
						   [ b2(r,s)_i + s * b3(r,s)_j ]
					j := 1, (K+1)(K+2)/2
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
		N      = rt.N
		R, S   = rt.R, rt.S
		Np     = (N + 1) * (N + 3)
		A      utils.Matrix
		p1, p2 []float64
	)
	// Add the edge and additional interior (duplicated) points to complete the RT geometry2D
	rt.R, rt.S = ExtendGeomToRT(N, R, S)
	/*
		Form the basis matrix by forming a dot product with unit vectors, matching the coordinate locations in R,S
	*/
	A = utils.NewMatrix(Np, Np)
	rowEdge := make([]float64, Np)
	oosr2 := 1 / math.Sqrt(2)
	_ = oosr2

	// Evaluate at geometric locations
	for ii, rr := range rt.R.Data() {
		ss := rt.S.Data()[ii]
		/*
			First, evaluate the polynomial at the (r,s) coordinates
			This is the same set that will be used for all dot products to form the basis matrix
		*/
		p1, p2 = rt.EvaluateRTBasis(rr, ss)
		// Implement dot product of (unit vector)_ii with each vector term in the polynomial evaluated at location ii
		switch rt.GetTermType(ii) {
		case InteriorR:
			// Unit vector is [1,0]
			A.M.SetRow(ii, p1)
		case InteriorS:
			// Unit vector is [0,1]
			A.M.SetRow(ii, p2)
		case Edge1:
			for i := range rowEdge {
				// Edge1: Unit vector is [1/sqrt(2), 1/sqrt(2)]
				rowEdge[i] = oosr2 * (p1[i] + p2[i])
			}
			A.M.SetRow(ii, rowEdge)
		case Edge2:
			for i := range rowEdge {
				// Edge2: Unit vector is [-1,0]
				rowEdge[i] = -p1[i]
			}
			A.M.SetRow(ii, rowEdge)
		case Edge3:
			for i := range rowEdge {
				// Edge3: // Unit vector is [0,-1]
				rowEdge[i] = -p2[i]
			}
			A.M.SetRow(ii, rowEdge)
		}
	}
	// Invert [ Ainv ] to obtain the coefficients (columns) of polynomials (rows), each row is a polynomial
	var Ainv utils.Matrix
	if Ainv, err = A.Inverse(); err != nil {
		panic(err)
	}
	/*
		Process the coefficient matrix to produce each direction of each polynomial (V1 and V2)
	*/
	rt.V[0], rt.V[1] = utils.NewMatrix(Np, Np), utils.NewMatrix(Np, Np)
	rt.Ainv = utils.NewMatrix(Np, Np)
	for j := 0; j < Np; j++ { // Process a column at a time, each column is a polynomial
		/*
			First, we extract the coefficients for each direction in this j-th polynomial
		*/
		coeffs := Ainv.Col(j).Data() // All coefficients of the j-th polynomial terms
		p1, p2 := make([]float64, Np), make([]float64, Np)
		p := make([]float64, Np)
		NGroup1 := (N+1)*N/2 + (N + 1) // The first set of polynomial terms
		for i, coeff := range coeffs {
			// one full set of polynomial coeffs twice, then a set for the last N+1 terms
			switch {
			case i < NGroup1:
				p1[i] = coeff
				p2[i] = 0
			case i >= NGroup1 && i < 2*NGroup1:
				p1[i] = 0
				p2[i] = coeff
			default:
				p1[i] = coeff
				p2[i] = coeff
			}
			p[i] = coeff
		}
		/*
			Load the evaluated polynomial into the j-th column of the Coefficient matrices
		*/
		rt.Ainv.SetCol(j, p)
		/*
			Evaluate the basis at this location, then multiply each term by its polynomial coefficient
		*/
		px, py := make([]float64, rt.R.Len()), make([]float64, rt.R.Len())
		for i, rr := range rt.R.Data() {
			ss := rt.S.Data()[i]
			b1, b2 := rt.EvaluateRTBasis(rr, ss) // All terms of the basis at this location
			for ii := range p1 {
				// Multiply each coefficient by its basis term to create the evaluated polynomial terms column
				px[i] += p1[ii] * b1[ii]
				py[i] += p2[ii] * b2[ii]
			}
		}
		/*
			Load the evaluated polynomial into the j-th column of the Vandermonde matrices
		*/
		rt.V[0].SetCol(j, px)
		rt.V[1].SetCol(j, py)
	}
	// Create derivative matrices, Vr and Vs
	rt.Vr[0], rt.Vr[1] = utils.NewMatrix(Np, Np), utils.NewMatrix(Np, Np)
	rt.Vs[0], rt.Vs[1] = utils.NewMatrix(Np, Np), utils.NewMatrix(Np, Np)
	for i := 0; i < Np; i++ { // Each geometric location
		rr, ss := rt.R.Data()[i], rt.S.Data()[i]
		for j := 0; j < Np; j++ { // Each polynomial
			p1r, p2r := rt.EvaluatePolynomial(j, rr, ss, Dr)
			rt.Vr[0].Set(i, j, p1r)
			rt.Vr[1].Set(i, j, p2r)
			p1s, p2s := rt.EvaluatePolynomial(j, rr, ss, Ds)
			rt.Vs[0].Set(i, j, p1s)
			rt.Vs[1].Set(i, j, p2s)
		}
	}
	rt.Np = Np
	return
}

func (rt *RTElement) Divergence(f1, f2 []float64) (div []float64) {
	if len(f1) != rt.Np || len(f2) != rt.Np {
		panic(fmt.Errorf("wrong input number of points, should be %d, is %d\n", rt.Np, len(f1)))
	}
	f1p, f2p := rt.ProjectFunctionOntoBasis(f1, f2)
	var (
		Npm = rt.Np
	)
	f1pV, f2pV := utils.NewMatrix(Npm, 1, f1p), utils.NewMatrix(Npm, 1, f2p)
	// Divergence is (Dr1 * f1pV) + (Ds2 * f2pV)
	divV := rt.Vr[0].Mul(f1pV).Add(rt.Vs[1].Mul(f2pV))
	div = divV.Data()
	return
}

func (rt *RTElement) DivergenceInterior(f1, f2 []float64) (div []float64) {
	if len(f1) != rt.Np || len(f2) != rt.Np {
		panic(fmt.Errorf("wrong input number of points, should be %d, is %d\n", rt.Np, len(f1)))
	}
	f1p, f2p := rt.ProjectFunctionOntoBasis(f1, f2)
	var (
		Npm       = rt.Np
		N         = rt.N
		NInterior = N * (N + 1) / 2 // one order less than RT element in (P_k)2
	)
	if rt.Vr0Int.IsEmpty() {
		// Restrict the derivative matrices to the interior points
		rt.Vr0Int, rt.Vs1Int = utils.NewMatrix(NInterior, Npm), utils.NewMatrix(NInterior, Npm)
		for i := 0; i < NInterior; i++ {
			rowDr1, rowDs2 := rt.Vr[0].Row(i).Data(), rt.Vs[1].Row(i).Data()
			rt.Vr0Int.M.SetRow(i, rowDr1)
			rt.Vs1Int.M.SetRow(i, rowDs2)
		}
	}
	f1pV, f2pV := utils.NewMatrix(Npm, 1, f1p), utils.NewMatrix(Npm, 1, f2p)
	// Divergence is (Dr1 * f1pV) + (Ds2 * f2pV)
	divV := rt.Vr0Int.Mul(f1pV).Add(rt.Vs1Int.Mul(f2pV))
	div = divV.Data()
	return
}

type DerivativeDirection uint8

const (
	None DerivativeDirection = iota
	Dr
	Ds
)

func (rt *RTElement) EvaluateRTBasis(r, s float64, derivO ...DerivativeDirection) (p1, p2 []float64) {
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
	DrONTerm2D := func(r, s float64, i, j int) (val float64) {
		val, _ = GradSimplex2DPTerm(r, s, i, j)
		//fmt.Printf("Vr r,s,i,j,val = %8.5f,%8.5f,%d,%d,%8.5f,", r, s, i, j, val)
		return
	}
	DsONTerm2D := func(r, s float64, i, j int) (val float64) {
		_, val = GradSimplex2DPTerm(r, s, i, j)
		//fmt.Printf("Vs r,s,i,j,val = %8.5f,%8.5f,%d,%d,%8.5f,", r, s, i, j, val)
		return
	}
	switch deriv {
	case None:
		tFunc = Simplex2DPTerm
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

func (rt *RTElement) EvaluatePolynomial(j int, r, s float64, derivO ...DerivativeDirection) (p1, p2 float64) {
	/*
		Get the coefficients for the j-th polynomial and compute:
			p(r,s) = sum(coeff_i*P_i(r,s))
		for each direction [1,2]
	*/
	var (
		deriv = None
	)
	if len(derivO) != 0 {
		deriv = derivO[0]
	}
	coeffs := rt.Ainv.Col(j).Data()
	b1, b2 := rt.EvaluateRTBasis(r, s, deriv)
	for i := range coeffs {
		//fmt.Printf("term(%d} = [%8.5f, %8.5f]\n", i, coeffs1[i]*b1[i], coeffs2[i]*b2[i])
		p1 += coeffs[i] * b1[i]
		p2 += coeffs[i] * b2[i]
	}
	return
}

func ExtendGeomToRT(N int, rInt, sInt utils.Vector) (r, s utils.Vector) {
	var (
		NpEdge       = N + 1
		rData, sData = rInt.Data(), sInt.Data()
	)
	/*
		Determine geometric locations of edge points, located at Gauss locations in 1D, projected onto the edges
	*/
	GQR, _ := DG1D.JacobiGQ(1, 1, N)
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
	GQRData := GQR.Data()
	rEdgeData := make([]float64, NpEdge*3)
	sEdgeData := make([]float64, NpEdge*3)
	for i := 0; i < NpEdge; i++ {
		gp := GQRData[i]
		// Edge 1 (hypotenuse)
		gpT := 0.5 * (gp + 1)
		rEdgeData[i] = 1 - 2*gpT
		sEdgeData[i] = -1 + 2*gpT
		// Edge 2
		rEdgeData[i+NpEdge] = -1
		sEdgeData[i+NpEdge] = gp
		// Edge 3
		rEdgeData[i+2*NpEdge] = gp
		sEdgeData[i+2*NpEdge] = -1
	}
	rData = append(rData, rEdgeData...)
	sData = append(sData, sEdgeData...)
	r = utils.NewVector(len(rData), rData)
	s = utils.NewVector(len(sData), sData)
	return
}

func NodesEpsilon(N int) (R, S utils.Vector) {
	/*
		From the 2017 paper "Ainv Direct Flux Reconstruction Scheme for Advection Diffusion Problems on Triangular Grids"

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
	// Cases 3,4 from Romero and Jameson, Others from Williams and Shun
	case 0:
		epsD = []float64{
			0.3333333333333333,
			0.3333333333333333,
			0.3333333333333333,
		}
	case 1:
		epsD = []float64{
			0.666666666666667, 0.166666666666667, 0.166666666666667,
			0.166666666666667, 0.666666666666667, 0.166666666666667,
			0.166666666666667, 0.166666666666667, 0.666666666666667,
		}
	case 2:
		epsD = []float64{
			0.816847572980440, 0.091576213509780, 0.091576213509780, 0.445948490915964, 0.445948490915964, 0.108103018168071,
			0.091576213509780, 0.816847572980440, 0.091576213509780, 0.445948490915964, 0.108103018168071, 0.445948490915964,
			0.091576213509780, 0.091576213509780, 0.816847572980440, 0.108103018168071, 0.445948490915964, 0.445948490915964,
		}
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
	case 5:
		epsD = []float64{
			0.943774095634672, 0.028112952182664, 0.028112952182664, 0.645721803061365, 0.177139098469317, 0.177139098469317, 0.405508595867433, 0.405508595867433, 0.188982808265134, 0.148565812270887, 0.148565812270887, 0.033533207700614, 0.817900980028499, 0.817900980028499, 0.033533207700614, 0.357196298615681, 0.357196298615681, 0.037824789609186, 0.604978911775132, 0.604978911775132, 0.037824789609186,
			0.028112952182664, 0.943774095634672, 0.028112952182664, 0.177139098469317, 0.645721803061365, 0.177139098469317, 0.405508595867433, 0.188982808265134, 0.405508595867433, 0.817900980028499, 0.033533207700614, 0.148565812270887, 0.148565812270887, 0.033533207700614, 0.817900980028499, 0.604978911775132, 0.037824789609186, 0.357196298615681, 0.357196298615681, 0.037824789609186, 0.604978911775132,
			0.028112952182664, 0.028112952182664, 0.943774095634672, 0.177139098469317, 0.177139098469317, 0.645721803061365, 0.188982808265134, 0.405508595867433, 0.405508595867433, 0.033533207700614, 0.817900980028499, 0.817900980028499, 0.033533207700614, 0.148565812270887, 0.148565812270887, 0.037824789609186, 0.604978911775132, 0.604978911775132, 0.037824789609186, 0.357196298615681, 0.357196298615681,
		}
	case 6:
		epsD = []float64{
			0.960045625755613, 0.019977187122193, 0.019977187122193, 0.736556464940005, 0.131721767529998, 0.131721767529998, 0.333333333333333, 0.485135346793461, 0.485135346793461, 0.029729306413079, 0.107951981846011, 0.107951981846011, 0.024136808036039, 0.867911210117951, 0.867911210117951, 0.024136808036039, 0.270840772921567, 0.270840772921567, 0.028286656697710, 0.700872570380723, 0.700872570380723, 0.028286656697710, 0.316549598844617, 0.316549598844617, 0.146795716949245, 0.536654684206138, 0.536654684206138, 0.146795716949245,
			0.019977187122193, 0.960045625755613, 0.019977187122193, 0.131721767529998, 0.736556464940005, 0.131721767529998, 0.333333333333333, 0.485135346793461, 0.029729306413079, 0.485135346793461, 0.867911210117951, 0.024136808036039, 0.107951981846011, 0.107951981846011, 0.024136808036039, 0.867911210117951, 0.700872570380723, 0.028286656697710, 0.270840772921567, 0.270840772921567, 0.028286656697710, 0.700872570380723, 0.536654684206138, 0.146795716949245, 0.316549598844617, 0.316549598844617, 0.146795716949245, 0.536654684206138,
			0.019977187122193, 0.019977187122193, 0.960045625755613, 0.131721767529998, 0.131721767529998, 0.736556464940005, 0.333333333333333, 0.029729306413079, 0.485135346793461, 0.485135346793461, 0.024136808036039, 0.867911210117951, 0.867911210117951, 0.024136808036039, 0.107951981846011, 0.107951981846011, 0.028286656697710, 0.700872570380723, 0.700872570380723, 0.028286656697710, 0.270840772921567, 0.270840772921567, 0.146795716949245, 0.536654684206138, 0.536654684206138, 0.146795716949245, 0.316549598844617, 0.316549598844617,
		}
	case 7:
		epsD = []float64{
			0.957657154441070, 0.021171422779465, 0.021171422779465, 0.798831205208225, 0.100584397395888, 0.100584397395888, 0.457923384576135, 0.271038307711932, 0.271038307711932, 0.440191258403832, 0.440191258403832, 0.119617483192335, 0.101763679498021, 0.101763679498021, 0.018256679074748, 0.879979641427232, 0.879979641427232, 0.018256679074748, 0.394033271669987, 0.394033271669987, 0.023404705466341, 0.582562022863673, 0.582562022863673, 0.023404705466341, 0.226245530909229, 0.226245530909229, 0.022223854547989, 0.751530614542782, 0.751530614542782, 0.022223854547989, 0.635737183263105, 0.635737183263105, 0.115183589115563, 0.249079227621332, 0.249079227621332, 0.115183589115563,
			0.021171422779465, 0.957657154441070, 0.021171422779465, 0.100584397395888, 0.798831205208225, 0.100584397395888, 0.271038307711932, 0.457923384576135, 0.271038307711932, 0.440191258403832, 0.119617483192335, 0.440191258403832, 0.879979641427232, 0.018256679074748, 0.101763679498021, 0.101763679498021, 0.018256679074748, 0.879979641427232, 0.582562022863673, 0.023404705466341, 0.394033271669987, 0.394033271669987, 0.023404705466341, 0.582562022863673, 0.751530614542782, 0.022223854547989, 0.226245530909229, 0.226245530909229, 0.022223854547989, 0.751530614542782, 0.249079227621332, 0.115183589115563, 0.635737183263105, 0.635737183263105, 0.115183589115563, 0.249079227621332,
			0.021171422779465, 0.021171422779465, 0.957657154441070, 0.100584397395888, 0.100584397395888, 0.798831205208225, 0.271038307711932, 0.271038307711932, 0.457923384576135, 0.119617483192335, 0.440191258403832, 0.440191258403832, 0.018256679074748, 0.879979641427232, 0.879979641427232, 0.018256679074748, 0.101763679498021, 0.101763679498021, 0.023404705466341, 0.582562022863673, 0.582562022863673, 0.023404705466341, 0.394033271669987, 0.394033271669987, 0.022223854547989, 0.751530614542782, 0.751530614542782, 0.022223854547989, 0.226245530909229, 0.226245530909229, 0.115183589115563, 0.249079227621332, 0.249079227621332, 0.115183589115563, 0.635737183263105, 0.635737183263105,
		}
	default:
		panic(fmt.Errorf("Epsilon nodes not defined for N = %v\n", N))
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

func CombineBasis(N int, V1, V2 utils.Matrix) (V1m, V2m utils.Matrix) {
	// Merge the R and S terms for the interior points, where there is duplication in the basis
	var (
		Npm  = (N + 6) * (N + 1) / 2 // Number of points in basis with one set of interior points with a vector for R and S
		Np   = (N + 3) * (N + 1)
		Nint = N * (N + 1) / 2
	)
	trimFloat := func(in float64) (out float64) {
		thresh := 1.e-12
		if math.Abs(in) > thresh {
			out = in
		}
		return
	}
	V1m, V2m = utils.NewMatrix(Npm, Npm), utils.NewMatrix(Npm, Npm)
	var j1m, j2m int
	for j := 0; j < Np; j++ {
		var i1m, i2m int
		for i := 0; i < Np; i++ {
			V1m.M.Set(i1m, j1m, trimFloat(V1.At(i, j)))
			V2m.M.Set(i2m, j2m, trimFloat(V2.At(i, j)))
			if i < Nint || i >= 2*Nint {
				i1m++
			}
			if i >= Nint {
				i2m++
			}
		}
		if j < Nint || j >= 2*Nint {
			j1m++
		}
		if j >= Nint {
			j2m++
		}
	}
	return
}

func (rt *RTElement) GetInternalLocations(F utils.Vector) (Finternal []float64) {
	var (
		Nint = rt.Nint
	)
	Finternal = make([]float64, Nint)
	for i := 0; i < Nint; i++ {
		Finternal[i] = F.Data()[i]
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
		Fedge[i] = F.Data()[i+2*Nint]
	}
	return
}
