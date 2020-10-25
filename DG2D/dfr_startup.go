package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
)

type DFR2D struct {
	N, Np            int
	SolutionElement  *LagrangeElement2D
	FluxElement      *RTElement
	FluxInterpMatrix utils.Matrix // A pre-calculated interpolation matrix covering all Flux (edge) points in K elements
	// Mesh Parameters
	K                    int             // Number of elements (triangles) in mesh
	VX, VY               utils.Vector    // X,Y vertex points in mesh (vertices)
	BCType               utils.Matrix    // Mapping of elements to vertices, each element has three integer vertex coordinates
	FluxX, FluxY         utils.Matrix    // Flux Element local coordinates
	SolutionX, SolutionY utils.Matrix    // Solution Element local coordinates
	Tris                 *Triangulation  // Triangle mesh and edge/face structures
	J, Jinv, Jdet        utils.Matrix    // Mesh Transform Jacobian, Kx[4], each K element has a 2x2 matrix, det is |J|
	FaceNorm             [2]utils.Matrix // Magnitude of each face normal, 3xK, used in projection of flux to RT element
}

func NewDFR2D(N int, meshFileO ...string) (dfr *DFR2D) {
	if N < 1 {
		panic(fmt.Errorf("Polynomial order must be >= 1, have %d", N))
	}
	le := NewLagrangeElement2D(N, Epsilon)
	rt := NewRTElement(N+1, le.R, le.S)
	RFlux := utils.NewVector(rt.Nedge*3, rt.GetEdgeLocations(rt.R)) // For the Interpolation matrix across three edges
	SFlux := utils.NewVector(rt.Nedge*3, rt.GetEdgeLocations(rt.S)) // For the Interpolation matrix across three edges
	dfr = &DFR2D{
		SolutionElement:  le,
		FluxElement:      rt,
		FluxInterpMatrix: le.Simplex2DInterpolatingPolyMatrix(RFlux, SFlux), // Interpolation matrix across three edges
	}
	if len(meshFileO) != 0 {
		var EToV utils.Matrix
		dfr.K, dfr.VX, dfr.VY, EToV, dfr.BCType = ReadGambit2d(meshFileO[0])
		dfr.Tris = NewTriangulation(EToV, dfr.BCType)
		// Build connectivity matrices
		dfr.FluxX, dfr.FluxY =
			CalculateElementLocalGeometry(dfr.Tris.EToV, dfr.VX, dfr.VY, dfr.FluxElement.R, dfr.FluxElement.S)
		dfr.SolutionX, dfr.SolutionY =
			CalculateElementLocalGeometry(dfr.Tris.EToV, dfr.VX, dfr.VY, dfr.SolutionElement.R, dfr.SolutionElement.S)
		dfr.FluxX.SetReadOnly("FluxX")
		dfr.FluxY.SetReadOnly("FluxY")
		dfr.SolutionX.SetReadOnly("SolutionX")
		dfr.SolutionY.SetReadOnly("SolutionY")
		dfr.CalculateJacobian()
		dfr.CalculateFaceNorms()
	}
	return
}

func (dfr *DFR2D) CalculateJacobian() {
	Jd := make([]float64, 4*dfr.K)
	Jdetd := make([]float64, dfr.K)
	JdInv := make([]float64, 4*dfr.K)
	for k := 0; k < dfr.K; k++ {
		tri := dfr.Tris.EToV.Row(k).Data()
		v := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
		v1x, v2x, v3x := dfr.VX.AtVec(v[0]), dfr.VX.AtVec(v[1]), dfr.VX.AtVec(v[2])
		v1y, v2y, v3y := dfr.VY.AtVec(v[0]), dfr.VY.AtVec(v[1]), dfr.VY.AtVec(v[2])
		xr, yr := 0.5*(v2x-v1x), 0.5*(v2y-v1y)
		xs, ys := 0.5*(v3x-v1x), 0.5*(v3y-v1y)
		// Jacobian is [xr, yr]
		//             [xs, ys]
		jd := Jd[k*4:]
		//jd[0], jd[1], jd[2], jd[3] = xr, xs, yr, ys
		jd[0], jd[1], jd[2], jd[3] = xr, yr, xs, ys
		Jdetd[k] = xr*ys - xs*yr
		oodet := 1. / Jdetd[k]
		jdInv := JdInv[k*4:]
		//jdInv[0], jdInv[1], jdInv[2], jdInv[3] = oodet*ys, -oodet*xs, -oodet*yr, oodet*xr
		jdInv[0], jdInv[1], jdInv[2], jdInv[3] = oodet*ys, -oodet*yr, -oodet*xs, oodet*xr
	}
	dfr.J, dfr.Jinv = utils.NewMatrix(dfr.K, 4, Jd), utils.NewMatrix(dfr.K, 4, JdInv)
	dfr.Jdet = utils.NewMatrix(dfr.K, 1, Jdetd)
	dfr.J.SetReadOnly("J")
	dfr.Jdet.SetReadOnly("Jdet")
	dfr.Jinv.SetReadOnly("Jinv")
}

func (dfr *DFR2D) CalculateFaceNorms() {
	dfr.FaceNorm[0], dfr.FaceNorm[1] = utils.NewMatrix(dfr.K, 3), utils.NewMatrix(dfr.K, 3)
	for en, e := range dfr.Tris.Edges {
		for ii, triNum := range e.ConnectedTris {
			//fn := dfr.FaceNorm.Row(int(triNum)).Data()[0:3]
			k := int(triNum)
			fnD1, fnD2 := dfr.FaceNorm[0].Data(), dfr.FaceNorm[1].Data()
			x1, x2 := dfr.Tris.GetEdgeCoordinates(en, e, triNum, dfr.VX, dfr.VY)
			dx, dy := x2[0]-x1[0], x2[1]-x1[1]
			nx, ny := -dy, dx
			switch e.ConnectedTriEdgeNumber[ii] {
			case First:
				fnD1[0+k*3], fnD2[0+k*3] = nx, ny
			case Second:
				fnD1[1+k*3], fnD2[1+k*3] = nx, ny
			case Third:
				fnD1[2+k*3], fnD2[2+k*3] = nx, ny
			}
		}
	}
	dfr.FaceNorm[0].SetReadOnly("FaceNorm1")
	dfr.FaceNorm[1].SetReadOnly("FaceNorm2")
}

// Build reference element matrices
/*
		We build the mixed elements for the DFR scheme with:

		Solution Points: We use points within a reference triangle, excluding the edges, for a Lagrangian element
		of O(K) to store the solution. If we need derivatives, or interpolated quantities (Flux), we use the
		solution points element.

		Flux Points: We use a customized Raviart-Thomas (RT) vector element of O(K+1) to store the vector Flux function
	    computed from the solution values. The RT element is of order O(K+1) and is a combination of the points from
		the solution element for the interior, and points along the three triangle edges. The custom RT basis is
		established using a procedure outlined in: "Ainv Direct Flux Reconstruction Scheme for Advection-Diffusion
		Problems on Triangular Grids" by Romero, Witherden and Jameson (2017). Ainv complete RT basis, [ B ], is used
		together with unit basis vectors, [ w ], to satisfy the following:
				[ B_j(r_i) dot w_i ] [ C ] = [ delta_i_j ]
				=> solve for [ C ], the coefficients defining the custom RT basis

		[ C ] is the vector of coefficients defining the basis using the basis vectors [ w ] and [ B ].

		The [ w ] directions of the custom RT element basis are defined such that:
			w([r]) = w(edge_locations) = unit normals on each of three edges
			w([r]) = w(interior) = unit normals in the two primary geometry2D directions (r and s)

		For order K there are:
			- (K+1) locations on each edge, for a total of 3(K+1) edge basis functions.
			- (K)(K+1) locations in the interior, half for the w_r direction and half for the w_s direction
			- Total: (K+3)(K+1) basis functions for the custom RT_K element

		Notes:
			1) The number of interior points matches the Lagrangian element in 2D at order (K-1). Ainv Lagrange element
			at order (K) has N_p = (K+1)(K+2)/2 degrees of freedom, so an order (K-1) element has (K)(K+1)/2 DOF.
			Considering that we need a term for each of the two interior directions at each interior point, we need
			exactly 2*N_p DOF at order (K-1) for the interior of the custom RT element, resulting in (K)(K+1) terms.
			2) Note (1) confirms that the custom element requires exactly the same number of interior points
			(K)(K+1)/2 as a Lagrange element of order (K-1), which means we can use the custom RT element for the
			DFR approach, which needs to provide a O(K+1) element to preserve the gradient at O(K). We will use the
			solution points from the Lagrange element at O(K) to construct the interior of the O(K+1) RT element
			without requiring interpolation of the solution points, as they already reside at the same geometric
			locations.
			(3) To create the custom RT element, we initialize the solution element, then define the custom RT element
			from the interior point locations of the solution element to ensure that they are colocated.
			(4) To use the custom RT element:
			a) calculate the solution, calculate the flux vector field from the solution at the solution points
			b) transfer the flux vector field values to the DFR element interior
			c) interpolate flux values at from the interior of the RT element to the locations on the triangle edges
			d) use the method of characteristics to calculate the corrected flux using the neighbor element's edge
			flux combined with the edge flux from this element
			e) calculate the gradient of the vector flux field using the custom RT element
			f) transfer the gradient values from the RT element to the solution element for use in advancing the
			solution in differential form (directly)

		By calculating the flux gradient in a way that yields an O(K) polynomial on the solution points, we can use
		the differential form of the equations directly for the solution, rather than using the traditional Galerkin
		approach of repeated integration by parts to obtain an equation with only first derivatives. This simplifies
		the solution process, resulting in a more efficient computational approach, in addition to making it easier
		to solve more complex equations with the identical formulation.
*/

func (dfr *DFR2D) PiolaTransform(J []float64, Jdet float64, f [2]float64) (fT [2]float64) {
	ooJdet := 1. / Jdet
	//fT[0], fT[1] = ooJdet*(J[0]*f[0]+J[2]*f[1]), ooJdet*(J[1]*f[0]+J[3]*f[1])
	fT[0], fT[1] = ooJdet*(J[0]*f[0]+J[1]*f[1]), ooJdet*(J[2]*f[0]+J[3]*f[1])
	return
}

func (dfr *DFR2D) ProjectFluxOntoRTSpace(Fx, Fy utils.Matrix) (Fxp, Fyp utils.Matrix) {
	var (
		Np    = dfr.FluxElement.Np
		K     = dfr.K
		rt    = dfr.FluxElement
		oosr2 = 1 / math.Sqrt(2)
	)
	Fxp, Fyp = utils.NewMatrix(K, Np), utils.NewMatrix(K, Np)
	for k := 0; k < K; k++ {
		var (
			edgeNorm1, edgeNorm2 = dfr.FaceNorm[0].Row(k).Data()[0:3], dfr.FaceNorm[1].Row(k).Data()[0:3]
			Jdet                 = dfr.Jdet.Row(k).Data()[0]
			Jinv                 = dfr.Jinv.Row(k).Data()[0:4]
			fxD, fyD             = Fx.Row(k).Data()[0:Np], Fy.Row(k).Data()[0:Np]
			fxpD, fypD           = Fxp.Data(), Fyp.Data()
		)
		getEdgeQuants := func(fxT, fyT float64, edgeNum int) (fdot, edgeNormMag float64) {
			en1, en2 := edgeNorm1[edgeNum], edgeNorm2[edgeNum]
			fdot = fxT*en1 + fyT*en2
			fmt.Printf("edge norm[%d] = [%8.5f,%8.5f], fdot[%8.5f,%8.5f] = %8.5f\n", edgeNum, en1, en2, fxT, fyT, fdot)
			edgeNormMag = math.Sqrt(en1*en1 + en2*en2)
			return
		}
		for n := 0; n < Np; n++ {
			ind := n + k*Np
			fxT := Jdet * (Jinv[0]*fxD[n] + Jinv[1]*fyD[n])
			fyT := Jdet * (Jinv[2]*fxD[n] + Jinv[3]*fyD[n])
			switch rt.GetTermType(n) {
			case All:
				panic("bad input")
			case InteriorR:
				// Unit vector is [1,0]
				fxpD[ind] = fxT
			case InteriorS:
				// Unit vector is [0,1]
				fypD[ind] = fyT
			case Edge1:
				// Edge1: Unit vector is [1/sqrt(2), 1/sqrt(2)]
				fdot, edgeNormMag := getEdgeQuants(fxT, fyT, 0)
				//edgeNormMag = 1
				//fxpD[ind] = oosr2 * fdot * edgeNormMag
				fxpD[ind] = oosr2 * fdot * edgeNormMag
				fypD[ind] = fxpD[ind]
			case Edge2:
				// Edge2: Unit vector is [-1,0]
				fdot, edgeNormMag := getEdgeQuants(fxT, fyT, 1)
				//edgeNormMag = 1
				fxpD[ind] = -fdot * edgeNormMag
			case Edge3:
				// Edge3: // Unit vector is [0,-1]
				fdot, edgeNormMag := getEdgeQuants(fxT, fyT, 2)
				//edgeNormMag = 1
				fypD[ind] = -fdot * edgeNormMag
			}
		}
	}
	return
}
