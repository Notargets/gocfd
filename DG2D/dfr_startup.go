package DG2D

import (
	"fmt"
	"math"

	graphics2D "github.com/notargets/avs/geometry"
	"github.com/notargets/gocfd/geometry2D"

	"github.com/notargets/gocfd/utils"
)

type DFR2D struct {
	N                    int
	SolutionElement      *LagrangeElement2D
	FluxElement          *RTElement
	FluxInterpMatrix     utils.Matrix
	FluxEdgeInterpMatrix utils.Matrix // A pre-calculated interpolation matrix covering all Flux (edge) points in K elements
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

func NewDFR2D(N int, plotMesh bool, meshFileO ...string) (dfr *DFR2D) {
	if N < 1 {
		panic(fmt.Errorf("Polynomial order must be >= 1, have %d", N))
	}
	le := NewLagrangeElement2D(N, Epsilon)
	rt := NewRTElement(N+1, le.R, le.S)
	RFlux := utils.NewVector(rt.Nedge*3, rt.GetEdgeLocations(rt.R)) // For the Interpolation matrix across three edges
	SFlux := utils.NewVector(rt.Nedge*3, rt.GetEdgeLocations(rt.S)) // For the Interpolation matrix across three edges
	dfr = &DFR2D{
		N:                    N,
		SolutionElement:      le,
		FluxElement:          rt,
		FluxInterpMatrix:     le.Simplex2DInterpolatingPolyMatrix(rt.R, rt.S),   // Interpolation matrix for flux nodes
		FluxEdgeInterpMatrix: le.Simplex2DInterpolatingPolyMatrix(RFlux, SFlux), // Interpolation matrix across three edges
	}
	if len(meshFileO) != 0 {
		var EToV utils.Matrix
		dfr.K, dfr.VX, dfr.VY, EToV, dfr.BCType = ReadGambit2d(meshFileO[0], false)
		dfr.Tris = NewTriangulation(dfr.VX, dfr.VY, EToV, dfr.BCType)
		// Build connectivity matrices
		dfr.FluxX, dfr.FluxY =
			CalculateElementLocalGeometry(dfr.Tris.EToV, dfr.VX, dfr.VY, dfr.FluxElement.R, dfr.FluxElement.S)
		dfr.SolutionX, dfr.SolutionY =
			CalculateElementLocalGeometry(dfr.Tris.EToV, dfr.VX, dfr.VY, dfr.SolutionElement.R, dfr.SolutionElement.S)
		if plotMesh {
			PlotMesh(dfr.VX, dfr.VY, EToV, dfr.BCType, dfr.SolutionX, dfr.SolutionY, true)
			utils.SleepFor(50000)
		}
		dfr.FluxX.SetReadOnly("FluxX")
		dfr.FluxY.SetReadOnly("FluxY")
		dfr.SolutionX.SetReadOnly("SolutionX")
		dfr.SolutionY.SetReadOnly("SolutionY")
		dfr.CalculateJacobian()
		dfr.CalculateFaceNorms()
	}
	return
}

func (dfr *DFR2D) GetJacobian(k int) (J, Jinv []float64, Jdet float64) {
	J = dfr.J.Row(k).Data()[0:4]
	Jinv = dfr.Jinv.Row(k).Data()[0:4]
	Jdet = dfr.Jdet.Row(k).Data()[0]
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
		//xr, yr := 0.5*(v2x-v1x), 0.5*(v2y-v1y)
		//xs, ys := 0.5*(v3x-v1x), 0.5*(v3y-v1y)
		xr, yr := 0.5*(v2x-v1x), 0.5*(v2y-v1y)
		xs, ys := 0.5*(v3x-v1x), 0.5*(v3y-v1y)
		// Jacobian is [xr, xs]
		//             [yr, ys]
		ind := k * 4
		Jd[ind+0], Jd[ind+1], Jd[ind+2], Jd[ind+3] = xr, xs, yr, ys
		Jdetd[k] = xr*ys - xs*yr
		// Inverse Jacobian is:
		// [ ys,-xs] * (1/(xr*ys-xs*yr))
		// [-yr, xr]
		JdInv[ind+0], JdInv[ind+1], JdInv[ind+2], JdInv[ind+3] = ys, -xs, -yr, xr
		for i := 0; i < 4; i++ {
			JdInv[ind+i] /= Jdetd[k]
		}
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
		for connNum, triNum := range e.ConnectedTris {
			k := int(triNum)
			fnD1, fnD2 := dfr.FaceNorm[0].Data(), dfr.FaceNorm[1].Data()
			x1, x2 := GetEdgeCoordinates(en, bool(e.ConnectedTriDirection[connNum]), dfr.VX, dfr.VY)
			dx, dy := x2[0]-x1[0], x2[1]-x1[1]
			nx, ny := -dy, dx
			edgeNum := e.ConnectedTriEdgeNumber[connNum].Index() // one of [0,1,2], aka "First", "Second", "Third"
			fnD1[edgeNum+k*3], fnD2[edgeNum+k*3] = nx, ny
		}
	}
	dfr.FaceNorm[0].SetReadOnly("FaceNorm1")
	dfr.FaceNorm[1].SetReadOnly("FaceNorm2")
}

func (dfr *DFR2D) ProjectFluxOntoRTSpace(Fx, Fy utils.Matrix) (Fp utils.Matrix) {
	var (
		Np = dfr.FluxElement.Np
		K  = dfr.K
		rt = dfr.FluxElement
	)
	Fp = utils.NewMatrix(K, Np)
	for k := 0; k < K; k++ {
		var (
			Jdet     = dfr.Jdet.Row(k).Data()[0]
			Jinv     = dfr.Jinv.Row(k).Data()[0:4]
			fxD, fyD = Fx.Data(), Fy.Data()
			fpD      = Fp.Data()
		)
		for n := 0; n < Np; n++ {
			ind := n + k*Np
			fT := [2]float64{Jdet * (Jinv[0]*fxD[ind] + Jinv[1]*fyD[ind]), Jdet * (Jinv[2]*fxD[ind] + Jinv[3]*fyD[ind])}
			oosr2 := 1 / math.Sqrt(2)
			switch rt.GetTermType(n) {
			case All:
				panic("bad input")
			case InteriorR:
				// Unit vector is [1,0]
				fpD[ind] = fT[0]
			case InteriorS:
				// Unit vector is [0,1]
				fpD[ind] = fT[1]
			case Edge1:
				// Edge3:
				fpD[ind] = -fT[1]
			case Edge2:
				// Edge1:
				fpD[ind] = oosr2 * (fT[0] + fT[1])
			case Edge3:
				// Edge2:
				fpD[ind] = -fT[0]
			}
		}
	}
	return
}

func (dfr *DFR2D) ConvertScalarToOutputMesh(f utils.Matrix) (fI []float32) {
	/*
				Input f contains the function data to be associated with the output mesh
				The input dimensions of f are: f(Np, K), where Np is the number of RT nodes and K is the element count

				Output fI contains the function data in the same order as the vertices of the output mesh
		    	The corners of each element are formed by averaging the nearest two edge values
	*/
	var (
		fD     = f.Data()
		Kmax   = dfr.K
		Nint   = dfr.FluxElement.Nint
		Nedge  = dfr.FluxElement.Nedge
		NpFlux = dfr.FluxElement.Np
		Np     = NpFlux - Nint + 3 // Subtract Nint to remove the dup pts and add 3 for the verts
	)
	Ind := func(k, i, Kmax int) (ind int) {
		ind = k + i*Kmax
		return
	}
	fI = make([]float32, Kmax*Np)
	for k := 0; k < Kmax; k++ {
		var (
			edge [3][2]float32
		)
		for ii := 0; ii < 3; ii++ {
			beg := 2*Nint + ii*Nedge
			end := beg + Nedge - 1
			ie0, ie1 := Ind(k, beg, Kmax), Ind(k, end, Kmax)
			// [ii][0] is the first point on the edge, [ii][1] is the second
			edge[ii][0], edge[ii][1] = float32(fD[ie0]), float32(fD[ie1])
		}
		for ii := 0; ii < Np; ii++ {
			ind := Ind(k, ii, Kmax)
			switch {
			case ii < 3:
				// Create values for each corner by averaging the nodes opposite each
				fI[ind] = 0.5 * (edge[(ii+2)%3][1] + edge[ii][0])
			case ii >= 3:
				indFlux := Ind(k, ii-3+Nint, Kmax) // Refers to the nodes, skipping the first Nint repeated points
				fI[ind] = float32(fD[indFlux])
			}
		}
	}
	return
}

func (dfr *DFR2D) OutputMesh() (gm *graphics2D.TriMesh) {
	/*
				For each of K elements, the layout of the output mesh is:
				Indexed geometry:
					Three vertices of the base triangle, followed by RT node points, excluding duplicated Nint pts
		    	Triangles:
					Some number of triangles, defined by indexing into the Indexed Geometry / function values

				******************************************************************************************************
				Given there is no function data for the corners of the RT element, these will have to be supplied when
				constructing the indexed function data to complement this output mesh
	*/
	// Triangulate the unit RT triangle: start with the bounding triangle, which includes the corners to constrain
	// the Delaunay triangulation
	var (
		Kmax   = dfr.K
		Nint   = dfr.FluxElement.Nint
		NpFlux = dfr.FluxElement.Np
	)
	Ind := func(k, i, Kmax int) (ind int) {
		ind = k + i*Kmax
		return
	}
	R := []float64{-1, 1, -1} // Vertices of unit triangle
	S := []float64{-1, -1, 1}
	tm := geometry2D.NewTriMesh(R, S)
	tri := &geometry2D.Tri{}
	tri.AddEdge(tm.NewEdge([2]int{0, 1}, true))
	e2 := tm.NewEdge([2]int{1, 2}, true)
	tri.AddEdge(e2)
	tri.AddEdge(tm.NewEdge([2]int{2, 0}, true))
	tm.AddBoundingTriangle(tri)
	// Now we add points to incrementally define the triangulation
	for i := Nint; i < NpFlux; i++ {
		r := dfr.FluxElement.R.Data()[i]
		s := dfr.FluxElement.S.Data()[i]
		tm.AddPoint(r, s)
	}
	gmB := tm.ToGraphMesh()
	gm = &gmB

	// Build the X,Y coordinates to support the triangulation index
	Np := NpFlux - Nint + 3 // Subtract Nint to remove the dup pts and add 3 for the verts
	VX, VY := utils.NewMatrix(Np, Kmax), utils.NewMatrix(Np, Kmax)
	vxd, vyd := VX.Data(), VY.Data()
	for k := 0; k < Kmax; k++ {
		verts := dfr.Tris.GetTriVerts(uint32(k))
		for ii := 0; ii < Np; ii++ {
			ind := Ind(k, ii, Kmax)
			switch {
			case ii < 3:
				vxd[ind], vyd[ind] = dfr.VX.Data()[verts[ii]], dfr.VY.Data()[verts[ii]]
			case ii >= 3:
				indFlux := Ind(k, ii-3+Nint, Kmax) // Refers to the nodes, skipping the first Nint repeated points
				vxd[ind], vyd[ind] = dfr.FluxX.Data()[indFlux], dfr.FluxY.Data()[indFlux]
			}
		}
	}
	//fmt.Println(VX.Transpose().Print("VX_Out"))
	//fmt.Println(VY.Transpose().Print("VY_Out"))

	// Now replicate the triangle mesh for all triangles
	baseTris := gm.Triangles
	gm.Triangles = make([]graphics2D.Triangle, Kmax*len(baseTris))
	for k := 0; k < Kmax; k++ {
		for i, tri := range baseTris {
			newTri := graphics2D.Triangle{Nodes: tri.Nodes}
			for ii := 0; ii < 3; ii++ {
				newTri.Nodes[ii] = int32(Ind(k, int(newTri.Nodes[ii]), Kmax))
			}
			gm.Triangles[Ind(k, i, Kmax)] = newTri
		}
	}
	gm.BaseGeometryClass.Geometry = make([]graphics2D.Point, Kmax*Np)
	for k := 0; k < Kmax; k++ {
		for ii := 0; ii < Np; ii++ {
			ind := Ind(k, ii, Kmax)
			gm.BaseGeometryClass.Geometry[ind] = graphics2D.Point{X: [2]float32{float32(vxd[ind]), float32(vyd[ind])}}
		}
	}
	//gm.Attributes = make([][]float32, len(gm.Triangles)) // Empty attributes
	gm.Attributes = nil
	return
}
