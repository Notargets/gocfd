package DG2D

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strings"
	"time"

	"github.com/notargets/gocfd/InputParameters"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/readfiles"

	graphics2D "github.com/notargets/avs/geometry"
	"github.com/notargets/gocfd/geometry2D"

	"github.com/notargets/gocfd/utils"
)

type DFR2D struct {
	N               int
	SolutionElement *LagrangeElement2D
	// FluxElement        *RTBasis2DSimplexLegacy
	FluxElement        *RTElement
	FluxInterp         utils.Matrix // Interpolates from the interior (solution) points to all of the flux points
	FluxEdgeInterp     utils.Matrix // Interpolates only from interior to the edge points in the flux element
	FluxDr, FluxDs     utils.Matrix // Derivatives from the interior (solution) points to all of the flux points
	DXMetric, DYMetric utils.Matrix // X and Y derivative metrics, multiply by scalar field values to create DOF for RT
	// Mesh Parameters
	K                    int          // Number of elements (triangles) in mesh
	VX, VY               utils.Vector // X,Y vertex points in mesh (vertices)
	BCEdges              types.BCMAP
	FluxX, FluxY         utils.Matrix    // Flux Element local coordinates
	SolutionX, SolutionY utils.Matrix    // Solution Element local coordinates
	Tris                 *Triangulation  // Triangle mesh and edge/face structures
	J, Jinv, Jdet        utils.Matrix    // Mesh Transform Jacobian, Kx[4], each K element has a 2x2 matrix, det is |J|
	FaceNorm             [2]utils.Matrix // Face normal (normalized), Kx3
	IInII                utils.Matrix    // Mag face normal divided by unit triangle face norm mag, Kx3 dimension
	EdgeNumber           []types.EdgeKey // Edge number for each edge, used to index into edge structures, Kx3 dimension
	SolutionBasis        *JacobiBasis2D
}

func NewDFR2D(N int, pm *InputParameters.PlotMeta, verbose bool, meshFileO ...string) (dfr *DFR2D) {
	if N < 0 {
		panic(fmt.Errorf("Polynomial order must be >= 0, have %d", N))
	}
	le := NewLagrangeElement2D(N, Epsilon)
	rt := NewRTElement(N+1, SimplexRTBasis)
	RFlux := utils.NewVector(rt.NpEdge*3, rt.GetEdgeLocations(rt.R.DataP)) // For the Interpolation matrix across three edges
	SFlux := utils.NewVector(rt.NpEdge*3, rt.GetEdgeLocations(rt.S.DataP)) // For the Interpolation matrix across three edges
	dfr = &DFR2D{
		N:               N,
		SolutionElement: le,
		FluxElement:     rt,
	}
	// Get interpolation matrix for edges using a basis on solution points at polynomial degree le.N
	if verbose {
		fmt.Printf("Using 2D Jacobi OrthoNormal Polynomial Basis\n")
	}
	dfr.SolutionBasis = NewJacobiBasis2D(le.N, le.R, le.S, 0, 0)
	// Get the interpolation matrices that interpolate the whole RT element and just the edges using solution points
	dfr.FluxInterp = dfr.SolutionBasis.GetInterpMatrix(rt.R, rt.S)       // Interpolation matrix for flux nodes
	dfr.FluxEdgeInterp = dfr.SolutionBasis.GetInterpMatrix(RFlux, SFlux) // Interpolation matrix across three edges
	dfr.FluxDr, dfr.FluxDs = le.GetDerivativeMatrices(rt.R, rt.S)
	if len(meshFileO) != 0 {
		var EToV utils.Matrix
		t := getFileTypeFromExtension(meshFileO[0])
		switch t {
		case GAMBIT_FILE:
			dfr.K, dfr.VX, dfr.VY, EToV, dfr.BCEdges =
				readfiles.ReadGambit2d(meshFileO[0], verbose)
		case SU2_FILE:
			dfr.K, dfr.VX, dfr.VY, EToV, dfr.BCEdges =
				readfiles.ReadSU2(meshFileO[0], verbose)
		}
		// dfr.BCEdges.Print()
		dfr.Tris = NewTriangulation(dfr.VX, dfr.VY, EToV, dfr.BCEdges)
		// Build connectivity matrices
		dfr.FluxX, dfr.FluxY =
			CalculateElementLocalGeometry(dfr.Tris.EToV, dfr.VX, dfr.VY, dfr.FluxElement.R, dfr.FluxElement.S)
		dfr.SolutionX, dfr.SolutionY =
			CalculateElementLocalGeometry(dfr.Tris.EToV, dfr.VX, dfr.VY, dfr.SolutionElement.R, dfr.SolutionElement.S)
		if pm.PlotMesh {
			readfiles.PlotMesh(dfr.VX, dfr.VY, EToV, dfr.SolutionX, dfr.SolutionY, true, pm)
			fmt.Println("Number of elements in mesh K = ", dfr.K)
			fmt.Println("Press 'Enter' to exit...")
			bufio.NewReader(os.Stdin).ReadBytes('\n')
			time.Sleep(30 * time.Second)
			os.Exit(0)
		}
		// Calculate RT based derivative metrics for use in calculating Dx and Dy using the RT element
		dfr.CalculateJacobian()
		dfr.CalculateFaceNorms()
		dfr.CalculateRTBasedDerivativeMetrics()
		dfr.DXMetric.SetReadOnly("DXMetric")
		dfr.DYMetric.SetReadOnly("DYMetric")
		dfr.J.SetReadOnly("GeometricJacobian")
		dfr.Jdet.SetReadOnly("GeometricJacobianDeterminant")
		dfr.Jinv.SetReadOnly("GeometricJacobianInverse")
		dfr.FluxX.SetReadOnly("FluxX")
		dfr.FluxY.SetReadOnly("FluxY")
		dfr.SolutionX.SetReadOnly("SolutionX")
		dfr.SolutionY.SetReadOnly("SolutionY")
	}
	return
}

type MeshFileType uint8

const (
	GAMBIT_FILE MeshFileType = iota
	SU2_FILE
)

func getFileTypeFromExtension(fileName string) (t MeshFileType) {
	var (
		err error
	)
	fileName = strings.Trim(fileName, " ")
	l := len(fileName)
	if l < 4 || fileName[l-4] != '.' {
		err = fmt.Errorf("unable to determine file type from name: %s", fileName)
		panic(err)
	}
	ext := fileName[l-3:]
	switch ext {
	case "neu": // Gambit neutral file
		return GAMBIT_FILE
	case "su2": // SU2 file
		return SU2_FILE
	default:
		err = fmt.Errorf("unsupported file type: %s", fileName)
		panic(err)
	}
}

func (dfr *DFR2D) GetJacobian(k int) (J, Jinv []float64, Jdet float64) {
	J = dfr.J.DataP[k*4 : (k+1)*4]
	Jinv = dfr.Jinv.DataP[k*4 : (k+1)*4]
	Jdet = dfr.Jdet.At(k, 0)
	return
}

func (dfr *DFR2D) CalculateRTBasedDerivativeMetrics() {
	/*
		Multiply either of these matrices by a scalar solution value at NpFlux points, then multiply
		that result by rt.Div to get either the X or Y derivative of that field
		Note that you can specify arbitrary solution values at the edges to establish continuity and you will
		still have accuracy of the X or Y derivative equal to the interior element polynomial degree
	*/
	var (
		NpFlux, NpInt, Kmax = dfr.FluxElement.Np, dfr.FluxElement.NpInt, dfr.K
		NpEdge              = dfr.FluxElement.NpEdge
		Jinv, Jdet          = dfr.Jinv, dfr.Jdet
	)
	dfr.DXMetric, dfr.DYMetric = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
	RTDXmd, RTDYmd := dfr.DXMetric.DataP, dfr.DYMetric.DataP
	for k := 0; k < Kmax; k++ {
		var (
			JinvD = Jinv.DataP[4*k : 4*(k+1)]
		)
		for i := 0; i < NpInt; i++ {
			ind := k + i*Kmax
			ind2 := k + (i+NpInt)*Kmax
			RTDXmd[ind], RTDXmd[ind2] = JinvD[0], JinvD[2]
			RTDYmd[ind], RTDYmd[ind2] = JinvD[1], JinvD[3]
		}
	}
	var (
		fnd0, fnd1 = dfr.FaceNorm[0].DataP, dfr.FaceNorm[1].DataP
	)
	for k := 0; k < Kmax; k++ {
		ooJd := 1. / Jdet.DataP[k]
		for fn := 0; fn < 3; fn++ {
			faceInd := k + Kmax*fn
			norm := [2]float64{fnd0[faceInd], fnd1[faceInd]}
			IInII := dfr.IInII.DataP[faceInd]
			for i := 0; i < NpEdge; i++ {
				shift := fn * NpEdge
				indFull := k + (2*NpInt+i+shift)*Kmax
				RTDXmd[indFull], RTDYmd[indFull] = ooJd*norm[0]*IInII, ooJd*norm[1]*IInII
			}
		}
	}
}

func (dfr *DFR2D) CalculateJacobian() {
	Jd := make([]float64, 4*dfr.K)
	Jdetd := make([]float64, dfr.K)
	JdInv := make([]float64, 4*dfr.K)
	for k := 0; k < dfr.K; k++ {
		tri := dfr.Tris.EToV.Row(k).DataP
		v := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
		v1x, v2x, v3x := dfr.VX.AtVec(v[0]), dfr.VX.AtVec(v[1]), dfr.VX.AtVec(v[2])
		v1y, v2y, v3y := dfr.VY.AtVec(v[0]), dfr.VY.AtVec(v[1]), dfr.VY.AtVec(v[2])
		// xr, yr := 0.5*(v2x-v1x), 0.5*(v2y-v1y)
		// xs, ys := 0.5*(v3x-v1x), 0.5*(v3y-v1y)
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
	var (
		Kmax = dfr.K
	)
	dfr.FaceNorm[0], dfr.FaceNorm[1] = utils.NewMatrix(3, Kmax), utils.NewMatrix(3, Kmax)
	dfr.IInII = utils.NewMatrix(3, Kmax)
	dfr.EdgeNumber = make([]types.EdgeKey, 3*Kmax)
	for en, e := range dfr.Tris.Edges {
		for connNum := 0; connNum < int(e.NumConnectedTris); connNum++ {
			k := int(e.ConnectedTris[connNum])
			edgeNum := e.ConnectedTriEdgeNumber[connNum].Index() // one of [0,1,2], aka "First", "Second", "Third"
			ind := k + Kmax*edgeNum
			dfr.IInII.DataP[ind] = e.IInII[connNum]
			dfr.EdgeNumber[ind] = en
			fnD1, fnD2 := dfr.FaceNorm[0].DataP, dfr.FaceNorm[1].DataP
			x1, x2 := GetEdgeCoordinates(en, bool(e.ConnectedTriDirection[connNum]), dfr.VX, dfr.VY)
			dx, dy := x2[0]-x1[0], x2[1]-x1[1]
			oonorm := 1. / math.Sqrt(dx*dx+dy*dy)
			nx, ny := -dy*oonorm, dx*oonorm
			fnD1[ind], fnD2[ind] = nx, ny
		}
	}
	dfr.FaceNorm[0].SetReadOnly("FaceNorm1")
	dfr.FaceNorm[1].SetReadOnly("FaceNorm2")
	dfr.IInII.SetReadOnly("IInII")
}

func (dfr *DFR2D) ProjectFluxOntoRTSpace(Fx, Fy utils.Matrix) (Fp utils.Matrix) {
	var (
		Np       = dfr.FluxElement.Np
		K        = dfr.K
		rt       = dfr.FluxElement
		JdetD    = dfr.Jdet.DataP
		JinvD    = dfr.Jinv.DataP
		fxD, fyD = Fx.DataP, Fy.DataP
	)
	Fp = utils.NewMatrix(K, Np)
	fpD := Fp.DataP
	for k := 0; k < K; k++ {
		var (
			Jdet = JdetD[k]
			Jinv = JinvD[4*k : 4*k+4]
		)
		for n := 0; n < Np; n++ {
			ind := n + k*Np
			fT := [2]float64{Jdet * (Jinv[0]*fxD[ind] + Jinv[1]*fyD[ind]), Jdet * (Jinv[2]*fxD[ind] + Jinv[3]*fyD[ind])}
			oosr2 := 1 / math.Sqrt(2)
			switch rt.getLocationType(n) {
			case All:
				panic("bad input")
			case E4:
				// Unit vector is [1,0]
				fpD[ind] = fT[0]
			case E5:
				// Unit vector is [0,1]
				fpD[ind] = fT[1]
			case E1:
				// E3:
				fpD[ind] = -fT[1]
			case E2:
				// E1:
				fpD[ind] = oosr2 * (fT[0] + fT[1])
			case E3:
				// E2:
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
		fD     = f.DataP
		Kmax   = dfr.K
		Nint   = dfr.FluxElement.NpInt
		Nedge  = dfr.FluxElement.NpEdge
		NpFlux = dfr.FluxElement.Np
		Np     = NpFlux - Nint + 3 // Subtract NpInt to remove the dup pts and add 3 for the verts
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
				indFlux := Ind(k, ii-3+Nint, Kmax) // Refers to the nodes, skipping the first NpInt repeated points
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
					Three vertices of the base triangle, followed by RT node points, excluding duplicated NpInt pts
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
		Nint   = dfr.FluxElement.NpInt
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
		r := dfr.FluxElement.R.DataP[i]
		s := dfr.FluxElement.S.DataP[i]
		tm.AddPoint(r, s)
	}
	gmB := tm.ToGraphMesh()
	gm = &gmB

	// Build the X,Y coordinates to support the triangulation index
	Np := NpFlux - Nint + 3 // Subtract NpInt to remove the dup pts and add 3 for the verts
	//	fmt.Printf("Size of constructed element: Np=%d\n", Np)
	//	fmt.Printf("Size of original flux element: NpFlux=%d\n", NpFlux)
	//	fmt.Printf("Size of original flux element interior: NpInt=%d\n", Nint)
	VX, VY := utils.NewMatrix(Np, Kmax), utils.NewMatrix(Np, Kmax)
	vxd, vyd := VX.DataP, VY.DataP
	for k := 0; k < Kmax; k++ {
		verts := dfr.Tris.GetTriVerts(uint32(k))
		for ii := 0; ii < Np; ii++ {
			ind := Ind(k, ii, Kmax)
			switch {
			case ii < 3:
				vxd[ind], vyd[ind] = dfr.VX.DataP[verts[ii]], dfr.VY.DataP[verts[ii]]
			case ii >= 3:
				indFlux := Ind(k, ii-3+Nint, Kmax) // Refers to the nodes, skipping the first NpInt repeated points
				vxd[ind], vyd[ind] = dfr.FluxX.DataP[indFlux], dfr.FluxY.DataP[indFlux]
			}
		}
	}
	// fmt.Println(VX.Transpose().Print("VX_Out"))
	// fmt.Println(VY.Transpose().Print("VY_Out"))

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
	// gm.Attributes = make([][]float32, len(gm.Triangles)) // Empty attributes
	gm.Attributes = nil
	return
}
