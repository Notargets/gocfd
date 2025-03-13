package DG2D

import (
	"fmt"
	"math"
	"strings"

	"github.com/notargets/avs/geometry"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/readfiles"

	"github.com/notargets/gocfd/utils"
)

type DFR2D struct {
	N               int
	SolutionElement *LagrangeElement2D
	// FluxElement        *RTBasis2DSimplexLegacy
	FluxElement *RTElement
	// Interpolates from the interior (solution) points to graph points
	GraphMesh      geometry.TriMesh
	GraphInterp    utils.Matrix
	GraphInterpMod utils.Matrix // Modulated
	FluxEdgeInterp utils.Matrix // Interpolates only from interior to the edge points in the flux element
	// FilterMod          utils.Matrix // Solution Modulation filter damps high modes
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

func NewDFR2D(N int, verbose bool, meshFileO ...string) (dfr *DFR2D) {
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
	GraphR, GraphS := dfr.GetRSForGraphMesh()
	Nu, p := 0.1, 3.29
	// dfr.GraphInterp = dfr.SolutionBasis.GetModInterpMatrix(GraphR, GraphS, Nu, p, 1)
	dfr.GraphInterpMod = dfr.SolutionBasis.GetModInterpMatrix(GraphR, GraphS,
		Nu, p, 1)
	dfr.GraphInterp = dfr.SolutionBasis.GetInterpMatrix(GraphR, GraphS)
	// dfr.FluxEdgeInterp = dfr.SolutionBasis.GetInterpMatrix(RFlux, SFlux) // Interpolation matrix across three edges
	dfr.FluxEdgeInterp = dfr.SolutionBasis.GetModInterpMatrix(RFlux, SFlux, Nu, p, 1) // Interpolation matrix across three edges
	// dfr.FluxEdgeInterp = dfr.SolutionBasis.GetModInterpMatrix(RFlux, SFlux, Nu, p, 1) // Interpolation matrix across three edges
	// dfr.FilterMod = dfr.SolutionBasis.GetModInterpMatrix(dfr.SolutionElement.R, dfr.SolutionElement.S, Nu, p, 1)
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
		dfr.ProcessGeometry(dfr.VX, dfr.VY, EToV, dfr.BCEdges)
		dfr.GraphMesh = dfr.CreateAVSGraphMesh()
	}
	return
}

func (dfr *DFR2D) ProcessGeometry(VX, VY utils.Vector,
	EToV utils.Matrix, BCEdges types.BCMAP) {
	// dfr.BCEdges.Print()
	dfr.Tris = NewTriangulation(VX, VY, EToV, BCEdges)
	// Build connectivity matrices
	dfr.FluxX, dfr.FluxY =
		CalculateElementLocalGeometry(EToV, VX, VY,
			dfr.FluxElement.R, dfr.FluxElement.S)
	dfr.SolutionX, dfr.SolutionY =
		CalculateElementLocalGeometry(EToV, VX, VY,
			dfr.SolutionElement.R, dfr.SolutionElement.S)
	// Calculate RT based derivative metrics for use in calculating Dx and Dy
	// using the RT element
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
