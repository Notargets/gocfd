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
		dfr.Tris = NewTriangulation(dfr.VX, dfr.VY, EToV, dfr.BCType)
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
		//xr, yr := 0.5*(v2x-v1x), 0.5*(v2y-v1y)
		//xs, ys := 0.5*(v3x-v1x), 0.5*(v3y-v1y)
		xr, yr := 0.5*(v2x-v1x), 0.5*(v2y-v1y)
		xs, ys := 0.5*(v3x-v1x), 0.5*(v3y-v1y)
		// Jacobian is [xr, xs]
		//             [yr, ys]
		jd := Jd[k*4:]
		jd[0], jd[1], jd[2], jd[3] = xr, xs, yr, ys
		Jdetd[k] = xr*ys - xs*yr
		oodet := 1. / Jdetd[k]
		jdInv := JdInv[k*4:]
		// Inverse Jacobian is:
		// (1/(xr*ys-xs*yr)) *
		//             [ ys,-xs]
		//             [-yr, xr]
		jdInv[0], jdInv[1], jdInv[2], jdInv[3] = ys, -xs, -yr, xr
		for i := 0; i < 4; i++ {
			jdInv[i] *= oodet
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
			//fn := dfr.FaceNorm.Row(int(triNum)).Data()[0:3]
			k := int(triNum)
			fnD1, fnD2 := dfr.FaceNorm[0].Data(), dfr.FaceNorm[1].Data()
			x1, x2 := dfr.Tris.GetEdgeCoordinates(en, bool(e.ConnectedTriDirection[connNum]), dfr.VX, dfr.VY)
			dx, dy := x2[0]-x1[0], x2[1]-x1[1]
			nx, ny := -dy, dx
			switch e.ConnectedTriEdgeNumber[connNum] {
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
				// Edge1:
				fpD[ind] = oosr2 * (fT[0] + fT[1])
			case Edge2:
				// Edge2:
				//fpD[ind] = fdot2 * edgeNormMag2
				fpD[ind] = -fT[0]
			case Edge3:
				// Edge3:
				//fpD[ind] = fdot3 * edgeNormMag3
				fpD[ind] = -fT[1]
			}
		}
	}
	return
}
