package DG2D

import (
	"math"

	"github.com/notargets/gocfd/utils"
)

type Elements2D struct {
	K, Np, NFaces                     int
	R, VX, VY, VZ, FMask              utils.Vector
	EToV, EToE, EToF                  utils.Matrix
	BCType                            utils.Matrix
	X, Dr, Rx, FScale, NX, LIFT       utils.Matrix
	V, Vinv                           utils.Matrix
	VmapM, VmapP, VmapB, VmapI, VmapO utils.Index
	MapB, MapI, MapO                  utils.Index
	Cub                               *Cubature
}

func NewElements2D(N, K int, meshFile string, plotMesh bool) (el *Elements2D) {
	var (
		NFaces = 3
		// choose order to integrate exactly
		CubatureOrder = int(math.Floor(2.0 * float64(N+1) * 3.0 / 2.0))
		NGauss        = int(math.Floor(2.0 * float64(N+1)))
	)
	_ = NGauss
	/*
	  // build cubature node data for all elements
	  CubatureVolumeMesh2D(CubatureOrder);

	  // build Gauss node data for all element faces
	  GaussFaceMesh2D(NGauss);

	  Resize_cub();           // resize cubature arrays
	  MapGaussFaceData();     // {nx = gauss.nx}, etc.
	  PreCalcBdryData();      // gmapB = concat(mapI, mapO), etc.
	*/
	// N is the polynomial degree, Np is the number of interpolant points = N+1
	el = &Elements2D{
		K:      K,
		Np:     (N + 1) * (N + 2) / 2,
		NFaces: NFaces,
	}
	el.ReadGambit2d(meshFile, plotMesh)
	el.NewCube2D(CubatureOrder)
	el.Startup2D()
	return
}
