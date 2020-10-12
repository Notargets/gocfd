package DG2D

import "github.com/notargets/gocfd/utils"

type DFR2D struct {
	N, Np        int
	R, S         utils.Vector // Solution point locations in unit triangle
	FluxR, FluxS utils.Vector // Flux (face) point locations in unit triangle
	V, Vinv      utils.Matrix // Vandermonde and inverse for solution polynomial (order N), defined on solution points
	Dr, Ds       utils.Matrix // Differentiation matrices in R,S directions for solution points
}

func NewDFR2D(N int) (dfr *DFR2D) {
	var (
		err error
	)
	dfr = &DFR2D{
		N:  N,
		Np: (N + 1) * (N + 2) / 2,
	}
	// Compute nodal set
	dfr.R, dfr.S = NodesEpsilon(N)
	// Build reference element matrices
	dfr.V = Vandermonde2D(N, dfr.R, dfr.S)
	if dfr.Vinv, err = dfr.V.Inverse(); err != nil {
		panic(err)
	}
	// Initialize the (r,s) differentiation matrices on the simplex, evaluated at (r,s) at order N
	Vr, Vs := GradVandermonde2D(N, dfr.R, dfr.S)
	dfr.Dr = Vr.Mul(dfr.Vinv)
	dfr.Ds = Vs.Mul(dfr.Vinv)

	// Mark fields read only
	dfr.Dr.SetReadOnly("Dr")
	dfr.Ds.SetReadOnly("Ds")
	dfr.V.SetReadOnly("V")
	dfr.Vinv.SetReadOnly("Vinv")
	return
}
