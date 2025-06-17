package gonudg

// INDEXING NOTE: Original C++ code uses 1-based indexing to emulate Matlab behavior.
// This Go port uses standard 0-based indexing. Example conversions:
//   C++: sk = 1; V3D(All,sk) = ...    ->    Go: sk = 0; V3D.SetCol(sk, ...)
//   C++: Fmask[1] (first face)        ->    Go: Fmask[0] (first face)
// The indexing has been correctly translated throughout this port.

import (
	"github.com/notargets/gocfd/utils"
)

// Vandermonde3D initializes the 3D Vandermonde Matrix V_{ij} = phi_j(r_i, s_i, t_i)
// This is the 0-based index version of the C++ Vandermonde3D function
func Vandermonde3D(N int, r, s, t []float64) utils.Matrix {
	Np := len(r)
	Ncol := (N + 1) * (N + 2) * (N + 3) / 6

	// Initialize the Vandermonde matrix
	V3D := utils.NewMatrix(Np, Ncol)

	// Transfer to (a,b,c) coordinates
	a, b, c := RSTtoABC(r, s, t)

	// Build the Vandermonde matrix
	sk := 0 // 0-based column index
	for i := 0; i <= N; i++ {
		for j := 0; j <= N-i; j++ {
			for k := 0; k <= N-i-j; k++ {
				// Evaluate basis function at all points
				col := Simplex3DP(a, b, c, i, j, k)

				// Copy to matrix column
				V3D.SetCol(sk, col)
				sk++
			}
		}
	}

	return V3D
}

// GradVandermonde3D builds the gradient Vandermonde matrices
// Returns Vr, Vs, Vt where (Vr)_{ij} = dphi_j/dr at point i
func GradVandermonde3D(N int, r, s, t []float64) (Vr, Vs, Vt utils.Matrix) {
	Np := len(r)
	Ncol := (N + 1) * (N + 2) * (N + 3) / 6

	// Initialize the gradient matrices
	Vr = utils.NewMatrix(Np, Ncol)
	Vs = utils.NewMatrix(Np, Ncol)
	Vt = utils.NewMatrix(Np, Ncol)

	// Build the gradient Vandermonde matrices
	sk := 0 // 0-based column index
	for i := 0; i <= N; i++ {
		for j := 0; j <= N-i; j++ {
			for k := 0; k <= N-i-j; k++ {
				// Evaluate gradient of basis function at all points
				dr, ds, dt := GradSimplex3DP(r, s, t, i, j, k)

				// Copy to matrix columns
				Vr.SetCol(sk, dr)
				Vs.SetCol(sk, ds)
				Vt.SetCol(sk, dt)
				sk++
			}
		}
	}

	return Vr, Vs, Vt
}

// Dmatrices3D computes the differentiation matrices Dr, Ds, Dt
// Given the Vandermonde matrix V and points (r,s,t)
func Dmatrices3D(N int, r, s, t []float64, V utils.Matrix) (Dr, Ds, Dt utils.Matrix) {
	// Get gradient Vandermonde matrices
	Vr, Vs, Vt := GradVandermonde3D(N, r, s, t)

	// Compute V inverse
	Vinv := V.InverseWithCheck()

	// Dr = Vr * V^{-1}, etc.
	Dr = Vr.Mul(Vinv)
	Ds = Vs.Mul(Vinv)
	Dt = Vt.Mul(Vinv)

	return Dr, Ds, Dt
}
