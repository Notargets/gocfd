package gonudg

// INDEXING NOTE: Original C++ code uses 1-based indexing to emulate Matlab behavior.
// This Go port uses standard 0-based indexing. Example conversions:
//   C++: sk = 1; V3D(All,sk) = ...    ->    Go: sk = 0; V3D.SetCol(sk, ...)
//   C++: Fmask[1] (first face)        ->    Go: Fmask[0] (first face)
// The indexing has been correctly translated throughout this port.

import (
	"fmt"
	"github.com/notargets/gocfd/utils"
	"math"
)

// StartUp3D sets up the 3D discontinuous Galerkin operators
// This is the 0-based index version of the C++ StartUp3D function
type DG3D struct {
	// Polynomial order
	N int

	// Number of nodes
	Np      int     // Total nodes per element
	Nfp     int     // Nodes per face
	Nfaces  int     // Number of faces (4 for tetrahedron)
	NODETOL float64 // Node tolerance

	// Node coordinates in reference element
	r, s, t []float64

	// Physical coordinates
	X, Y, Z utils.Matrix // [Np][K] matrices

	// Vandermonde matrix and inverse
	V, Vinv utils.Matrix

	// Mass matrix
	MassMatrix utils.Matrix

	// Differentiation matrices
	Dr, Ds, Dt utils.Matrix

	// Weak differentiation matrices
	Drw, Dsw, Dtw utils.Matrix

	// Face mask
	Fmask [][]int // [Nfp][Nfaces]

	// Face coordinates
	Fx, Fy, Fz utils.Matrix

	// Lift matrix
	LIFT utils.Matrix

	// Geometric factors
	Rx, Ry, Rz utils.Matrix
	Sx, Sy, Sz utils.Matrix
	Tx, Ty, Tz utils.Matrix
	J          utils.Matrix

	// Surface geometric factors
	Nx, Ny, Nz utils.Matrix
	SJ         utils.Matrix
	Fscale     utils.Matrix

	// Connectivity
	EToV   [][]int // Element to vertex
	EToE   [][]int // Element to element
	EToF   [][]int // Element to face
	BCType []int   // Boundary condition types

	// Mesh info
	K          int       // Number of elements
	VX, VY, VZ []float64 // Vertex coordinates

	// Maps
	VmapM, VmapP []int // Volume to surface node maps
	MapM, MapP   []int // Surface node maps
	VmapB, MapB  []int // Boundary node maps
}

// NewDG3D creates and initializes a new 3D DG solver
func NewDG3D(N int, VX, VY, VZ []float64, EToV [][]int) (*DG3D, error) {
	dg := &DG3D{
		N:       N,
		VX:      VX,
		VY:      VY,
		VZ:      VZ,
		EToV:    EToV,
		K:       len(EToV),
		NODETOL: 1e-7,
		Nfaces:  4,
	}

	// Initialize
	if err := dg.StartUp3D(); err != nil {
		return nil, err
	}

	return dg, nil
}

// StartUp3D initializes the 3D DG operators
func (dg *DG3D) StartUp3D() error {
	// Definition of constants
	dg.Np = (dg.N + 1) * (dg.N + 2) * (dg.N + 3) / 6
	dg.Nfp = (dg.N + 1) * (dg.N + 2) / 2

	// Compute nodal set
	x1, y1, z1 := Nodes3D(dg.N)
	dg.r, dg.s, dg.t = XYZtoRST(x1, y1, z1)

	// Build reference element matrices
	dg.V = Vandermonde3D(dg.N, dg.r, dg.s, dg.t)
	dg.Vinv = dg.V.InverseWithCheck()
	dg.MassMatrix = dg.Vinv.Transpose().Mul(dg.Vinv)

	// Build differentiation matrices
	dg.Dr, dg.Ds, dg.Dt = Dmatrices3D(dg.N, dg.r, dg.s, dg.t, dg.V)

	// Build coordinates of all the nodes
	dg.X = utils.NewMatrix(dg.Np, dg.K)
	dg.Y = utils.NewMatrix(dg.Np, dg.K)
	dg.Z = utils.NewMatrix(dg.Np, dg.K)

	// Map from reference to physical elements
	for k := 0; k < dg.K; k++ {
		va := dg.EToV[k][0]
		vb := dg.EToV[k][1]
		vc := dg.EToV[k][2]
		vd := dg.EToV[k][3]

		for i := 0; i < dg.Np; i++ {
			dg.X.Set(i, k, 0.5*(-(1.0+dg.r[i]+dg.s[i]+dg.t[i])*dg.VX[va]+
				(1.0+dg.r[i])*dg.VX[vb]+
				(1.0+dg.s[i])*dg.VX[vc]+
				(1.0+dg.t[i])*dg.VX[vd]))

			dg.Y.Set(i, k, 0.5*(-(1.0+dg.r[i]+dg.s[i]+dg.t[i])*dg.VY[va]+
				(1.0+dg.r[i])*dg.VY[vb]+
				(1.0+dg.s[i])*dg.VY[vc]+
				(1.0+dg.t[i])*dg.VY[vd]))

			dg.Z.Set(i, k, 0.5*(-(1.0+dg.r[i]+dg.s[i]+dg.t[i])*dg.VZ[va]+
				(1.0+dg.r[i])*dg.VZ[vb]+
				(1.0+dg.s[i])*dg.VZ[vc]+
				(1.0+dg.t[i])*dg.VZ[vd]))
		}
	}

	// Find all the nodes that lie on each face
	dg.BuildFmask()

	// Extract face coordinates
	dg.ExtractFaceCoordinates()

	// Create surface integral terms
	if err := dg.Lift3D(); err != nil {
		return fmt.Errorf("Lift3D failed: %v", err)
	}

	// Calculate geometric factors and normals
	if err := dg.Normals3D(); err != nil {
		return fmt.Errorf("Normals3D failed: %v", err)
	}

	// Compute Fscale = SJ ./ J(Fmask,:)
	dg.ComputeFscale()

	// Build connectivity matrix
	dg.tiConnect3D()

	// Build connectivity maps
	dg.BuildMaps3D()

	// Compute weak operators
	dg.ComputeWeakOperators()

	return nil
}

// BuildFmask finds all nodes that lie on each face
func (dg *DG3D) BuildFmask() {
	dg.Fmask = make([][]int, 4)

	// Face 1: t = -1
	for i := 0; i < dg.Np; i++ {
		if math.Abs(1.0+dg.t[i]) < dg.NODETOL {
			dg.Fmask[0] = append(dg.Fmask[0], i)
		}
	}

	// Face 2: s = -1
	for i := 0; i < dg.Np; i++ {
		if math.Abs(1.0+dg.s[i]) < dg.NODETOL {
			dg.Fmask[1] = append(dg.Fmask[1], i)
		}
	}

	// Face 3: r+s+t = -1
	for i := 0; i < dg.Np; i++ {
		if math.Abs(1.0+dg.r[i]+dg.s[i]+dg.t[i]) < dg.NODETOL {
			dg.Fmask[2] = append(dg.Fmask[2], i)
		}
	}

	// Face 4: r = -1
	for i := 0; i < dg.Np; i++ {
		if math.Abs(1.0+dg.r[i]) < dg.NODETOL {
			dg.Fmask[3] = append(dg.Fmask[3], i)
		}
	}
}

// ExtractFaceCoordinates extracts coordinates at face nodes
func (dg *DG3D) ExtractFaceCoordinates() {
	// Create face coordinate matrices
	dg.Fx = utils.NewMatrix(dg.Nfp*4, dg.K)
	dg.Fy = utils.NewMatrix(dg.Nfp*4, dg.K)
	dg.Fz = utils.NewMatrix(dg.Nfp*4, dg.K)

	// Extract coordinates for each face
	for face := 0; face < 4; face++ {
		for k := 0; k < dg.K; k++ {
			for i, nodeIdx := range dg.Fmask[face] {
				row := face*dg.Nfp + i
				dg.Fx.Set(row, k, dg.X.At(nodeIdx, k))
				dg.Fy.Set(row, k, dg.Y.At(nodeIdx, k))
				dg.Fz.Set(row, k, dg.Z.At(nodeIdx, k))
			}
		}
	}
}

// ComputeFscale computes Fscale = SJ ./ J(Fmask,:)
func (dg *DG3D) ComputeFscale() {
	dg.Fscale = utils.NewMatrix(dg.Nfp*4, dg.K)

	for face := 0; face < 4; face++ {
		for k := 0; k < dg.K; k++ {
			for i, nodeIdx := range dg.Fmask[face] {
				row := face*dg.Nfp + i
				dg.Fscale.Set(row, k, dg.SJ.At(row, k)/dg.J.At(nodeIdx, k))
			}
		}
	}
}

// ComputeWeakOperators computes weak differentiation operators
func (dg *DG3D) ComputeWeakOperators() {
	// Compute V*V^T
	VVT := dg.V.Mul(dg.V.Transpose())
	VVTinv := VVT.InverseWithCheck()

	// Get gradient Vandermonde matrices
	Vr, Vs, Vt := GradVandermonde3D(dg.N, dg.r, dg.s, dg.t)

	// Drw = (V*Vr^T)/(V*V^T), etc.
	dg.Drw = dg.V.Mul(Vr.Transpose()).Mul(VVTinv)
	dg.Dsw = dg.V.Mul(Vs.Transpose()).Mul(VVTinv)
	dg.Dtw = dg.V.Mul(Vt.Transpose()).Mul(VVTinv)
}
