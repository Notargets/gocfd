package tetrahedra

import (
	"fmt"
	"math"
	"sort"

	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/DG3D/mesh/readers"
	"github.com/notargets/gocfd/DG3D/tetrahedra/gonudg"
	"github.com/notargets/gocfd/utils"
)

type Element3D struct {
	K          int            // Number of elements
	VX, VY, VZ utils.Vector   // Vertex geometry, input from mesh, dimension K*4
	X, Y, Z    utils.Matrix   // Physical coordinates at each node (Np x K)
	Fx, Fy, Fz []utils.Matrix // Physical coordinates on each face
	// Connectivity
	EToV       [][]int // Element to Vertex map from mesh file
	EToP       []int   // Element to Partition map (optionally nil) [K]
	BCType     []int   // Boundary condition types per face from mesh file
	VmapM      []int   // Maps face nodes to volume nodes
	VmapP      []int   // Maps face nodes to neighbor volume nodes
	MapM       []int   // Maps face nodes to unique face nodes
	MapP       []int   // Maps to neighbor face nodes
	MapB       []int   // Boundary nodes map
	VmapB      []int   // Boundary nodes map
	tetIndices []int   // Original indices of tet elements in full mesh
	*TetBasis          // Basis functions and matrices
	*GeometricFactors
	*FaceGeometricFactors
	*ConnectivityArrays
	Mesh           *mesh.Mesh
	SplitElement3D []*Element3D
}

type ConnectivityArrays struct {
	EToE [][]int // Element to element connectivity
	EToF [][]int // Element to face connectivity
}

type GeometricFactors struct {
	X, Y, Z    utils.Matrix // Physical coordinates at each node
	Rx, Ry, Rz utils.Matrix // dr/dx, dr/dy, dr/dz
	Sx, Sy, Sz utils.Matrix // ds/dx, ds/dy, ds/dz
	Tx, Ty, Tz utils.Matrix // dt/dx, dt/dy, dt/dz
	J          utils.Matrix // Jacobian determinant
}

type FaceGeometricFactors struct {
	Nx, Ny, Nz utils.Matrix // Outward normal components
	SJ         utils.Matrix // Surface Jacobian (area scaling)
	Fscale     utils.Matrix // Face integration scaling = sJ/J(face)
}

// TetBasis represents the basis functions on reference tetrahedron
// Now uses gonudg functions internally
type TetBasis struct {
	N       int
	Np      int // Number of basis functions = (P+1)(P+2)(P+3)/6
	R, S, T utils.Vector
	V, VInv utils.Matrix // Vandermonde matrix and inverse
	// Mass and differentiation matrices
	M, MInv    utils.Matrix // Mass matrix, Inverse mass matrix
	Dr, Ds, Dt utils.Matrix // Differentiation matrices
	LIFT       utils.Matrix // Lift matrix for surface integrals
	Fmask      [][]int      // Face node indices
	Nfp        int          // Number of face points per face
}

// NewTetBasis creates a new tetrahedral basis using gonudg functions
func NewTetBasis(N int) (tb *TetBasis) {
	tb = &TetBasis{
		N:  N,
		Np: (N + 1) * (N + 2) * (N + 3) / 6,
	}

	// Get nodes in equilateral tetrahedron coordinates using gonudg
	XVec, YVec, ZVec := gonudg.Nodes3D(N)

	// Convert to reference tetrahedron coordinates using gonudg
	r, s, t := gonudg.XYZtoRST(XVec, YVec, ZVec)
	tb.R = utils.NewVector(len(r), r)
	tb.S = utils.NewVector(len(s), s)
	tb.T = utils.NewVector(len(t), t)

	// Build Vandermonde matrix using gonudg
	tb.V = gonudg.Vandermonde3D(N, r, s, t)
	tb.VInv = tb.V.InverseWithCheck()

	// Build differentiation matrices using gonudg
	tb.Dr, tb.Ds, tb.Dt = gonudg.Dmatrices3D(N, r, s, t, tb.V)

	// Build mass matrix
	tb.M = tb.VInv.Transpose().Mul(tb.VInv)
	tb.MInv = tb.M.InverseWithCheck()

	// Build face mask and lift matrix
	tb.buildFaceMask()
	tb.LIFT = gonudg.Lift3D(tb.N, tb.R.DataP, tb.S.DataP, tb.T.DataP, tb.V, tb.Fmask)

	return
}

// buildFaceMask finds nodes on each face
func (tb *TetBasis) buildFaceMask() {
	NODETOL := 1e-7
	tb.Fmask = make([][]int, 4)

	r := tb.R.DataP
	s := tb.S.DataP
	t := tb.T.DataP

	// Face 1: t = -1
	for i := 0; i < tb.Np; i++ {
		if math.Abs(1.0+t[i]) < NODETOL {
			tb.Fmask[0] = append(tb.Fmask[0], i)
		}
	}

	// Face 2: s = -1
	for i := 0; i < tb.Np; i++ {
		if math.Abs(1.0+s[i]) < NODETOL {
			tb.Fmask[1] = append(tb.Fmask[1], i)
		}
	}

	// Face 3: r+s+t = -1
	for i := 0; i < tb.Np; i++ {
		if math.Abs(1.0+r[i]+s[i]+t[i]) < NODETOL {
			tb.Fmask[2] = append(tb.Fmask[2], i)
		}
	}

	// Face 4: r = -1
	for i := 0; i < tb.Np; i++ {
		if math.Abs(1.0+r[i]) < NODETOL {
			tb.Fmask[3] = append(tb.Fmask[3], i)
		}
	}

	tb.Nfp = (tb.N + 1) * (tb.N + 2) / 2
}

// NewElement3D creates an Element3D from a mesh file
func NewElement3D(order int, meshFile string) (el *Element3D, err error) {
	// Read mesh file
	m, err := readers.ReadMeshFile(meshFile)
	if err != nil {
		return nil, err
	}

	// Use the common constructor
	return NewElement3DFromMesh(order, m)
}

// NewElement3DFromMesh creates an Element3D from an existing mesh
func NewElement3DFromMesh(order int, m *mesh.Mesh) (el *Element3D, err error) {
	el = &Element3D{
		TetBasis: NewTetBasis(order),
		Mesh:     m,
	}

	// Extract tetrahedral elements from mesh
	if err = el.extractTetElements(); err != nil {
		return nil, err
	}

	// Build physical coordinates from vertices and reference nodes
	el.buildPhysicalCoordinates()

	// Build connectivity and geometric factors
	el.ConnectivityArrays = el.Connect3D()
	el.GeometricFactors = el.GeometricFactors3D()
	el.FaceGeometricFactors = el.CalcFaceGeometry()

	el.BuildMaps3D()

	// Extract boundary conditions from mesh
	if err = el.extractBoundaryConditions(); err != nil {
		return nil, fmt.Errorf("failed to extract boundary conditions: %v", err)
	}

	// Extract periodic boundaries if present
	if err = el.extractPeriodicBoundaries(); err != nil {
		return nil, fmt.Errorf("failed to extract periodic boundaries: %v", err)
	}

	// Split mesh by partition if EToP is present
	if el.EToP != nil {
		if err = el.SplitByPartition(); err != nil {
			return nil, fmt.Errorf("failed to split mesh by partition: %v", err)
		}
	}

	return el, nil
}

// extractTetElements verifies all elements are tetrahedral
func (el *Element3D) extractTetElements() error {
	// Verify ALL elements are tetrahedral - no filtering!
	for i := 0; i < el.Mesh.NumElements; i++ {
		elemType := el.Mesh.ElementTypes[i]

		// Check if element is 3D
		if elemType.GetDimension() != 3 {
			return fmt.Errorf("element %d is not 3D (type=%v)", i, elemType)
		}

		// Check if element is tetrahedral
		if elemType != mesh.Tet && elemType != mesh.Tet10 {
			return fmt.Errorf("element %d is not tetrahedral (type=%v) - only tetrahedral meshes supported", i, elemType)
		}
	}

	// Store number of elements
	el.K = el.Mesh.NumElements

	// Build EToV connectivity - direct copy since all are tets
	el.EToV = make([][]int, el.K)
	for i := 0; i < el.K; i++ {
		// Get corner nodes only (first 4 nodes for tet)
		nodes := el.Mesh.EtoV[i]
		if len(nodes) >= 4 {
			el.EToV[i] = nodes[:4]
		} else {
			return fmt.Errorf("tetrahedral element %d has insufficient nodes", i)
		}
	}

	// Copy partition data if available
	if el.Mesh.EToP != nil {
		el.EToP = make([]int, el.K)
		copy(el.EToP, el.Mesh.EToP)
	}

	return nil
}

// buildPhysicalCoordinates builds the physical coordinates by transforming reference nodes
func (el *Element3D) buildPhysicalCoordinates() {
	// Get dimensions
	K := el.K
	Np := el.Np

	// Allocate vertex coordinates for all elements
	el.VX = utils.NewVector(K * 4) // 4 vertices per tet
	el.VY = utils.NewVector(K * 4)
	el.VZ = utils.NewVector(K * 4)

	// Extract vertex coordinates from mesh
	for k := 0; k < K; k++ {
		for v := 0; v < 4; v++ {
			nodeIdx := el.EToV[k][v] // This is already an array index, NOT a node ID!

			// Direct access using the index
			coord := el.Mesh.Vertices[nodeIdx]
			idx := k*4 + v
			el.VX.Set(idx, coord[0])
			el.VY.Set(idx, coord[1])
			el.VZ.Set(idx, coord[2])
		}
	}

	// Initialize physical coordinate matrices
	el.X = utils.NewMatrix(Np, K)
	el.Y = utils.NewMatrix(Np, K)
	el.Z = utils.NewMatrix(Np, K)

	// Reference nodes in [-1,1]^3
	r := el.R.DataP
	s := el.S.DataP
	t := el.T.DataP

	// Map reference to physical for each element
	for k := 0; k < K; k++ {
		// Get vertex indices for this element
		va := k*4 + 0
		vb := k*4 + 1
		vc := k*4 + 2
		vd := k*4 + 3

		// Transform each node
		for i := 0; i < Np; i++ {
			// Barycentric coordinates
			L1 := 0.5 * (1.0 + t[i])
			L2 := 0.5 * (1.0 + s[i])
			L3 := -0.5 * (1.0 + r[i] + s[i] + t[i])
			L4 := 0.5 * (1.0 + r[i])

			// Physical coordinates
			el.X.Set(i, k, L3*el.VX.At(va)+L4*el.VX.At(vb)+L2*el.VX.At(vc)+L1*el.VX.At(vd))
			el.Y.Set(i, k, L3*el.VY.At(va)+L4*el.VY.At(vb)+L2*el.VY.At(vc)+L1*el.VY.At(vd))
			el.Z.Set(i, k, L3*el.VZ.At(va)+L4*el.VZ.At(vb)+L2*el.VZ.At(vc)+L1*el.VZ.At(vd))
		}
	}

	// Build face coordinate arrays
	el.buildFaceCoordinates()
}

// buildFaceCoordinates extracts coordinates at face nodes
func (el *Element3D) buildFaceCoordinates() {
	Nfaces := 4
	el.Fx = make([]utils.Matrix, Nfaces)
	el.Fy = make([]utils.Matrix, Nfaces)
	el.Fz = make([]utils.Matrix, Nfaces)

	for face := 0; face < Nfaces; face++ {
		Nfp := len(el.Fmask[face])
		el.Fx[face] = utils.NewMatrix(Nfp, el.K)
		el.Fy[face] = utils.NewMatrix(Nfp, el.K)
		el.Fz[face] = utils.NewMatrix(Nfp, el.K)

		for k := 0; k < el.K; k++ {
			for i, nodeIdx := range el.Fmask[face] {
				el.Fx[face].Set(i, k, el.X.At(nodeIdx, k))
				el.Fy[face].Set(i, k, el.Y.At(nodeIdx, k))
				el.Fz[face].Set(i, k, el.Z.At(nodeIdx, k))
			}
		}
	}
}

// GeometricFactors3D computes metric terms using gonudg-based differentiation
func (el *Element3D) GeometricFactors3D() (gf *GeometricFactors) {
	// Get physical coordinate matrices
	x, y, z := el.X, el.Y, el.Z

	// Get differentiation matrices (already computed in TetBasis)
	Dr, Ds, Dt := el.Dr, el.Ds, el.Dt

	// Compute derivatives
	xr := Dr.Mul(x)
	xs := Ds.Mul(x)
	xt := Dt.Mul(x)
	yr := Dr.Mul(y)
	ys := Ds.Mul(y)
	yt := Dt.Mul(y)
	zr := Dr.Mul(z)
	zs := Ds.Mul(z)
	zt := Dt.Mul(z)

	// Compute Jacobian determinant
	J := utils.NewMatrix(xr.Rows(), xr.Cols())
	for i := 0; i < xr.Rows(); i++ {
		for j := 0; j < xr.Cols(); j++ {
			J.Set(i, j, xr.At(i, j)*(ys.At(i, j)*zt.At(i, j)-zs.At(i, j)*yt.At(i, j))-
				yr.At(i, j)*(xs.At(i, j)*zt.At(i, j)-zs.At(i, j)*xt.At(i, j))+
				zr.At(i, j)*(xs.At(i, j)*yt.At(i, j)-ys.At(i, j)*xt.At(i, j)))
		}
	}

	// Compute inverse metric terms
	rx := utils.NewMatrix(xr.Rows(), xr.Cols())
	ry := utils.NewMatrix(xr.Rows(), xr.Cols())
	rz := utils.NewMatrix(xr.Rows(), xr.Cols())
	sx := utils.NewMatrix(xr.Rows(), xr.Cols())
	sy := utils.NewMatrix(xr.Rows(), xr.Cols())
	sz := utils.NewMatrix(xr.Rows(), xr.Cols())
	tx := utils.NewMatrix(xr.Rows(), xr.Cols())
	ty := utils.NewMatrix(xr.Rows(), xr.Cols())
	tz := utils.NewMatrix(xr.Rows(), xr.Cols())

	for i := 0; i < xr.Rows(); i++ {
		for j := 0; j < xr.Cols(); j++ {
			rx.Set(i, j, (ys.At(i, j)*zt.At(i, j)-zs.At(i, j)*yt.At(i, j))/J.At(i, j))
			ry.Set(i, j, -(xs.At(i, j)*zt.At(i, j)-zs.At(i, j)*xt.At(i, j))/J.At(i, j))
			rz.Set(i, j, (xs.At(i, j)*yt.At(i, j)-ys.At(i, j)*xt.At(i, j))/J.At(i, j))

			sx.Set(i, j, -(yr.At(i, j)*zt.At(i, j)-zr.At(i, j)*yt.At(i, j))/J.At(i, j))
			sy.Set(i, j, (xr.At(i, j)*zt.At(i, j)-zr.At(i, j)*xt.At(i, j))/J.At(i, j))
			sz.Set(i, j, -(xr.At(i, j)*yt.At(i, j)-yr.At(i, j)*xt.At(i, j))/J.At(i, j))

			tx.Set(i, j, (yr.At(i, j)*zs.At(i, j)-zr.At(i, j)*ys.At(i, j))/J.At(i, j))
			ty.Set(i, j, -(xr.At(i, j)*zs.At(i, j)-zr.At(i, j)*xs.At(i, j))/J.At(i, j))
			tz.Set(i, j, (xr.At(i, j)*ys.At(i, j)-yr.At(i, j)*xs.At(i, j))/J.At(i, j))
		}
	}

	return &GeometricFactors{
		X: x, Y: y, Z: z,
		Rx: rx, Ry: ry, Rz: rz,
		Sx: sx, Sy: sy, Sz: sz,
		Tx: tx, Ty: ty, Tz: tz,
		J: J,
	}
}

// CalcFaceGeometry computes face normals and surface Jacobians
func (el *Element3D) CalcFaceGeometry() *FaceGeometricFactors {
	fmask := el.Fmask
	K := el.K
	Nfp := el.Nfp
	Nfaces := 4

	nx := utils.NewMatrix(Nfp*Nfaces, K)
	ny := utils.NewMatrix(Nfp*Nfaces, K)
	nz := utils.NewMatrix(Nfp*Nfaces, K)
	sJ := utils.NewMatrix(Nfp*Nfaces, K)

	// Face 1: t = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[0][i]
			nx.Set(i, k, -el.Tx.At(vid, k))
			ny.Set(i, k, -el.Ty.At(vid, k))
			nz.Set(i, k, -el.Tz.At(vid, k))
		}
	}

	// Face 2: s = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[1][i]
			nx.Set(Nfp+i, k, -el.Sx.At(vid, k))
			ny.Set(Nfp+i, k, -el.Sy.At(vid, k))
			nz.Set(Nfp+i, k, -el.Sz.At(vid, k))
		}
	}

	// Face 3: r+s+t = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[2][i]
			nx.Set(2*Nfp+i, k, el.Rx.At(vid, k)+el.Sx.At(vid, k)+el.Tx.At(vid, k))
			ny.Set(2*Nfp+i, k, el.Ry.At(vid, k)+el.Sy.At(vid, k)+el.Ty.At(vid, k))
			nz.Set(2*Nfp+i, k, el.Rz.At(vid, k)+el.Sz.At(vid, k)+el.Tz.At(vid, k))
		}
	}

	// Face 4: r = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[3][i]
			nx.Set(3*Nfp+i, k, -el.Rx.At(vid, k))
			ny.Set(3*Nfp+i, k, -el.Ry.At(vid, k))
			nz.Set(3*Nfp+i, k, -el.Rz.At(vid, k))
		}
	}

	// Normalize and compute surface Jacobian
	for i := 0; i < Nfp*Nfaces; i++ {
		for k := 0; k < K; k++ {
			sJ.Set(i, k, math.Sqrt(nx.At(i, k)*nx.At(i, k)+
				ny.At(i, k)*ny.At(i, k)+
				nz.At(i, k)*nz.At(i, k)))

			if sJ.At(i, k) > 0 {
				nx.Set(i, k, nx.At(i, k)/sJ.At(i, k))
				ny.Set(i, k, ny.At(i, k)/sJ.At(i, k))
				nz.Set(i, k, nz.At(i, k)/sJ.At(i, k))
			}
		}
	}

	// Compute Fscale
	Fscale := utils.NewMatrix(Nfp*Nfaces, K)
	for f := 0; f < Nfaces; f++ {
		for i := 0; i < Nfp; i++ {
			for k := 0; k < K; k++ {
				vid := fmask[f][i]
				Fscale.Set(f*Nfp+i, k, sJ.At(f*Nfp+i, k)/el.J.At(vid, k))
			}
		}
	}

	return &FaceGeometricFactors{
		Nx: nx, Ny: ny, Nz: nz,
		SJ:     sJ,
		Fscale: Fscale,
	}
}

// BuildMaps3D builds connectivity maps for DG
func (el *Element3D) BuildMaps3D() {
	// Use the existing implementation that properly computes face mappings
	K := el.K
	Np := el.Np
	Nfp := el.Nfp
	Nfaces := 4
	x, y, z := el.X, el.Y, el.Z

	// Build node-to-node mappings
	NODETOL := 1e-7
	nodeIds := make(map[[3]float64]int)
	Nnodes := 0

	// Assign unique IDs to nodes
	for k := 0; k < K; k++ {
		for n := 0; n < Np; n++ {
			key := [3]float64{
				math.Round(x.At(n, k)/NODETOL) * NODETOL,
				math.Round(y.At(n, k)/NODETOL) * NODETOL,
				math.Round(z.At(n, k)/NODETOL) * NODETOL,
			}
			if _, exists := nodeIds[key]; !exists {
				nodeIds[key] = Nnodes
				Nnodes++
			}
		}
	}

	// Initialize mappings
	el.VmapM = make([]int, Nfp*Nfaces*K)
	el.VmapP = make([]int, Nfp*Nfaces*K)
	el.MapM = make([]int, Nfp*Nfaces*K)
	el.MapP = make([]int, Nfp*Nfaces*K)

	// Build VmapM - maps face nodes to volume nodes
	for k := 0; k < K; k++ {
		for f := 0; f < Nfaces; f++ {
			for i := 0; i < Nfp; i++ {
				faceIdx := k*Nfaces*Nfp + f*Nfp + i
				volNode := el.Fmask[f][i] + k*Np
				el.VmapM[faceIdx] = volNode
			}
		}
	}

	// Initialize VmapP to VmapM (for boundaries)
	copy(el.VmapP, el.VmapM)

	// Initialize MapM and MapP as identity
	for i := 0; i < len(el.MapM); i++ {
		el.MapM[i] = i
		el.MapP[i] = i
	}

	// Build face-to-face mappings for interior faces
	for k1 := 0; k1 < K; k1++ {
		for f1 := 0; f1 < Nfaces; f1++ {
			k2 := el.EToE[k1][f1]
			f2 := el.EToF[k1][f1]

			// Skip boundary faces
			if k2 == k1 && f2 == f1 {
				continue
			}

			// Skip invalid neighbors
			if k2 < 0 || k2 >= K {
				continue
			}

			// Match nodes between faces
			for i := 0; i < Nfp; i++ {
				// Get node on face f1 of element k1
				n1 := el.Fmask[f1][i]
				key1 := [3]float64{
					math.Round(x.At(n1, k1)/NODETOL) * NODETOL,
					math.Round(y.At(n1, k1)/NODETOL) * NODETOL,
					math.Round(z.At(n1, k1)/NODETOL) * NODETOL,
				}

				// Find matching node on face f2 of element k2
				for j := 0; j < Nfp; j++ {
					n2 := el.Fmask[f2][j]
					key2 := [3]float64{
						math.Round(x.At(n2, k2)/NODETOL) * NODETOL,
						math.Round(y.At(n2, k2)/NODETOL) * NODETOL,
						math.Round(z.At(n2, k2)/NODETOL) * NODETOL,
					}

					if key1 == key2 {
						// Found matching nodes
						id1 := k1*Nfaces*Nfp + f1*Nfp + i
						id2 := k2*Nfaces*Nfp + f2*Nfp + j

						// Set VmapP
						el.VmapP[id1] = k2*Np + n2
						el.VmapP[id2] = k1*Np + n1

						// Set MapP
						el.MapP[id1] = id2
						el.MapP[id2] = id1
						break
					}
				}
			}
		}
	}

	// Build boundary mappings
	el.MapB = []int{}
	el.VmapB = []int{}

	for i := 0; i < len(el.VmapP); i++ {
		if el.VmapP[i] == el.VmapM[i] {
			el.MapB = append(el.MapB, i)
			el.VmapB = append(el.VmapB, el.VmapM[i])
		}
	}
}

// Connect3D builds element connectivity arrays
func (el *Element3D) Connect3D() (ca *ConnectivityArrays) {
	EToV := el.EToV
	K := el.K
	Nfaces := 4

	// Build face to vertex mapping for tetrahedron
	vn := [][]int{
		{0, 1, 2}, // Face 0
		{0, 1, 3}, // Face 1
		{1, 2, 3}, // Face 2
		{0, 2, 3}, // Face 3
	}

	// Initialize connectivity arrays
	EToE := make([][]int, K)
	EToF := make([][]int, K)

	for k := 0; k < K; k++ {
		EToE[k] = make([]int, Nfaces)
		EToF[k] = make([]int, Nfaces)

		// Initialize with self-connectivity (boundary)
		for f := 0; f < Nfaces; f++ {
			EToE[k][f] = k
			EToF[k][f] = f
		}
	}

	// Build face map and connectivity simultaneously
	type faceInfo struct {
		elem int
		face int
	}
	faces := make(map[[3]int]faceInfo)

	for k := 0; k < K; k++ {
		for f := 0; f < Nfaces; f++ {
			// Get vertices of this face
			v := make([]int, 3)
			for i := 0; i < 3; i++ {
				v[i] = EToV[k][vn[f][i]]
			}

			// Sort vertices to create unique key
			sort.Ints(v)
			key := [3]int{v[0], v[1], v[2]}

			if match, exists := faces[key]; exists {
				// Found matching face - set RECIPROCAL connectivity
				// Element k, face f connects to match.elem, face match.face
				EToE[k][f] = match.elem
				EToF[k][f] = match.face

				// CRITICAL FIX: Set reverse connectivity
				// match.elem, face match.face connects back to element k, face f
				EToE[match.elem][match.face] = k
				EToF[match.elem][match.face] = f
			} else {
				// First time seeing this face
				faces[key] = faceInfo{elem: k, face: f}
			}
		}
	}

	ca = &ConnectivityArrays{
		EToE: EToE,
		EToF: EToF,
	}
	return
}
