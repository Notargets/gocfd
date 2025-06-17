package tetrahedra

import (
	"fmt"
	"github.com/notargets/gocfd/DG3D/mesh"
	"github.com/notargets/gocfd/DG3D/mesh/readers"
	"github.com/notargets/gocfd/utils"
	"math"
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
	tetIndices []int   // Original indices of tet elements in full mesh
	*TetBasis          // Hesthaven tetrahedron PKD basis
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

// In extractTetElements(), it should verify ALL elements are tets:
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

// NewElement3DFromMesh creates an Element3D from an existing mesh
// This is the core constructor that does all the work
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

			// Direct access using the index - no GetNodeIndex needed!
			coords := el.Mesh.Vertices[nodeIdx]
			el.VX.Set(k*4+v, coords[0])
			el.VY.Set(k*4+v, coords[1])
			el.VZ.Set(k*4+v, coords[2])
		}
	}

	// Allocate physical coordinate matrices (Np x K)
	el.X = utils.NewMatrix(Np, K)
	el.Y = utils.NewMatrix(Np, K)
	el.Z = utils.NewMatrix(Np, K)

	// Get reference nodes from basis
	r, s, t := el.R, el.S, el.T

	// Transform reference nodes to physical coordinates for each element
	for k := 0; k < K; k++ {
		// Get vertices for this element
		v1x, v1y, v1z := el.VX.At(k*4+0), el.VY.At(k*4+0), el.VZ.At(k*4+0)
		v2x, v2y, v2z := el.VX.At(k*4+1), el.VY.At(k*4+1), el.VZ.At(k*4+1)
		v3x, v3y, v3z := el.VX.At(k*4+2), el.VY.At(k*4+2), el.VZ.At(k*4+2)
		v4x, v4y, v4z := el.VX.At(k*4+3), el.VY.At(k*4+3), el.VZ.At(k*4+3)

		// Transform each node using affine mapping
		for n := 0; n < Np; n++ {
			// Convert reference coordinates to barycentric
			// Reference tet vertices: v1=(-1,-1,-1), v2=(1,-1,-1), v3=(-1,1,-1), v4=(-1,-1,1)
			rn, sn, tn := r.At(n), s.At(n), t.At(n)

			// Barycentric coordinates (matching vertex ordering)
			L1 := -(1.0 + rn + sn + tn) / 2.0
			L2 := (1.0 + rn) / 2.0
			L3 := (1.0 + sn) / 2.0
			L4 := (1.0 + tn) / 2.0

			// Affine transformation to physical coordinates
			x := L1*v1x + L2*v2x + L3*v3x + L4*v4x
			y := L1*v1y + L2*v2y + L3*v3y + L4*v4y
			z := L1*v1z + L2*v2z + L3*v3z + L4*v4z

			// Store in matrix format (Np x K)
			el.X.Set(n, k, x)
			el.Y.Set(n, k, y)
			el.Z.Set(n, k, z)
		}
	}

	// Build face coordinates
	el.buildFaceCoordinates()
}

// buildFaceCoordinates extracts physical coordinates on each face
func (el *Element3D) buildFaceCoordinates() {
	K := el.K
	Nfp := el.Nfp // From TetBasis
	Nfaces := 4

	el.Fx = make([]utils.Matrix, Nfaces)
	el.Fy = make([]utils.Matrix, Nfaces)
	el.Fz = make([]utils.Matrix, Nfaces)

	for f := 0; f < Nfaces; f++ {
		el.Fx[f] = utils.NewMatrix(Nfp, K)
		el.Fy[f] = utils.NewMatrix(Nfp, K)
		el.Fz[f] = utils.NewMatrix(Nfp, K)

		for k := 0; k < K; k++ {
			for i := 0; i < Nfp; i++ {
				vid := el.Fmask[f][i] // Volume node index for this face node
				el.Fx[f].Set(i, k, el.X.At(vid, k))
				el.Fy[f].Set(i, k, el.Y.At(vid, k))
				el.Fz[f].Set(i, k, el.Z.At(vid, k))
			}
		}
	}
}

func (el *Element3D) Connect3D() *ConnectivityArrays {
	// fmt.Printf("DEBUG Connect3D: Starting\n")
	// fmt.Printf("  el.K = %d\n", el.K)
	// fmt.Printf("  el.Mesh = %v\n", el.Mesh)
	// if el.Mesh != nil {
	// 	fmt.Printf("  el.Mesh.NumElements = %d\n", el.Mesh.NumElements)
	// 	fmt.Printf("  el.Mesh.EToE = %v\n", el.Mesh.EToE)
	// 	fmt.Printf("  el.Mesh.EToF = %v\n", el.Mesh.EToF)
	// }
	// fmt.Printf("  el.tetIndices = %v\n", el.tetIndices)

	// Build connectivity if not available
	if el.Mesh.EToE == nil || el.Mesh.EToF == nil {
		// fmt.Printf("DEBUG: Building connectivity\n")
		el.Mesh.BuildConnectivity()
		// fmt.Printf("DEBUG: After BuildConnectivity:\n")
		// fmt.Printf("  el.Mesh.EToE = %v\n", el.Mesh.EToE)
		// fmt.Printf("  el.Mesh.EToF = %v\n", el.Mesh.EToF)
	}

	// Special case: if the mesh connectivity arrays are empty or nil after building,
	// create default connectivity for boundary elements
	if len(el.Mesh.EToE) == 0 {
		// Create connectivity for all tets as boundary elements
		el.EToE = make([][]int, el.K)
		el.EToF = make([][]int, el.K)

		for k := 0; k < el.K; k++ {
			el.EToE[k] = make([]int, 4) // 4 faces per tet
			el.EToF[k] = make([]int, 4)

			// Initialize all as boundary (no neighbors)
			for f := 0; f < 4; f++ {
				el.EToE[k][f] = -1
				el.EToF[k][f] = -1
			}
		}

		return &ConnectivityArrays{
			EToE: el.EToE,
			EToF: el.EToF,
		}
	}

	// If we filtered elements, we need to extract the relevant connectivity
	if len(el.tetIndices) > 0 && len(el.tetIndices) < el.Mesh.NumElements {
		// fmt.Printf("DEBUG: Filtering connectivity for %d tets out of %d elements\n",
		// 	len(el.tetIndices), el.Mesh.NumElements)
		// Create filtered connectivity arrays
		filteredEToE := make([][]int, el.K)
		filteredEToF := make([][]int, el.K)

		// Map from original indices to filtered indices
		indexMap := make(map[int]int)
		for newIdx, origIdx := range el.tetIndices {
			indexMap[origIdx] = newIdx
		}

		// Extract connectivity for tet elements only
		for newIdx, origIdx := range el.tetIndices {
			filteredEToE[newIdx] = make([]int, 4) // 4 faces per tet
			filteredEToF[newIdx] = make([]int, 4)

			// Check if original index is valid
			if origIdx >= len(el.Mesh.EToE) {
				// Element doesn't have connectivity - treat as boundary
				for f := 0; f < 4; f++ {
					filteredEToE[newIdx][f] = -1
					filteredEToF[newIdx][f] = -1
				}
				continue
			}

			for f := 0; f < 4; f++ {
				// Check if face connectivity exists
				if f >= len(el.Mesh.EToE[origIdx]) {
					// Face doesn't have connectivity - treat as boundary
					filteredEToE[newIdx][f] = -1
					filteredEToF[newIdx][f] = -1
					continue
				}

				neighbor := el.Mesh.EToE[origIdx][f]
				if neighbor == -1 {
					// Boundary face
					filteredEToE[newIdx][f] = -1
					filteredEToF[newIdx][f] = -1
				} else if mappedIdx, ok := indexMap[neighbor]; ok {
					// Neighbor is also a tet in our filtered set
					filteredEToE[newIdx][f] = mappedIdx
					filteredEToF[newIdx][f] = el.Mesh.EToF[origIdx][f]
				} else {
					// Neighbor is not a tet or not in filtered set - treat as boundary
					filteredEToE[newIdx][f] = -1
					filteredEToF[newIdx][f] = -1
				}
			}
		}

		el.ConnectivityArrays = &ConnectivityArrays{
			EToE: filteredEToE,
			EToF: filteredEToF,
		}
	} else {
		// fmt.Printf("DEBUG: Using connectivity directly (all elements are tets)\n")
		// All elements are tets, use connectivity directly
		// Initialize ConnectivityArrays if nil
		if el.ConnectivityArrays == nil {
			el.ConnectivityArrays = &ConnectivityArrays{}
		}
		el.ConnectivityArrays.EToE = el.Mesh.EToE
		el.ConnectivityArrays.EToF = el.Mesh.EToF
	}

	// fmt.Printf("DEBUG: Final connectivity:\n")
	// fmt.Printf("  el.ConnectivityArrays = %v\n", el.ConnectivityArrays)
	// if el.ConnectivityArrays != nil {
	// 	fmt.Printf("  el.ConnectivityArrays.EToE = %v\n", el.ConnectivityArrays.EToE)
	// 	fmt.Printf("  el.ConnectivityArrays.EToF = %v\n", el.ConnectivityArrays.EToF)
	// }

	return el.ConnectivityArrays
}

func (el *Element3D) GeometricFactors3D() (gf *GeometricFactors) {
	// Get physical coordinate matrices
	x, y, z := el.X, el.Y, el.Z

	// Get differentiation matrices
	Dr, Ds, Dt := el.Dr, el.Ds, el.Dt

	xr := Dr.Mul(x)
	xs := Ds.Mul(x)
	xt := Dt.Mul(x)
	yr := Dr.Mul(y)
	ys := Ds.Mul(y)
	yt := Dt.Mul(y)
	zr := Dr.Mul(z)
	zs := Ds.Mul(z)
	zt := Dt.Mul(z)

	// Compute Jacobian
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

	// Normalize all faces
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp*Nfaces; i++ {
			// Compute magnitude
			mag := math.Sqrt(nx.At(i, k)*nx.At(i, k) + ny.At(i, k)*ny.At(i, k) + nz.At(i, k)*nz.At(i, k))
			sJ.Set(i, k, mag)

			// Normalize to unit vector
			nx.Set(i, k, nx.At(i, k)/mag)
			ny.Set(i, k, ny.At(i, k)/mag)
			nz.Set(i, k, nz.At(i, k)/mag)
		}
	}

	// Scale sJ by volume Jacobian at face nodes
	for f := 0; f < Nfaces; f++ {
		for i := 0; i < Nfp; i++ {
			for k := 0; k < K; k++ {
				vid := fmask[f][i]
				idx := f*Nfp + i
				sJ.Set(idx, k, sJ.At(idx, k)*el.J.At(vid, k))
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
func (el *Element3D) BuildMaps3D() {
	var (
		K       = el.K
		Np      = el.Np
		Nfp     = el.Nfp
		fmask   = el.Fmask
		x, y, z = el.X, el.Y, el.Z
	)

	Nfaces := 4
	NF := Nfp * Nfaces // Total face points per element

	// Initialize arrays
	vmapM := make([]int, Nfp*Nfaces*K)
	vmapP := make([]int, Nfp*Nfaces*K)
	mapM := make([]int, Nfp*Nfaces*K)
	mapP := make([]int, Nfp*Nfaces*K)

	// Initialize mapM and mapP as identity (following C++: mapM.range(1,Nfp*Nfaces*K); mapP = mapM;)
	for i := 0; i < len(mapM); i++ {
		mapM[i] = i
		mapP[i] = i
	}

	// Find index of face nodes with respect to volume node ordering
	// (C++ loop: for k1=1; k1<=K; ++k1)
	for k1 := 0; k1 < K; k1++ {
		// Define target range in vmapM for element k1
		iL1 := k1 * NF

		// Map face nodes in element k1
		for f := 0; f < Nfaces; f++ {
			for i := 0; i < Nfp; i++ {
				idsL := iL1 + f*Nfp + i           // Index in vmapM array
				volumeNode := fmask[f][i] + k1*Np // Global volume node index
				vmapM[idsL] = volumeNode
			}
		}
	}

	// Initialize vmapP to vmapM (for boundary faces)
	copy(vmapP, vmapM)

	// Create face to face mapping
	NODETOL := 1e-7

	for k1 := 0; k1 < K; k1++ {
		for f1 := 0; f1 < Nfaces; f1++ {
			// find neighbor
			k2 := el.EToE[k1][f1]
			f2 := el.EToF[k1][f1]

			// Skip boundary faces
			// Handle both conventions: -1 or self-reference
			if k2 < 0 || (k2 == k1 && f2 == f1) {
				continue
			}

			// Also skip if k2 is out of range
			if k2 >= K {
				continue
			}

			skM := k1 * NF // offset to element k1
			skP := k2 * NF // offset to element k2

			// Define index ranges (C++ uses 1-based, we use 0-based)
			// idsM.range((f1-1)*Nfp+1+skM, f1*Nfp+skM)
			// idsP.range((f2-1)*Nfp+1+skP, f2*Nfp+skP)

			// Build index lists for this face pair
			// C++ uses 1-based indexing: idsM.range((f1-1)*Nfp+1+skM, f1*Nfp+skM)
			// Go uses 0-based, so no +1 needed
			idsM := make([]int, Nfp)
			idsP := make([]int, Nfp)
			for i := 0; i < Nfp; i++ {
				idsM[i] = f1*Nfp + i + skM
				idsP[i] = f2*Nfp + i + skP
			}

			// Find volume node numbers of left and right nodes
			// C++: vidM = vmapM(idsM); vidP = vmapM(idsP)
			vidM := make([]int, Nfp)
			vidP := make([]int, Nfp)
			for i := 0; i < Nfp; i++ {
				vidM[i] = vmapM[idsM[i]]
				vidP[i] = vmapM[idsP[i]]
			}

			// Extract coordinates using volume node numbers
			// C++: x1=x(vidM); y1=y(vidM); z1=z(vidM)
			x1 := make([]float64, Nfp)
			y1 := make([]float64, Nfp)
			z1 := make([]float64, Nfp)
			x2 := make([]float64, Nfp)
			y2 := make([]float64, Nfp)
			z2 := make([]float64, Nfp)

			for i := 0; i < Nfp; i++ {
				// M side
				nodeM := vidM[i] % Np
				elemM := vidM[i] / Np
				x1[i] = x.At(nodeM, elemM)
				y1[i] = y.At(nodeM, elemM)
				z1[i] = z.At(nodeM, elemM)

				// P side
				nodeP := vidP[i] % Np
				elemP := vidP[i] / Np
				x2[i] = x.At(nodeP, elemP)
				y2[i] = y.At(nodeP, elemP)
				z2[i] = z.At(nodeP, elemP)
			}

			// Compute distance matrix and find matches
			// C++: D = sqr(xM-trans(xP)) + sqr(yM-trans(yP)) + sqr(zM-trans(zP))
			for i := 0; i < Nfp; i++ {
				for j := 0; j < Nfp; j++ {
					dx := x1[i] - x2[j]
					dy := y1[i] - y2[j]
					dz := z1[i] - z2[j]
					dist := math.Sqrt(dx*dx + dy*dy + dz*dz)

					if dist < NODETOL {
						// C++: idM += (f1-1)*Nfp + skM; vmapP(idM) = vidP(idP)
						idM := idsM[i]
						vmapP[idM] = vidP[j]

						// C++: idP += (f2-1)*Nfp + skP; mapP(idM) = idP
						idP := idsP[j]
						mapP[idM] = idP
					}
				}
			}
		}
	}

	// Assign to element
	el.VmapM = vmapM
	el.VmapP = vmapP
	el.MapM = mapM
	el.MapP = mapP
}
func (el *Element3D) BuildMaps3D_Diagnostic() {
	var (
		K       = el.K
		Np      = el.Np
		Nfp     = el.Nfp
		fmask   = el.Fmask
		x, y, z = el.X, el.Y, el.Z
	)

	Nfaces := 4
	NF := Nfp * Nfaces // Total face points per element

	// Initialize arrays
	vmapM := make([]int, Nfp*Nfaces*K)
	vmapP := make([]int, Nfp*Nfaces*K)
	mapM := make([]int, Nfp*Nfaces*K)
	mapP := make([]int, Nfp*Nfaces*K)

	// Initialize mapM and mapP as identity
	for i := 0; i < len(mapM); i++ {
		mapM[i] = i
		mapP[i] = i
	}

	// Find index of face nodes with respect to volume node ordering
	for k1 := 0; k1 < K; k1++ {
		iL1 := k1 * NF
		for f := 0; f < Nfaces; f++ {
			for i := 0; i < Nfp; i++ {
				idsL := iL1 + f*Nfp + i
				volumeNode := fmask[f][i] + k1*Np
				vmapM[idsL] = volumeNode
			}
		}
	}

	// Initialize vmapP to vmapM (for boundary faces)
	copy(vmapP, vmapM)

	// Create face to face mapping
	NODETOL := 1e-7

	// Track statistics
	matchesFound := 0
	facesProcessed := 0
	minDist := 1e10
	maxMatchDist := 0.0

	// Debug first few elements
	debugElements := 3

	for k1 := 0; k1 < K && k1 < debugElements; k1++ {
		for f1 := 0; f1 < Nfaces; f1++ {
			k2 := el.EToE[k1][f1]
			f2 := el.EToF[k1][f1]

			if k2 < 0 || (k2 == k1 && f2 == f1) || k2 >= K {
				continue
			}

			facesProcessed++
			fmt.Printf("\nProcessing face: elem %d face %d -> elem %d face %d\n", k1, f1, k2, f2)

			skM := k1 * NF
			skP := k2 * NF

			// Build index lists
			idsM := make([]int, Nfp)
			idsP := make([]int, Nfp)
			for i := 0; i < Nfp; i++ {
				idsM[i] = f1*Nfp + i + skM
				idsP[i] = f2*Nfp + i + skP
			}

			// Get volume node numbers
			vidM := make([]int, Nfp)
			vidP := make([]int, Nfp)
			for i := 0; i < Nfp; i++ {
				vidM[i] = vmapM[idsM[i]]
				vidP[i] = vmapM[idsP[i]]
			}

			// Extract coordinates
			x1 := make([]float64, Nfp)
			y1 := make([]float64, Nfp)
			z1 := make([]float64, Nfp)
			x2 := make([]float64, Nfp)
			y2 := make([]float64, Nfp)
			z2 := make([]float64, Nfp)

			fmt.Printf("  M side coordinates:\n")
			for i := 0; i < Nfp; i++ {
				nodeM := vidM[i] % Np
				elemM := vidM[i] / Np
				x1[i] = x.At(nodeM, elemM)
				y1[i] = y.At(nodeM, elemM)
				z1[i] = z.At(nodeM, elemM)
				fmt.Printf("    Point %d: vidM=%d (elem %d, node %d) -> (%.6f, %.6f, %.6f)\n",
					i, vidM[i], elemM, nodeM, x1[i], y1[i], z1[i])
			}

			fmt.Printf("  P side coordinates:\n")
			for i := 0; i < Nfp; i++ {
				nodeP := vidP[i] % Np
				elemP := vidP[i] / Np
				x2[i] = x.At(nodeP, elemP)
				y2[i] = y.At(nodeP, elemP)
				z2[i] = z.At(nodeP, elemP)
				fmt.Printf("    Point %d: vidP=%d (elem %d, node %d) -> (%.6f, %.6f, %.6f)\n",
					i, vidP[i], elemP, nodeP, x2[i], y2[i], z2[i])
			}

			// Compute distance matrix
			fmt.Printf("  Distance matrix:\n")
			fmt.Printf("        ")
			for j := 0; j < Nfp; j++ {
				fmt.Printf("   P%d     ", j)
			}
			fmt.Printf("\n")

			faceMatches := 0
			for i := 0; i < Nfp; i++ {
				fmt.Printf("    M%d: ", i)
				for j := 0; j < Nfp; j++ {
					dx := x1[i] - x2[j]
					dy := y1[i] - y2[j]
					dz := z1[i] - z2[j]
					dist := math.Sqrt(dx*dx + dy*dy + dz*dz)

					if dist < minDist {
						minDist = dist
					}

					if dist < NODETOL {
						fmt.Printf(" MATCH   ")
						idM := idsM[i]
						vmapP[idM] = vidP[j]
						idP := idsP[j]
						mapP[idM] = idP
						matchesFound++
						faceMatches++
						if dist > maxMatchDist {
							maxMatchDist = dist
						}
					} else {
						fmt.Printf(" %.2e ", dist)
					}
				}
				fmt.Printf("\n")
			}
			fmt.Printf("  Matches found on this face: %d/%d\n", faceMatches, Nfp)
		}
	}

	fmt.Printf("\n=== SUMMARY ===\n")
	fmt.Printf("Faces processed: %d\n", facesProcessed)
	fmt.Printf("Total matches found: %d\n", matchesFound)
	fmt.Printf("Minimum distance seen: %.2e\n", minDist)
	fmt.Printf("Maximum match distance: %.2e\n", maxMatchDist)
	fmt.Printf("NODETOL: %.2e\n", NODETOL)

	// Assign to element
	el.VmapM = vmapM
	el.VmapP = vmapP
	el.MapM = mapM
	el.MapP = mapP
}
