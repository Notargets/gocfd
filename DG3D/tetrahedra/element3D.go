package tetrahedra

import (
	"github.com/notargets/gocfd/utils"
	"math"
	"sort"
)

type Element3D struct {
	*TetBasis

	// Node information
	X, Y, Z utils.Vector // Physical coordinates

	// Surface integration
	Fx, Fy, Fz []utils.Matrix // Physical coordinates on each face

	// Connectivity
	VmapM  []int // Maps face nodes to volume nodes
	VmapP  []int // Maps face nodes to neighbor volume nodes
	MapM   []int // Maps face nodes to unique face nodes
	MapP   []int // Maps to neighbor face nodes
	BCType []int // Boundary condition types per face

	// Geometric factors
	GeomFactors     *GeometricFactors
	FaceGeomFactors *FaceGeometricFactors
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

func NewElement3D(order int) (el *Element3D) {
	el = &Element3D{
		TetBasis: NewTetBasis(order),
	}
	return
}

func (el *Element3D) GeometricFactors3D(x, y, z, Dr, Ds,
	Dt utils.Matrix) *GeometricFactors {
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

func (el *Element3D) Normals3D(geom *GeometricFactors, fmask [][]int,
	K int) *FaceGeometricFactors {
	Nfaces := 4
	Nfp := len(fmask[0])

	nx := utils.NewMatrix(Nfp*Nfaces, K)
	ny := utils.NewMatrix(Nfp*Nfaces, K)
	nz := utils.NewMatrix(Nfp*Nfaces, K)
	sJ := utils.NewMatrix(Nfp*Nfaces, K)

	// Face 1: t = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[0][i]
			nx.Set(i, k, -geom.Tx.At(vid, k))
			ny.Set(i, k, -geom.Ty.At(vid, k))
			nz.Set(i, k, -geom.Tz.At(vid, k))
			sJ.Set(i, k, geom.J.At(vid, k)*math.Sqrt(
				nx.At(i, k)*nx.At(i, k)+ny.At(i, k)*ny.At(i, k)+nz.At(i, k)*nz.At(i, k)))

			// Normalize
			nx.Set(i, k, nx.At(i, k)/sJ.At(i, k)*geom.J.At(vid, k))
			ny.Set(i, k, ny.At(i, k)/sJ.At(i, k)*geom.J.At(vid, k))
			nz.Set(i, k, nz.At(i, k)/sJ.At(i, k)*geom.J.At(vid, k))
		}
	}

	// Face 2: s = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[1][i]
			nx.Set(Nfp+i, k, -geom.Sx.At(vid, k))
			ny.Set(Nfp+i, k, -geom.Sy.At(vid, k))
			nz.Set(Nfp+i, k, -geom.Sz.At(vid, k))
			sJ.Set(Nfp+i, k, geom.J.At(vid, k)*math.Sqrt(
				nx.At(Nfp+i, k)*nx.At(Nfp+i, k)+ny.At(Nfp+i, k)*ny.At(Nfp+i, k)+nz.At(Nfp+i, k)*nz.At(Nfp+i, k)))

			// Normalize
			nx.Set(Nfp+i, k, nx.At(Nfp+i, k)/sJ.At(Nfp+i, k)*geom.J.At(vid, k))
			ny.Set(Nfp+i, k, ny.At(Nfp+i, k)/sJ.At(Nfp+i, k)*geom.J.At(vid, k))
			nz.Set(Nfp+i, k, nz.At(Nfp+i, k)/sJ.At(Nfp+i, k)*geom.J.At(vid, k))
		}
	}

	// Face 3: r+s+t = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[2][i]
			nx.Set(2*Nfp+i, k, geom.Rx.At(vid, k)+geom.Sx.At(vid, k)+geom.Tx.At(vid, k))
			ny.Set(2*Nfp+i, k, geom.Ry.At(vid, k)+geom.Sy.At(vid, k)+geom.Ty.At(vid, k))
			nz.Set(2*Nfp+i, k, geom.Rz.At(vid, k)+geom.Sz.At(vid, k)+geom.Tz.At(vid, k))
			sJ.Set(2*Nfp+i, k, geom.J.At(vid, k)*math.Sqrt(
				nx.At(2*Nfp+i, k)*nx.At(2*Nfp+i, k)+ny.At(2*Nfp+i, k)*ny.At(2*Nfp+i, k)+nz.At(2*Nfp+i, k)*nz.At(2*Nfp+i, k)))

			// Normalize
			nx.Set(2*Nfp+i, k, nx.At(2*Nfp+i, k)/sJ.At(2*Nfp+i, k)*geom.J.At(vid, k))
			ny.Set(2*Nfp+i, k, ny.At(2*Nfp+i, k)/sJ.At(2*Nfp+i, k)*geom.J.At(vid, k))
			nz.Set(2*Nfp+i, k, nz.At(2*Nfp+i, k)/sJ.At(2*Nfp+i, k)*geom.J.At(vid, k))
		}
	}

	// Face 4: r = -1
	for k := 0; k < K; k++ {
		for i := 0; i < Nfp; i++ {
			vid := fmask[3][i]
			nx.Set(3*Nfp+i, k, -geom.Rx.At(vid, k))
			ny.Set(3*Nfp+i, k, -geom.Ry.At(vid, k))
			nz.Set(3*Nfp+i, k, -geom.Rz.At(vid, k))
			sJ.Set(3*Nfp+i, k, geom.J.At(vid, k)*math.Sqrt(
				nx.At(3*Nfp+i, k)*nx.At(3*Nfp+i, k)+ny.At(3*Nfp+i, k)*ny.At(3*Nfp+i, k)+nz.At(3*Nfp+i, k)*nz.At(3*Nfp+i, k)))

			// Normalize
			nx.Set(3*Nfp+i, k, nx.At(3*Nfp+i, k)/sJ.At(3*Nfp+i, k)*geom.J.At(vid, k))
			ny.Set(3*Nfp+i, k, ny.At(3*Nfp+i, k)/sJ.At(3*Nfp+i, k)*geom.J.At(vid, k))
			nz.Set(3*Nfp+i, k, nz.At(3*Nfp+i, k)/sJ.At(3*Nfp+i, k)*geom.J.At(vid, k))
		}
	}

	// Compute Fscale
	Fscale := utils.NewMatrix(Nfp*Nfaces, K)
	for f := 0; f < Nfaces; f++ {
		for i := 0; i < Nfp; i++ {
			for k := 0; k < K; k++ {
				vid := fmask[f%4][i]
				Fscale.Set(f*Nfp+i, k, sJ.At(f*Nfp+i, k)/geom.J.At(vid, k))
			}
		}
	}

	return &FaceGeometricFactors{
		Nx: nx, Ny: ny, Nz: nz,
		SJ:     sJ,
		Fscale: Fscale,
	}
}

func BuildMaps3D(K int, Np int, Nfp int, fmask [][]int, x, y, z utils.Matrix,
	conn *ConnectivityArrays) (VmapM, VmapP, MapM, MapP []int) {

	Nfaces := 4

	// Build global face to node mapping
	nodeIds := make(map[[3]float64]int)
	NODETOL := 1e-7

	// Unique node identification
	Nnodes := 0
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

	// Volume to global node mapping
	vmapM := make([]int, Nfp*Nfaces*K)
	vmapP := make([]int, Nfp*Nfaces*K)
	mapM := make([]int, Nfp*Nfaces*K)
	mapP := make([]int, Nfp*Nfaces*K)

	// Initialize with identity mapping
	for i := range vmapM {
		vmapM[i] = i
		vmapP[i] = i
		mapM[i] = i
		mapP[i] = i
	}

	// Create face to face mapping
	for k1 := 0; k1 < K; k1++ {
		for f1 := 0; f1 < Nfaces; f1++ {
			k2 := conn.EToE[k1][f1]
			f2 := conn.EToF[k1][f1]

			if k2 < k1 || (k2 == k1 && f2 < f1) {
				continue // Only process each face pair once
			}

			// Find matching nodes
			for i := 0; i < Nfp; i++ {
				n1 := fmask[f1][i]
				key1 := [3]float64{
					math.Round(x.At(n1, k1)/NODETOL) * NODETOL,
					math.Round(y.At(n1, k1)/NODETOL) * NODETOL,
					math.Round(z.At(n1, k1)/NODETOL) * NODETOL,
				}

				for j := 0; j < Nfp; j++ {
					n2 := fmask[f2][j]
					key2 := [3]float64{
						math.Round(x.At(n2, k2)/NODETOL) * NODETOL,
						math.Round(y.At(n2, k2)/NODETOL) * NODETOL,
						math.Round(z.At(n2, k2)/NODETOL) * NODETOL,
					}

					if key1 == key2 {
						id1 := k1*Nfaces*Nfp + f1*Nfp + i
						id2 := k2*Nfaces*Nfp + f2*Nfp + j

						vmapM[id1] = k1*Np + n1
						vmapM[id2] = k2*Np + n2
						vmapP[id1] = k2*Np + n2
						vmapP[id2] = k1*Np + n1

						mapP[id1] = id2
						mapP[id2] = id1
						break
					}
				}
			}
		}
	}

	return vmapM, vmapP, mapM, mapP
}

func Connect3D(EToV [][]int) *ConnectivityArrays {
	K := len(EToV)
	Nfaces := 4

	// Build face to vertex mapping for tetrahedron
	vn := [][]int{
		{0, 1, 2}, // Face 1
		{0, 1, 3}, // Face 2
		{1, 2, 3}, // Face 3
		{0, 2, 3}, // Face 4
	}

	// Create face to element+face mapping
	faces := make(map[[3]int]struct{ elem, face int })

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
				// Found matching face
				// Update both elements
				if k != match.elem {
					// Different elements share this face
					// This is an internal face
				}
			} else {
				faces[key] = struct{ elem, face int }{k, f}
			}
		}
	}

	// Build connectivity arrays
	EToE := make([][]int, K)
	EToF := make([][]int, K)

	for k := 0; k < K; k++ {
		EToE[k] = make([]int, Nfaces)
		EToF[k] = make([]int, Nfaces)

		// Initialize with self-connectivity
		for f := 0; f < Nfaces; f++ {
			EToE[k][f] = k
			EToF[k][f] = f
		}

		// Find actual connections
		for f := 0; f < Nfaces; f++ {
			v := make([]int, 3)
			for i := 0; i < 3; i++ {
				v[i] = EToV[k][vn[f][i]]
			}
			sort.Ints(v)
			key := [3]int{v[0], v[1], v[2]}

			// Search all faces for match
			for fkey, fdata := range faces {
				if fkey == key && fdata.elem != k {
					EToE[k][f] = fdata.elem
					EToF[k][f] = fdata.face
					break
				}
			}
		}
	}

	return &ConnectivityArrays{
		EToE: EToE,
		EToF: EToF,
	}
}
