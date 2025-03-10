package DG2D

import (
	"encoding/binary"
	"encoding/gob"
	"os"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/avs/geometry"
	"github.com/notargets/gocfd/geometry2D"
)

var GitVersion string

// MeshMetadata is used to document mesh files
type MeshMetadata struct {
	Description     string // What type of method, e.g. Hybrid Lagrange/RT
	NDimensions     int    // Spatial dimensions
	Order           int    // Polynomial order
	NumElements     int    // Total number of elements in the mesh
	NumBaseElements int    // Number of elements, excluding the sub-element tris
	NumPerElement   int    // Elements are triangulated to approximate the poly
	LenXY           int    // Length of the XY coordinates in the mesh
	GitVersion      string
}

func (dfr *DFR2D) GetRSForGraphMesh() (R, S utils.Vector) {
	// Compose the edges of the triangle,
	// they include the vertices and the edge nodes on the RT element,
	// in counter-clockwise progression
	var (
		NpInt     = dfr.FluxElement.NpInt
		NpEdge    = dfr.FluxElement.NpEdge
		TriEdgeNp = 3 * (1 + NpEdge)
		VertR     = []float64{-1, 1, -1}
		VertS     = []float64{-1, -1, 1}
	)
	R, S = utils.NewVector(TriEdgeNp+NpInt), utils.NewVector(TriEdgeNp+NpInt)
	var ii int
	for n := 0; n < 3; n++ {
		R.DataP[ii], S.DataP[ii] = VertR[n], VertS[n]
		ii++
		offset := 2*NpInt + n*NpEdge
		for i := 0; i < NpEdge; i++ {
			R.DataP[ii], S.DataP[ii] = dfr.FluxElement.R.DataP[offset+i],
				dfr.FluxElement.S.DataP[offset+i]
			ii++
		}
	}
	for i := 0; i < NpInt; i++ {
		R.DataP[ii], S.DataP[ii] =
			dfr.FluxElement.R.DataP[i], dfr.FluxElement.S.DataP[i]
		ii++
	}
	return
}

func (dfr *DFR2D) CreateAVSGraphMesh() (gm geometry.TriMesh) {
	var (
		NpInt        = dfr.FluxElement.NpInt
		NpEdge       = dfr.FluxElement.NpEdge
		TriEdgeNp    = 3 * (1 + NpEdge)
		VX, VY       = dfr.VX, dfr.VY
		FluxX, FluxY = dfr.FluxX, dfr.FluxY
	)
	R, S := dfr.GetRSForGraphMesh()
	gmRS := geometry2D.TriangulateTriangle(
		R.DataP[:TriEdgeNp], S.DataP[:TriEdgeNp],
		R.DataP[TriEdgeNp:], S.DataP[TriEdgeNp:])

	// We need to convert the BCIndex, which lists vertices in VX,
	// VY into line segments

	// Output a 2025 AVS compatible TriMesh
	// The edges are first, including vertices, then the interior points
	// Each edge begins with a vertex, then edge points, in CCW order
	lenElement := 3*(NpEdge+1) + NpInt
	XY := make([]float32, 2*dfr.K*lenElement)
	Verts := make([][3]int64, dfr.K*len(gmRS.TriVerts))
	var sk int
	for k := 0; k < dfr.K; k++ {
		elVerts := dfr.Tris.GetTriVerts(uint32(k))
		// For each triangle, create a contiguous edge starting with the vertex
		for n := 0; n < 3; n++ {
			vi := elVerts[n]
			XY[2*sk+0], XY[2*sk+1] =
				float32(VX.DataP[vi]), float32(VY.DataP[vi])
			sk++
			for i := 0; i < NpEdge; i++ {
				// ind := k*NpFlux + 2*NpInt + n*NpEdge + i
				// ind := k + i*Kmax
				ind := k + (i+2*NpInt+n*NpEdge)*dfr.K
				XY[2*sk+0], XY[2*sk+1] =
					float32(FluxX.DataP[ind]), float32(FluxY.DataP[ind])
				sk++
			}
		}
		for i := 0; i < NpInt; i++ {
			// ind := k*NpFlux + i
			ind := k + i*dfr.K
			XY[2*sk+0], XY[2*sk+1] =
				float32(FluxX.DataP[ind]), float32(FluxY.DataP[ind])
			sk++
		}
		lenTriVerts := len(gmRS.TriVerts)
		for i := 0; i < lenTriVerts; i++ {
			for n := 0; n < 3; n++ {
				Verts[i+k*lenTriVerts][n] =
					gmRS.TriVerts[i][n] + int64(k*lenElement)
			}
		}
	}
	gm = geometry.NewTriMesh(XY, Verts)
	return
}

func WriteAVSSolutionField(field []float32, fileName string) {
	var (
		lenField = int64(len(field))
	)
	file, err := os.OpenFile(fileName, os.O_CREATE, os.ModeAppend)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// TODO: Amend this format to include the step number, maybe other meta info
	binary.Write(file, binary.LittleEndian, &lenField)
	binary.Write(file, binary.LittleEndian, &field)
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

func (dfr *DFR2D) OutputMesh(fileName string, BCXY map[string][][]float32) {
	var (
		err         error
		gm          = dfr.CreateAVSGraphMesh()
		lenXYCoords = len(gm.XY)
		lenTriVerts = len(gm.TriVerts)
	)
	md := &MeshMetadata{
		Description:     "Hybrid Lagrangian / Raviart Thomas elements",
		NDimensions:     2,
		Order:           dfr.N,
		NumBaseElements: dfr.K,
		NumPerElement:   lenTriVerts / dfr.K,
		NumElements:     lenTriVerts,
		LenXY:           lenXYCoords,
		GitVersion:      GitVersion, // Will be auto filled at build time
	}
	err = WriteMesh(fileName, md, gm, BCXY)
	if err != nil {
		panic(err)
	}
	return
}

// Function to write MeshMetadata and TriMesh sequentially
func WriteMesh(filename string, metadata *MeshMetadata,
	mesh geometry.TriMesh, BCXY map[string][][]float32) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	encoder := gob.NewEncoder(file)

	// Encode metadata first
	if err = encoder.Encode(metadata); err != nil {
		return err
	}

	// Encode mesh data
	if err = encoder.Encode(mesh); err != nil {
		return err
	}

	// Encode BCXY data
	if err = encoder.Encode(BCXY); err != nil {
		return err
	}

	return nil
}

// Function to read MeshMetadata and TriMesh sequentially
func ReadMesh(filename string) (md MeshMetadata, gm geometry.TriMesh,
	BCXY map[string][][]float32, err error) {
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	decoder := gob.NewDecoder(file)

	// Decode metadata
	if err = decoder.Decode(&md); err != nil {
		panic(err)
	}

	// Decode mesh data
	if err = decoder.Decode(&gm); err != nil {
		panic(err)
	}

	// Decode BCXY data
	if err = decoder.Decode(&BCXY); err != nil {
		panic(err)
	}
	return
}
