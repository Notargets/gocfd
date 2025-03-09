package DG2D

import (
	"encoding/binary"
	"fmt"
	"os"

	"github.com/notargets/avs/geometry"
	"github.com/notargets/gocfd/geometry2D"
)

func CreateAVSGraphMesh(dfr *DFR2D) (gm geometry.TriMesh) {
	var (
		NpFlux       = dfr.FluxElement.Np
		NpInt        = dfr.FluxElement.NpInt
		NpEdge       = dfr.FluxElement.NpEdge
		TriNp        = 3 + 3*NpEdge
		EdgeR, EdgeS = make([]float64, TriNp), make([]float64, TriNp)
		VX, VY       = dfr.VX, dfr.VY
		FluxX, FluxY = dfr.FluxX, dfr.FluxY
	)
	// Compose the edges of the triangle,
	// they include the vertices and the edge nodes on the RT element,
	// in counter-clockwise progression
	VertR := []float64{-1, 1, -1}
	VertS := []float64{-1, -1, 1}
	var ii int
	for n := 0; n < 3; n++ {
		EdgeR[ii], EdgeS[ii] = VertR[n], VertS[n]
		ii++
		offset := 2*NpInt + n*NpEdge
		for i := 0; i < NpEdge; i++ {
			EdgeR[ii], EdgeS[ii] = dfr.FluxElement.R.DataP[offset+i],
				dfr.FluxElement.S.DataP[offset+i]
			ii++
		}
	}
	gmRS := geometry2D.TriangulateTriangle(EdgeR, EdgeS,
		dfr.FluxElement.R.DataP[:NpInt], dfr.FluxElement.S.DataP[:NpInt])
	// Output a 2025 AVS compatible TriMesh
	// The edges are first, including vertices, then the interior points
	// Each edge begins with a vertex, then edge points, in CCW order
	lenElement := NpInt + 3*NpEdge + 3
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
				ind := k*NpFlux + n*NpEdge + 2*NpInt + i
				XY[2*sk+0], XY[2*sk+1] =
					float32(FluxX.DataP[ind]), float32(FluxY.DataP[ind])
				sk++
			}
		}
		for i := 0; i < NpInt; i++ {
			ind := k*NpFlux + i
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

func WriteAVSGraphMesh(gm geometry.TriMesh, fileName string) {
	var (
		err         error
		file        *os.File
		lenXYCoords = len(gm.XY)
		lenTriVerts = len(gm.TriVerts)
	)
	file, err = os.Create(fileName)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	fmt.Printf("Number of Coordinate Pairs: %d\n", lenXYCoords)
	fmt.Printf("Number of RT Triangle Elements: %d\n", lenTriVerts/3)
	nDimensions := int64(2) // 2D
	binary.Write(file, binary.LittleEndian, nDimensions)
	binary.Write(file, binary.LittleEndian, lenTriVerts)
	binary.Write(file, binary.LittleEndian, gm.TriVerts)
	binary.Write(file, binary.LittleEndian, lenXYCoords)
	binary.Write(file, binary.LittleEndian, gm.XY)
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
