package DG2D

import (
	"encoding/binary"
	"fmt"
	"os"

	"github.com/notargets/avs/geometry"
	"github.com/notargets/gocfd/geometry2D"
)

func CreateAVSGraphMesh(dfr *DFR2D) (gm geometry.TriMesh) {
	// TODO: Fix this to:
	// TODO: triangulate a single element in (r,s) space
	// TODO: compose a TriMesh with the actual (x,y) coordinates using the
	// TODO: indices from the triangulation
	var (
		NpInt        = dfr.FluxElement.NpInt
		NpEdge       = dfr.FluxElement.NpEdge
		TriNp        = 3 + 3*NpEdge
		EdgeX, EdgeY = make([]float64, TriNp), make([]float64, TriNp)
	)
	// Compose the edges of the triangle,
	// they include the vertices and the edge nodes on the RT element,
	// in counter-clockwise progression
	var ii int
	for n := 0; n < 3; n++ {
		EdgeX[ii], EdgeY[ii] = dfr.VX.DataP[n], dfr.VY.DataP[n]
		ii++
		offset := 2*NpInt + n*NpEdge
		for i := 0; i < NpEdge; i++ {
			EdgeX[ii], EdgeY[ii] = dfr.FluxX.DataP[offset+i], dfr.FluxY.DataP[offset+i]
			ii++
		}
	}
	// Output a 2025 AVS compatible TriMesh
	// The edges are first, including vertices, then the interior points
	// Each edge begins with a vertex, then edge points, in CCW order
	gm = geometry2D.TriangulateTriangle(EdgeX, EdgeY,
		dfr.FluxX.DataP[:NpInt], dfr.FluxY.DataP[:NpInt])
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
