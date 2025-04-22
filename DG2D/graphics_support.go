package DG2D

import (
	"encoding/gob"
	"fmt"
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

type FieldMetadata struct {
	NumFields        int // How many fields are in the [][]float32
	FieldNames       []string
	SolutionMetadata map[string]interface{} // Fields like ReynoldsNumber,
	// gamma...
	GitVersion string
}

type SingleFieldMetadata struct {
	Iterations int
	Time       float32
	Count      int // Number of fields
	Length     int // of each field, for skipping / readahead
}

type AVSFieldWriter struct {
	Dimensions        [2]int // NpGraph, KMax - Points per element and #elements
	NpGraph, KMax     int
	FieldMeta         *FieldMetadata
	MeshMetadata      *MeshMetadata
	FileName          string
	IsMetaDataWritten bool
	file              *os.File
	encoder           *gob.Encoder
	fields            map[string][]float32
}

func (dfr *DFR2D) NewAVSFieldWriter(fmd *FieldMetadata,
	fileName string, gm geometry.TriMesh) (fw *AVSFieldWriter) {
	var (
		err error
	)
	fw = &AVSFieldWriter{
		FieldMeta:    fmd,
		MeshMetadata: dfr.GetMeshMetadata(),
		FileName:     fileName,
		NpGraph:      len(gm.XY) / 2,
		KMax:         len(gm.TriVerts),
		fields:       make(map[string][]float32),
	}
	fw.file, err = os.OpenFile(fileName, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		panic(err)
	}
	fw.encoder = gob.NewEncoder(fw.file)
	if !fw.IsMetaDataWritten {
		fw.WriteFieldsMeta()
		fw.IsMetaDataWritten = true
	}
	return
}

func (fw *AVSFieldWriter) saveField(name string, field []float64) {
	if len(field) != fw.NpGraph {
		err := fmt.Errorf("Dimensions mismatch between field and mesh\n"+
			"Field length is: %d, Mesh Data Length: %d\n",
			len(field), fw.NpGraph)
		panic(err)
	}
	if _, exists := fw.fields[name]; !exists {
		fw.fields[name] = make([]float32, fw.NpGraph)
	}
	for i, f := range field {
		fw.fields[name][i] = float32(f)
	}
}

func (fw *AVSFieldWriter) Save(sfmd *SingleFieldMetadata, fields map[string][]float64) {
	for name, fld := range fields {
		fw.saveField(name, fld)
	}
	fw.appendFields(sfmd)
}

func (fw *AVSFieldWriter) WriteFieldsMeta() {
	var (
		err error
	)
	// Encode mesh metadata first
	if err = fw.encoder.Encode(fw.MeshMetadata); err != nil {
		panic(err)
	}
	// Encode field metadata
	if err = fw.encoder.Encode(fw.FieldMeta); err != nil {
		panic(err)
	}
}

func (fw *AVSFieldWriter) appendFields(sfmd *SingleFieldMetadata) {
	var (
		err error
	)
	// Encode single field metadata
	if err = fw.encoder.Encode(sfmd); err != nil {
		panic(err)
	}
	// Encode field data
	if err = fw.encoder.Encode(fw.fields); err != nil {
		panic(err)
	}
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

func (dfr *DFR2D) TriangulateRTElement() (gm geometry.TriMesh) {
	var (
		NpEdge    = dfr.FluxElement.NpEdge
		TriEdgeNp = 3 * (1 + NpEdge)
	)
	R, S := dfr.GetRSForGraphMesh()
	gm = geometry2D.TriangulateTriangle(
		R.DataP[:TriEdgeNp], S.DataP[:TriEdgeNp],
		R.DataP[TriEdgeNp:], S.DataP[TriEdgeNp:])
	return
}

func (dfr *DFR2D) CreateAVSGraphMesh() (gm geometry.TriMesh) {
	var (
		NpInt        = dfr.FluxElement.NpInt
		NpEdge       = dfr.FluxElement.NpEdge
		VX, VY       = dfr.VX, dfr.VY
		FluxX, FluxY = dfr.FluxX, dfr.FluxY
	)
	gmRS := dfr.TriangulateRTElement()

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

func (dfr *DFR2D) GetMeshMetadata() (md *MeshMetadata) {
	var (
		gm          = dfr.GraphMesh
		lenXYCoords = len(gm.XY)
		lenTriVerts = len(gm.TriVerts)
	)
	md = &MeshMetadata{
		Description:     "Hybrid Lagrangian / Raviart Thomas elements",
		NDimensions:     2,
		Order:           dfr.N,
		NumBaseElements: dfr.K,
		NumPerElement:   lenTriVerts / dfr.K,
		NumElements:     lenTriVerts,
		LenXY:           lenXYCoords,
		GitVersion:      GitVersion, // Will be autofilled at build time
	}
	return
}

func (dfr *DFR2D) OutputMesh(fileName string, BCXY map[string][][]float32) {
	var (
		err error
	)
	md := dfr.GetMeshMetadata()
	err = WriteMesh(fileName, md, dfr.GraphMesh, BCXY)
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
