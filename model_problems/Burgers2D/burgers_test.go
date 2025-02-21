package Burgers2D

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"strings"
	"testing"
	"time"

	"github.com/notargets/avs/chart2d"
	avsUtils "github.com/notargets/avs/utils"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/avs/geometry"

	"github.com/notargets/gocfd/DG2D"
)

func TestNewBurgers2D(t *testing.T) {
	P := 2
	b2d := &Burgers2D{}
	mesh := "../../DG2D/test_data/test_tris_9.neu"
	b2d.DFR = DG2D.NewDFR2D(P, true, mesh)
	if testing.Verbose() {
		plotMesh(b2d.DFR, 30*time.Second)
	}
	return
}

func plotMesh(dfr *DG2D.DFR2D, waitTime time.Duration) {
	var xmin, ymin, xmax, ymax float64
	xmin, xmax = math.MaxFloat64, -math.MaxFloat64
	ymin, ymax = math.MaxFloat64, -math.MaxFloat64

	pMesh := makePlotMesh(dfr.VX, dfr.VY, dfr.Tris)

	for i, x := range dfr.VX.DataP {
		y := dfr.VY.AtVec(i)
		if x < xmin {
			xmin = x
		}
		if x > xmax {
			xmax = x
		}
		if y < ymin {
			ymin = y
		}
		if y > ymax {
			ymax = y
		}
	}
	xMin, xMax, yMin, yMax := getSquareBoundingBox(float32(xmin),
		float32(xmax), float32(ymin), float32(ymax))
	cc := chart2d.NewChart2D(xMin, xMax, yMin, yMax, 1920, 1080,
		avsUtils.WHITE, avsUtils.BLACK, 0.9)
	// cc.NewWindow("Geom", 0.9, screen.AUTO)
	cc.AddTriMesh(pMesh)
	time.Sleep(waitTime)
}

func makePlotMesh(VX, VY utils.Vector, Tris *DG2D.Triangulation) (tMesh geometry.TriMesh) {
	var (
		K = len(Tris.EtoE)
	)
	tMesh = geometry.TriMesh{
		XY:       make([]float32, 2*VX.Len()),
		TriVerts: make([][3]int64, K),
	}
	for i, x := range VX.DataP {
		y := VY.DataP[i]
		tMesh.XY[2*i] = float32(x)
		tMesh.XY[2*i+1] = float32(y)
	}
	for k := 0; k < K; k++ {
		for n := 0; n < 3; n++ {
			tMesh.TriVerts[k][n] = int64(Tris.EToV.At(k, n))
		}
	}
	return
}

func getSurfaceLines(XY []float32, edges []*geometry.EdgeGroup) (
	XYSurf []float32) {
	var x1, y1, x2, y2 float32
	for _, edgeGroup := range edges {
		// fmt.Printf("BC Group Name: [%s]\n", edgeGroup.GroupName)
		// if edgeGroup.GroupName == "wall" {
		for _, edgeXY := range edgeGroup.EdgeXYs {
			x1, y1 = edgeXY[0], edgeXY[1]
			x2, y2 = edgeXY[2], edgeXY[3]
			XYSurf = append(XYSurf, x1, y1, x2, y2)
		}
		// }
	}
	return
}

func getSquareBoundingBox(xMin, xMax, yMin, yMax float32) (xBMin,
	xBMax, yBMin, yBMax float32) {
	xRange := xMax - xMin
	yRange := yMax - yMin
	if yRange > xRange {
		yBMin = yMin
		yBMax = yMax
		xCent := xRange/2. + xMin
		xBMin = xCent - yRange/2.
		xBMax = xCent + yRange/2.
	} else {
		xBMin = xMin
		xBMax = xMax
		yCent := yRange/2. + yMin
		yBMin = yCent - xRange/2.
		yBMax = yCent + xRange/2.
	}
	return
}

func ReadGoCFDMesh(path string, verbose bool) (tMesh geometry.TriMesh,
	BCEdges []*geometry.EdgeGroup) {
	var (
		file        *os.File
		lenTriVerts int64
		lenXYCoords int64
	)

	// Read triVerts array
	triVerts := make([]int64, lenTriVerts)
	binary.Read(file, binary.LittleEndian, &triVerts)

	// Read XY coordinates
	binary.Read(file, binary.LittleEndian, &lenXYCoords)
	fmt.Printf("Number of Coordinates: %d\n", lenXYCoords)
	xy := make([]float64, lenXYCoords*2)
	binary.Read(file, binary.LittleEndian, &xy)

	xy32 := make([]float32, lenXYCoords*2)
	for i := range xy {
		xy32[i] = float32(xy[i])
	}

	verts := make([][3]int64, lenTriVerts/3)
	for i := range verts {
		verts[i][0] = triVerts[3*i]
		verts[i][1] = triVerts[3*i+1]
		verts[i][2] = triVerts[3*i+2]
	}
	tMesh = *geometry.NewTriMesh(xy32, verts)

	var nBCs int64
	binary.Read(file, binary.LittleEndian, &nBCs)
	BCEdges = make([]*geometry.EdgeGroup, nBCs)
	fmt.Printf("Number of BCs: %d\n", nBCs)
	// Read boundary conditions
	for n := 0; n < int(nBCs); n++ {
		var fString [16]byte
		binary.Read(file, binary.LittleEndian, &fString)
		bcName := strings.TrimRight(string(fString[:]), "\x00 ")
		var bcLen int64
		binary.Read(file, binary.LittleEndian, &bcLen)
		BCEdges[n] = geometry.NewEdgeGroup(bcName, int(bcLen))
		for i, _ := range BCEdges[n].EdgeXYs {
			binary.Read(file, binary.LittleEndian, &BCEdges[n].EdgeXYs[i])
		}
	}
	return
}
