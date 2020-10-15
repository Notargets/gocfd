package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/utils"
)

type Triangulation struct {
	EToV  utils.Matrix         // K x 3 matrix mapping vertices to triangles
	Edges map[EdgeNumber]*Edge // map of edges, key is the edge number, an int packed with the two vertices of each edge
}

type EdgeNumber uint64

func NewEdgeNumber(verts [2]int) (packed EdgeNumber) {
	// This packs two index coordinates into two 32 bit unsigned integers to act as a hash and an indirect access method
	var (
		limit = math.MaxUint32
	)
	for _, vert := range verts {
		if vert < 0 || vert > limit {
			panic(fmt.Errorf("unable to pack two ints into a uint64, have %d and %d as inputs",
				verts[0], verts[1]))
		}
	}
	var i1, i2 int
	if verts[0] < verts[1] {
		i1, i2 = verts[0], verts[1]
	} else {
		i1, i2 = verts[1], verts[0]
	}
	packed = EdgeNumber(i1 + i2<<32)
	return
}

func (en EdgeNumber) GetVertices() (verts [2]int) {
	var (
		enTmp EdgeNumber
	)
	enTmp = en >> 32
	verts[1] = int(enTmp)
	verts[0] = int(en - enTmp*(1<<32))
	return
}

func NewTriangulation(EToV utils.Matrix) (tmesh *Triangulation) {
	tmesh = &Triangulation{
		EToV: EToV,
	}
	return
}

func (tmesh *Triangulation) NewEdge(verts [2]int, connectedElementNumber int) (e *Edge) {
	return
}

type Edge struct {
	ConnectedTris          []uint32 // Index number of triangles connected to this edge
	ConnectedTriDirection  []bool   // If false, the edge runs from smaller index to larger within the connected triangle
	ConnectedTriEdgeNumber []uint8  // For the connected triangles, what is the edge number (one of 0, 1 or 2)
}

func ConstructEdges() {

}
