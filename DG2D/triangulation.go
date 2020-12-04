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

func NewTriangulation(VX, VY utils.Vector, EToV, BCType utils.Matrix) (tmesh *Triangulation) {
	tmesh = &Triangulation{
		EToV:  EToV,
		Edges: make(map[EdgeNumber]*Edge),
	}
	K, _ := EToV.Dims()
	// Create edges map
	for k := 0; k < K; k++ {
		tri := EToV.Row(k).Data()
		verts := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
		bcs := BCType.Row(k).Data()
		bcFaces := [3]int{int(bcs[0]), int(bcs[1]), int(bcs[2])}
		// Create / store the edges for this triangle
		tmesh.NewEdge(VX, VY, [2]int{verts[0], verts[1]}, k, First, bcFaces[0])
		tmesh.NewEdge(VX, VY, [2]int{verts[1], verts[2]}, k, Second, bcFaces[1])
		tmesh.NewEdge(VX, VY, [2]int{verts[2], verts[0]}, k, Third, bcFaces[2])
	}
	return
}

func (tmesh *Triangulation) GetTriVerts(k uint32) (verts [3]int) {
	tri := tmesh.EToV.Row(int(k)).Data()
	verts = [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
	return
}

func (tmesh *Triangulation) NewEdge(VX, VY utils.Vector,
	verts [2]int, connectedElementNumber int, intEdgeNumber InternalEdgeNumber, bcFace int) (e *Edge) {
	var (
		ok bool
	)
	/*
		The input vertices are ordered as the normal traversal within the triangle
	*/
	// Determine edge direction
	var dir InternalEdgeDirection
	if verts[0] > verts[1] {
		dir = Reversed
	}
	// Check if edge is already stored, allocate new one if not
	en := NewEdgeNumber(verts)
	conn := 1 // If edge exists, this will be the second (max) connection
	if e, ok = tmesh.Edges[en]; !ok {
		e = &Edge{}
		tmesh.Edges[en] = e
		conn = 0
	} else {
		// Check to ensure that we aren't adding more connections than are possible
		if e.NumConnectedTris > 1 {
			panic("incorrect edge construction, more than two connected triangles")
		}
	}
	e.AddTri(en, connectedElementNumber, conn, bcFace, intEdgeNumber, dir, VX, VY)
	return
}

func (e *Edge) AddTri(en EdgeNumber, k, conn, bcFace int,
	intEdgeNumber InternalEdgeNumber, direction InternalEdgeDirection,
	VX, VY utils.Vector) {
	e.ConnectedTris[conn] = uint32(k)
	e.ConnectedTriDirection[conn] = direction
	e.ConnectedTriEdgeNumber[conn] = intEdgeNumber
	e.NumConnectedTris++
	if bcFace != 0 {
		e.BCType = BCFLAG(bcFace)
	}
	// Calculate ||n|| scaling factor for each edge
	norm := func(vec [2]float64) (n float64) {
		n = math.Sqrt(vec[0]*vec[0] + vec[1]*vec[1])
		return
	}
	revDir := bool(e.ConnectedTriDirection[conn])
	edgeNumber := e.ConnectedTriEdgeNumber[conn]
	x1, x2 := GetEdgeCoordinates(en, revDir, VX, VY)
	dx := [2]float64{x2[0] - x1[0], x2[1] - x1[1]}
	edgeNorm := norm([2]float64{-dx[1], dx[0]}) // Norm of the edge normal
	// ||n|| = untransformed_edge_length / unit_tri_edge_length
	switch edgeNumber {
	case First, Third:
		e.IInII[conn] = edgeNorm / 2.
	case Second:
		e.IInII[conn] = edgeNorm / (2. * math.Sqrt(2))
	}
}

/*
Note that we do not use a slice for any of the fields inside of an edge - why?
	- each slice ([]type) has an overhead of pointer plus length and capacity ints, total of 24 bytes for an empty slice
	- adding to a slice using append() is expensive
	-> by using fixed allocations with the max of 2 ([2]type), we use less memory (total of 14 bytes) and compute
*/
type Edge struct {
	// Total Storage: 32 bytes (14 bytes, 64-bit aligned to 16, plus 2x8 bytes for ||n||)
	NumConnectedTris       uint8                    // Either 1 or 2
	ConnectedTris          [2]uint32                // Index numbers of triangles connected to this edge
	ConnectedTriDirection  [2]InternalEdgeDirection // If false(default), the edge runs from smaller to larger within the connected tri
	ConnectedTriEdgeNumber [2]InternalEdgeNumber    // For the connected triangles, what is the edge number (one of 0, 1 or 2)
	BCType                 BCFLAG                   // If not connected to two tris, this field will be used
	IInII                  [2]float64               // ||n|| scale factor to pre-multiply normal values prior to transforming into unit tri
}

func (e *Edge) GetEdgeLength() (len float64) {
	edgeNumber := e.ConnectedTriEdgeNumber[0]
	switch edgeNumber {
	case First, Third:
		len = 2. * e.IInII[0]
	case Second:
		len = (2. * math.Sqrt(2)) * e.IInII[0]
	}
	return
}

func (e *Edge) Print() (p string) {
	//for i, triNum := range e.ConnectedTris {
	for i := 0; i < int(e.NumConnectedTris); i++ {
		triNum := e.ConnectedTris[i]
		pp := fmt.Sprintf("Tri[%d] Edge[%d] BC=%s Reversed?%v,",
			triNum, e.ConnectedTriEdgeNumber[i], e.BCType.String(), e.ConnectedTriDirection[i])
		p += pp
	}
	return
}

func GetEdgeCoordinates(en EdgeNumber, rev bool, VX, VY utils.Vector) (x1, x2 [2]float64) {
	ev := en.GetVertices(!rev) // oriented for outward facing normals
	x1[0], x1[1] = VX.AtVec(ev[0]), VY.AtVec(ev[0])
	x2[0], x2[1] = VX.AtVec(ev[1]), VY.AtVec(ev[1])
	return
}

type InternalEdgeNumber uint8

const (
	First InternalEdgeNumber = iota
	Second
	Third
)

func (ien InternalEdgeNumber) Index() int {
	return int(ien)
}

func (ien InternalEdgeNumber) String() string {
	switch ien {
	case First:
		return "First"
	case Second:
		return "Second"
	case Third:
		return "Third"
	default:
		panic("unknown option")
	}
}

type InternalEdgeDirection bool

const (
	SmallestToLargest InternalEdgeDirection = false // Edge runs smallest vertex index to largest within triangle
	Reversed          InternalEdgeDirection = true
)

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

func (en EdgeNumber) GetVertices(rev bool) (verts [2]int) {
	var (
		enTmp EdgeNumber
	)
	enTmp = en >> 32
	verts[1] = int(enTmp)
	verts[0] = int(en - enTmp*(1<<32))
	if rev {
		verts[0], verts[1] = verts[1], verts[0]
	}
	return
}
