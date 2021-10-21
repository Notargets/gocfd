package DG2D

import (
	"fmt"
	"math"

	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/utils"
)

type Triangulation struct {
	EToV  utils.Matrix            // K x 3 matrix mapping vertices to triangles
	Edges map[types.EdgeKey]*Edge // map of edges, key is the edge number, an int packed with the two vertices of each edge
	EtoE  [][3]int                // For each element in [k], the connected element for each of the three edges [0,1,2]
}

func NewTriangulation(VX, VY utils.Vector, EToV utils.Matrix, BCEdges types.BCMAP) (tmesh *Triangulation) {
	var (
		K, _ = EToV.Dims()
	)
	tmesh = &Triangulation{
		EToV:  EToV,
		Edges: make(map[types.EdgeKey]*Edge),
		EtoE:  make([][3]int, K),
	}
	// Create edges map
	for k := 0; k < K; k++ {
		tri := EToV.Row(k).DataP
		verts := [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
		// Create / store the edges for this triangle
		tmesh.NewEdge(VX, VY, [2]int{verts[0], verts[1]}, k, First)
		tmesh.NewEdge(VX, VY, [2]int{verts[1], verts[2]}, k, Second)
		tmesh.NewEdge(VX, VY, [2]int{verts[2], verts[0]}, k, Third)
	}
	// Initialize EtoE with -1
	for k := 0; k < K; k++ {
		tmesh.EtoE[k] = [3]int{-1, -1, -1}
	}
	// Traverse edge map, filling in EtoE connections
	for _, edges := range tmesh.Edges {
		if edges.NumConnectedTris == 2 {
			en1, en2 := edges.ConnectedTriEdgeNumber[0], edges.ConnectedTriEdgeNumber[1]
			k1, k2 := edges.ConnectedTris[0], edges.ConnectedTris[1]
			tmesh.EtoE[k1][en1] = int(k2)
			tmesh.EtoE[k2][en2] = int(k1)
		}
	}
	// Insert BCs into edges map
	var err error
	var reversedLastEdge bool
	for key, edges := range BCEdges {
		flag := key.GetFLAG()
		switch flag {
		case types.BC_Far, types.BC_IVortex, types.BC_Wall, types.BC_In, types.BC_Out, types.BC_Cyl:
			for _, e := range edges {
				ee := tmesh.Edges[e.GetKey()]
				ee.BCType = flag
			}
		case types.BC_PeriodicReversed:
			fmt.Printf("reversing periodic edge\n")
			reversedLastEdge = true
			fallthrough
		case types.BC_Periodic:
			l := len(edges)
			if l%2 != 0 {
				err = fmt.Errorf("periodic boundaries must be present in pairs, have dimension %d, not even", l)
				panic(err)
			}
			l2 := l / 2
			// First, treat the edge lists as connected curves and order them
			var periodicEdges [2]types.Curve
			for i := 0; i < 2; i++ {
				periodicEdges[i] = make(types.Curve, l2)
			}
			for i := 0; i < l2; i++ {
				periodicEdges[0][i] = edges[i]
				periodicEdges[1][i] = edges[i+l2]
			}
			var reverseEdge bool
			for i := 0; i < 2; i++ {
				if reversedLastEdge && i == 1 { // Reverse one edge if requested
					reverseEdge = true
				}
				periodicEdges[i], _ = periodicEdges[i].ReOrder(reverseEdge)
			}
			for i := 0; i < l2; i++ {
				e1, e2 := periodicEdges[0][i], periodicEdges[1][i] // paired edges
				ee1, ee2 := tmesh.Edges[e1.GetKey()], tmesh.Edges[e2.GetKey()]
				ee1.BCType = flag
				ee2.BCType = flag
				k2 := int(ee2.ConnectedTris[0])
				en2 := ee2.ConnectedTriEdgeNumber[0]
				dir := ee2.ConnectedTriDirection[0] // Direction is relative to order of edge traversal in element/tri
				ee1.AddTri(e2.GetKey(), k2, 1, en2, dir, VX, VY)
			}
		default:
			err = fmt.Errorf("BC type %s not implemented yet", flag.String())
			panic(err)
		}
	}
	return
}

func (tmesh *Triangulation) GetTriVerts(k uint32) (verts [3]int) {
	tri := tmesh.EToV.Row(int(k)).DataP
	verts = [3]int{int(tri[0]), int(tri[1]), int(tri[2])}
	return
}

func (tmesh *Triangulation) NewEdge(VX, VY utils.Vector, verts [2]int, connectedElementNumber int, intEdgeNumber InternalEdgeNumber) (e *Edge) {
	var (
		ok bool
	)
	/*
		The input vertices are ordered as the normal traversal within the triangle
		The reversal boolean checks if this "natural order" is opposite to the "always least first"
		orientation of the edge key
	*/
	// Determine edge direction
	var dir InternalEdgeDirection
	if verts[0] > verts[1] {
		dir = Reversed
	}
	// Check if edge is already stored, allocate new one if not
	en := types.NewEdgeKey(verts)
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
	e.AddTri(en, connectedElementNumber, conn, intEdgeNumber, dir, VX, VY)
	return
}

func (e *Edge) AddTri(en types.EdgeKey, k, conn int, intEdgeNumber InternalEdgeNumber, direction InternalEdgeDirection, VX, VY utils.Vector) {
	e.ConnectedTris[conn] = uint32(k)
	e.ConnectedTriDirection[conn] = direction
	e.ConnectedTriEdgeNumber[conn] = intEdgeNumber
	e.NumConnectedTris++
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
	BCType                 types.BCFLAG             // If not connected to two tris, this field will be used
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

func (e *Edge) GetEdgeNormal(conn int, en types.EdgeKey, dfr *DFR2D) (normal [2]float64) {
	norm := func(vec [2]float64) (n float64) {
		n = math.Sqrt(vec[0]*vec[0] + vec[1]*vec[1])
		return
	}
	normalize := func(vec [2]float64) (normed [2]float64) {
		n := norm(vec)
		for i := 0; i < 2; i++ {
			normed[i] = vec[i] / n
		}
		return
	}
	revDir := bool(e.ConnectedTriDirection[conn])
	x1, x2 := GetEdgeCoordinates(en, revDir, dfr.VX, dfr.VY)
	dx := [2]float64{x2[0] - x1[0], x2[1] - x1[1]}
	normal = normalize([2]float64{-dx[1], dx[0]})
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

func GetEdgeCoordinates(en types.EdgeKey, rev bool, VX, VY utils.Vector) (x1, x2 [2]float64) {
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
