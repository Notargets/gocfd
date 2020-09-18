package geometry2D

import (
	"fmt"
	"math"

	graphics2D "github.com/notargets/avs/geometry"
)

type Point struct {
	X [2]float64
}

type Edge struct {
	Tris        map[*Tri]struct{} // Associated triangles
	Verts       [2]int
	IsImmovable bool // Is a BC edge
}

func NewEdge(verts [2]int, isImmovableO ...bool) (e *Edge) {
	e = &Edge{
		Verts: verts,
	}
	if len(isImmovableO) != 0 {
		e.IsImmovable = isImmovableO[0]
	}
	return
}

func getUniqueVertices(edges ...*Edge) (verts []int) {
	vMap := make(map[int]struct{})
	for _, e := range edges {
		for _, v := range e.Verts {
			vMap[v] = struct{}{}
		}
	}
	for key, _ := range vMap {
		verts = append(verts, key)
	}
	return
}

func NewEdgeFromEdges(edges []*Edge, isImmovableO ...bool) (ee *Edge) {
	// Create a new edge by connecting two other edges
	vMap := make(map[int]int)
	for _, e := range edges {
		for _, v := range e.Verts {
			vMap[v]++
		}
	}
	var verts [2]int
	var ii int
	for key, val := range vMap {
		if val == 1 {
			verts[ii] = key
			ii++
		}
	}
	ee = &Edge{
		Verts: verts,
	}
	if len(isImmovableO) != 0 {
		ee.IsImmovable = isImmovableO[0]
	}
	return
}

func (e *Edge) DeleteTri(tri *Tri) {
	delete(e.Tris, tri)
}

func (e *Edge) AddTri(tri *Tri) {
	if e.Tris == nil {
		e.Tris = make(map[*Tri]struct{})
	}
	e.Tris[tri] = struct{}{}
}

func (e *Edge) GetOpposingVertices() (pts []int) {
	// Get vertices opposing this edge from all tris with this edge, could be 1, 2 or 0 points (unconnected edge)
	ptMap := make(map[int]struct{})
	ptMap[e.Verts[0]] = struct{}{}
	ptMap[e.Verts[1]] = struct{}{}
	for tri := range e.Tris {
		for _, ee := range tri.Edges {
			for _, ptI := range ee.Verts {
				if _, ok := ptMap[ptI]; !ok {
					pts = append(pts, ptI)
					goto NEXT_TRI
				}
			}
		}
	NEXT_TRI:
	}
	return
}

type TriMesh struct {
	Tris     []*Tri
	Points   []Point
	TriGraph *TriGraphNode // Graph used to determine which triangle a new point is inside
}

type TriGraphNode struct { // Tracks nested triangles during formation by using "point inside triangle" approach
	Triangle *Tri
	Children []*TriGraphNode
}

type Tri struct {
	Edges     []*Edge
	Reversals [3]int
	TGN       *TriGraphNode
}

func (tm *TriMesh) NewTri(edgesO ...*Edge) (tri *Tri) {
	if len(edgesO) < 2 || len(edgesO) > 3 {
		panic("need at least two and no more than three connected edges to make a triangle")
	}
	if len(edgesO) == 2 {
		edgesO = append(edgesO, NewEdgeFromEdges(edgesO))
	}
	// Check whether triangle has the correct number of unique vertices
	numVerts := len(getUniqueVertices(edgesO[0], edgesO[1], edgesO[2]))
	if numVerts != 3 {
		panic(fmt.Errorf("wrong number of vertices, need 3, have %d\n", numVerts))
	}
	tri = &Tri{}
	for _, e := range edgesO {
		tri.AddEdge(e)
	}
	// Orient the edges CCW
	//fmt.Println(tm.PrintTri(tri))
	tm.orientEdges(tri)
	return
}

func (tm *TriMesh) orientEdges(tri *Tri) {
	var (
		edges = tri.Edges
		pts   = tm.Points
	)
	// First connect each edge, reversing when needed to get the correct matching vertex
	orientFirstTwoEdges := func(e1, e2 *Edge) (r1, r2 int) {
		//Orient edges so that they connect the right of e1 to the left of e2
		switch {
		case e1.Verts[1] == e2.Verts[0]:
			// No reversals needed
			return
		case e1.Verts[1] == e2.Verts[1]:
			r2 = 1
			return
		case e1.Verts[0] == e2.Verts[0]:
			r1 = 1
			return
		case e1.Verts[0] == e2.Verts[1]:
			r1, r2 = 1, 1
			return
		default:
			err := fmt.Errorf("unable to connect edges: e1 Verts, e2 Verts = %v, %v\n", e1.Verts, e2.Verts)
			panic(err)
		}
		return
	}
	orientLastEdge := func(e2, e3 *Edge, e2Reversal int) (r3 int) {
		// If e2 must be reversed to connect it's 0 vertex to e1's 1 vertex, output a "1", otherwise "0"
		switch {
		case e2.Verts[1-e2Reversal] == e2.Verts[0]:
			return
		case e2.Verts[1-e2Reversal] == e2.Verts[1]:
			r3 = 1
			return
		default:
			err := fmt.Errorf("unable to connect edges: e2 Verts, e3 Verts, e2Reversal = %v, %v, %d\n", e2.Verts, e3.Verts, e2Reversal)
			panic(err)
		}
		return
	}
	tri.Reversals[0], tri.Reversals[1] = orientFirstTwoEdges(edges[0], edges[1])
	tri.Reversals[2] = orientLastEdge(edges[1], edges[2], tri.Reversals[1])
	// The triangle is now connected simply, let's check the orientation (CCW or CW)
	signF := func(p [3]Point) float64 {
		return (p[0].X[0]-p[2].X[0])*(p[1].X[1]-p[2].X[1]) - (p[1].X[0]-p[2].X[0])*(p[0].X[1]-p[2].X[1])
	}
	var verts [3]Point
	verts[0] = pts[edges[0].Verts[0]]
	verts[1] = pts[edges[0].Verts[1]]
	verts[2] = pts[edges[1].Verts[1-tri.Reversals[1]]]
	if math.Signbit(signF(verts)) {
		for i := 0; i < 3; i++ {
			// Reverse all edges in tri to orient it CCW
			tri.Reversals[i] = 1 - tri.Reversals[i]
		}
		// Reverse edge traversal to match orientation
		var (
			newR [3]int
		)
		newE := make([]*Edge, 3)
		for i := 0; i < 3; i++ {
			newE[2-i] = tri.Edges[i]
			newR[2-i] = tri.Reversals[i]
		}
		tri.Edges = newE
		tri.Reversals = newR
	}
}

func (tri *Tri) GetVertices() (verts [3]int) {
	if len(tri.Edges) != 3 {
		panic("not enough edges")
	}
	verts[0] = tri.Edges[0].Verts[0+tri.Reversals[0]]
	verts[1] = tri.Edges[0].Verts[1-tri.Reversals[0]]
	verts[2] = tri.Edges[1].Verts[1-tri.Reversals[1]]
	return
}

func (tri *Tri) AddEdge(e *Edge) {
	e.AddTri(tri)
	tri.Edges = append(tri.Edges, e)
	return
}

func NewTriMesh(X, Y []float64) (tris *TriMesh) {
	pts := make([]Point, len(X))
	for i, x := range X {
		y := Y[i]
		pts[i].X[0] = x
		pts[i].X[1] = y
	}
	tris = &TriMesh{
		Tris:   nil,
		Points: pts,
	}
	return
}

func (tm *TriMesh) PrintTri(tri *Tri, labelO ...string) string {
	if tri == nil {
		return "nil triangle"
	}
	var (
		pts   = tm.Points
		v     = tri.GetVertices()
		label = "triangle"
	)
	if len(labelO) != 0 {
		label = labelO[0]
	}
	return fmt.Sprintf("%s: Edges: [%d,%d],[%d,%d],[%d,%d] - v1: %8.5f, v2: %8.5f, v3: %8.5f\n",
		label,
		tri.Edges[0].Verts[0], tri.Edges[0].Verts[1],
		tri.Edges[1].Verts[0], tri.Edges[1].Verts[1],
		tri.Edges[2].Verts[0], tri.Edges[2].Verts[1],
		pts[v[0]], pts[v[1]], pts[v[2]])
}

func (tm *TriMesh) AddBoundingTriangle(tri *Tri) {
	tm.orientEdges(tri)
	tm.TriGraph = &TriGraphNode{Triangle: tri}
	tri.TGN = tm.TriGraph
}

func (tm *TriMesh) TriContainsPoint(tri *Tri, pt Point) (contains bool) {
	var (
		verts = tri.GetVertices()
		pts   = tm.Points
		//v1, v2, v3 = pts[verts[0]], pts[verts[1]], pts[verts[2]]
		vPts = []Point{pts[verts[0]], pts[verts[1]], pts[verts[2]]}
	)
	// Fast no - bounding box check
	xmin, xmax := vPts[0].X[0], vPts[0].X[0]
	ymin, ymax := vPts[0].X[1], vPts[0].X[1]
	for _, vpt := range vPts {
		xmin = math.Min(xmin, vpt.X[0])
		ymin = math.Min(ymin, vpt.X[1])
		xmax = math.Max(xmax, vpt.X[0])
		ymax = math.Max(ymax, vpt.X[1])
	}
	if pt.X[0] < xmin || pt.X[0] > xmax || pt.X[1] < ymin || pt.X[1] > ymax {
		return false
	}
	// Edge test - is point on a triangle edge?
	if tm.WhichEdgeIsPointOn(pt.X[0], pt.X[1], tri) != -1 {
		// Point is on an edge of this tri
		return true
	}
	// Interior point test
	signF := func(p1, p2, p3 Point) float64 {
		return (p1.X[0]-p3.X[0])*(p2.X[1]-p3.X[1]) - (p2.X[0]-p3.X[0])*(p1.X[1]-p3.X[1])
	}
	b1 := math.Signbit(signF(pt, vPts[0], vPts[1]))
	b2 := math.Signbit(signF(pt, vPts[1], vPts[2]))
	b3 := math.Signbit(signF(pt, vPts[2], vPts[0]))
	return (b1 == b2) && (b2 == b3)
}

func (tm *TriMesh) getLeafTri(tgn *TriGraphNode, pt Point) (triLeaf *Tri, leafNode *TriGraphNode) {
	for _, tgnDown := range tgn.Children {
		tri := tgnDown.Triangle
		if tm.TriContainsPoint(tri, pt) {
			triLeaf, leafNode = tm.getLeafTri(tgnDown, pt)
			return
		}
	}
	if len(tgn.Children) == 0 && tm.TriContainsPoint(tgn.Triangle, pt) {
		triLeaf = tgn.Triangle
		leafNode = tgn
	}
	return
}

func (tm *TriMesh) IsPointOnEdge(X, Y float64, e *Edge) (onEdge bool) {
	var (
		pts        = tm.Points
		e1x, e1y   = pts[e.Verts[0]].X[0], pts[e.Verts[0]].X[1]
		e2x, e2y   = pts[e.Verts[1]].X[0], pts[e.Verts[1]].X[1]
		xmin, ymin = math.Min(e1x, e2x), math.Min(e1y, e2y)
		xmax, ymax = math.Max(e1x, e2x), math.Max(e1y, e2y)
		tol        = 1.e-6
	)
	// Fast no bounding box test
	if X < xmin || X > xmax {
		return
	}
	if Y < ymin || Y > ymax {
		return
	}
	length := func(x1, y1, x2, y2 float64) (l float64) {
		l = math.Sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1))
		return
	}
	eLength := length(e1x, e1y, e2x, e2y)
	leftLength := length(e1x, e1y, X, Y)
	rightLength := length(X, Y, e2x, e2y)
	if math.Abs(leftLength+rightLength-eLength) < tol*eLength {
		return true
	}
	return
}

func (tm *TriMesh) WhichEdgeIsPointOn(X, Y float64, tri *Tri) (edgeNumber int) {
	// Returns -1 if point is not on triangle edge
	edgeNumber = -1
	for i, e := range tri.Edges {
		if tm.IsPointOnEdge(X, Y, e) {
			return i
		}
	}
	return
}

func (tm *TriMesh) AddPoint(X, Y float64) {
	pt := Point{X: [2]float64{X, Y}}
	tm.Points = append(tm.Points, pt)
	ptI := len(tm.Points) - 1
	//Find a triangle containing the point
	baseTri, leafNode := tm.getLeafTri(tm.TriGraph, pt)
	if baseTri == nil {
		err := fmt.Errorf(
			"unable to add point to triangulation, point %v is outside %v",
			pt, tm.PrintTri(tm.TriGraph.Triangle, "bounding triangle"))
		panic(err)
	}
	eNumber := tm.WhichEdgeIsPointOn(X, Y, baseTri)
	if eNumber == -1 { // Point is inside base tri, make 3 new triangles by connecting vertices of base tri to pt
		// Remove base tri edges
		for _, e := range baseTri.Edges {
			e.DeleteTri(baseTri)
		}
		v := baseTri.GetVertices()
		e1 := NewEdge([2]int{ptI, v[0]})
		e2 := NewEdge([2]int{ptI, v[1]})
		e3 := NewEdge([2]int{ptI, v[2]})

		tri := tm.NewTri(e1, baseTri.Edges[0], e2)
		tm.AddTriToGraph(tri, leafNode)
		tri = tm.NewTri(e2, baseTri.Edges[1], e3)
		tm.AddTriToGraph(tri, leafNode)
		tri = tm.NewTri(e3, baseTri.Edges[2], e1)
		tm.AddTriToGraph(tri, leafNode)
		// Legalize edges opposing ptI
		tm.LegalizeEdge(baseTri.Edges[0], ptI)
		tm.LegalizeEdge(baseTri.Edges[1], ptI)
		tm.LegalizeEdge(baseTri.Edges[2], ptI)
	} else {
		e := baseTri.Edges[eNumber]
		oppoVerts := e.GetOpposingVertices()
		fmt.Printf("Opposing Vertices = %v\n", oppoVerts)
	}
}

func (tm *TriMesh) GetEdgeTriWithoutPoint(e *Edge, ptI int) (oppoTri *Tri) {
	// Get the triangle on the other side of the edge, relative to ptI
	for tri := range e.Tris {
		var found bool
		for _, vertI := range tri.GetVertices() {
			if vertI == ptI {
				found = true
			}
		}
		if !found {
			return tri
		}
	}
	return
}

func IsIllegalEdge(prX, prY, piX, piY, pjX, pjY, pkX, pkY float64) bool {
	/*
		pr is a new point for candidate triangle pi-pj-pr
		pi-pj is a shared edge between pi-pj-pk and pi-pj-pr
		if pr lies inside the circle defined by pi-pj-pk:
			- The edge pi-pj should be swapped with pr-pk to make two new triangles:
				pi-pr-pk and pj-pk-pr
	*/
	inCircle := func(ax, ay, bx, by, cx, cy, dx, dy float64) (inside bool) {
		// Calculate handedness, counter-clockwise is (positive) and clockwise is (negative)
		signBit := math.Signbit((bx-ax)*(cy-ay) - (cx-ax)*(by-ay))
		ax_ := ax - dx
		ay_ := ay - dy
		bx_ := bx - dx
		by_ := by - dy
		cx_ := cx - dx
		cy_ := cy - dy
		det := (ax_*ax_+ay_*ay_)*(bx_*cy_-cx_*by_) -
			(bx_*bx_+by_*by_)*(ax_*cy_-cx_*ay_) +
			(cx_*cx_+cy_*cy_)*(ax_*by_-bx_*ay_)
		if signBit {
			return det < 0
		} else {
			return det > 0
		}
	}
	return inCircle(piX, piY, pjX, pjY, pkX, pkY, prX, prY)
}

func (tm *TriMesh) LegalizeEdge(e *Edge, testPtI int) {
	var (
		pts = tm.Points
		tri = tm.GetEdgeTriWithoutPoint(e, testPtI)
	)
	if tri == nil || e.IsImmovable { // There is no opposing tri, edge may be a boundary
		return
	}
	prX, prY := pts[testPtI].X[0], pts[testPtI].X[1]
	v := tri.GetVertices()
	p1x, p2x, p3x := pts[v[0]].X[0], pts[v[1]].X[0], pts[v[2]].X[0]
	p1y, p2y, p3y := pts[v[0]].X[1], pts[v[1]].X[1], pts[v[2]].X[1]
	if IsIllegalEdge(prX, prY, p1x, p1y, p2x, p2y, p3x, p3y) {
		//fmt.Printf("illegal edge, flipping\n")
		//flip edge, update leaves of trigraph
		tm.flipEdge(e)
		// Call recursively for the other two edges in the (formerly) opposing tri
		for _, eee := range tri.Edges {
			if eee != e {
				tm.LegalizeEdge(eee, testPtI)
			}
		}
	}
	return
}

func (tm *TriMesh) flipEdge(e *Edge) {
	// Reformulate the pair of triangles adjacent to edge into two new triangles connecting the opposing vertices
	if len(e.Tris) != 2 || e.IsImmovable { // Not able to flip edge
		fmt.Printf("unable to flip edge, #tris = %d, isImmovable = %v\n", len(e.Tris), e.IsImmovable)
		return
	}
	vv := [2][3]int{}
	tris := [2]*Tri{}
	var ii int
	for tri := range e.Tris {
		vv[ii] = tri.GetVertices()
		tris[ii] = tri
		ii++
	}
	eMap := make(map[*Edge]struct{}) // a bucket of edges for use in forming two new triangles
	for tri := range e.Tris {
		for _, ee := range tri.Edges {
			eMap[ee] = struct{}{}
		}
	}
	delete(eMap, e) // Remove "illegal" edge from bucket prior to triangle formation
	// Get opposing points
	getPtsExclEdge := func(verts [3]int) (op1 int) {
		for _, val := range verts {
			if val != e.Verts[0] && val != e.Verts[1] {
				op1 = val
				return
			}
		}
		panic("unable to find opposing point")
	}
	findConnectedEdge := func(ptI int) (ee *Edge) {
		for ee = range eMap {
			for _, vI := range ee.Verts {
				if vI == ptI {
					delete(eMap, ee)                  // remove edge from bucket
					ee.Tris = make(map[*Tri]struct{}) // reset connected tris prior to reuse
					return ee
				}
			}
		}
		panic("unable to find connected edge")
	}
	// Form new edge from points opposing "illegal" edge
	eNew := NewEdge([2]int{getPtsExclEdge(vv[0]), getPtsExclEdge(vv[1])})

	// Form first (of 2) new triangles
	pt1 := eNew.Verts[0]
	e1 := findConnectedEdge(pt1)
	var pt2 int
	if e1.Verts[0] == pt1 {
		pt2 = e1.Verts[1]
	} else {
		pt2 = e1.Verts[0]
	}
	e2 := findConnectedEdge(pt2)
	triNew1 := tm.NewTri(eNew, e1, e2)
	tm.AddTriToGraph(triNew1, tris[0].TGN)
	/*
		tri := tm.NewTri(e1, baseTri.Edges[0], e2)
		tm.AddTriToGraph(tri, leafNode)
	*/
	// Form second (of 2) new triangles
	pt1 = eNew.Verts[1]
	e1 = findConnectedEdge(pt1)
	if e1.Verts[0] == pt1 {
		pt2 = e1.Verts[1]
	} else {
		pt2 = e1.Verts[0]
	}
	e2 = findConnectedEdge(pt2)
	triNew2 := tm.NewTri(eNew, e1, e2)
	tm.AddTriToGraph(triNew2, tris[1].TGN)
}

func (tm *TriMesh) extractFinishedTris() {
	var (
		extractTris func(tgn *TriGraphNode)
	)
	// Recurse the Trigraph to extract leaves and load finished tri slice
	extractTris = func(tgn *TriGraphNode) {
		if len(tgn.Children) == 0 {
			tm.Tris = append(tm.Tris, tgn.Triangle)
			return
		} else {
			for _, node := range tgn.Children {
				extractTris(node)
			}
		}
	}
	extractTris(tm.TriGraph)
}

func (tm *TriMesh) ToGraphMesh() (trisOut graphics2D.TriMesh) {
	if len(tm.Tris) == 0 {
		tm.extractFinishedTris()
		if len(tm.Tris) == 0 {
			return
		}
	}
	pts := make([]graphics2D.Point, len(tm.Points))
	for i, pt := range tm.Points {
		pts[i].X[0] = float32(pt.X[0])
		pts[i].X[1] = float32(pt.X[1])
	}
	tris := make([]graphics2D.Triangle, len(tm.Tris))
	for i, tri := range tm.Tris {
		pts := tri.GetVertices()
		tris[i].Nodes[0] = int32(pts[0])
		tris[i].Nodes[1] = int32(pts[1])
		tris[i].Nodes[2] = int32(pts[2])
	}
	Attributes := make([][]float32, len(tm.Tris))
	for i, tri := range tm.Tris {
		for _, e := range tri.Edges {
			value := float32(0.5)
			if e.IsImmovable {
				value = float32(1.0)
			}
			Attributes[i] = append(Attributes[i], value)
		}
	}
	trisOut = graphics2D.TriMesh{
		BaseGeometryClass: graphics2D.BaseGeometryClass{
			Geometry: pts,
		},
		Triangles:  tris,
		Attributes: Attributes,
	}
	return
}

func (tm *TriMesh) AddTriToGraph(tri *Tri, leaf *TriGraphNode) {
	if leaf != nil { // Root node
		tgn := &TriGraphNode{Triangle: tri}
		leaf.Children = append(leaf.Children, tgn)
		tri.TGN = tgn
	}
}
