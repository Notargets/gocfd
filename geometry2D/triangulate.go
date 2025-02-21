package geometry2D

import (
	"fmt"
	"math"
)

type PointG struct {
	X [2]float64
}

func (pt *PointG) Print() (output string) {
	output = fmt.Sprintf("[%8.5f,%8.5f]", pt.X[0], pt.X[1])
	return
}

type Edge struct {
	Tris        map[*Tri]struct{} // Associated triangles
	Verts       [2]int
	IsImmovable bool // Is a BC edge
}

func (e *Edge) containsIndex(pI int) (contains bool) {
	// fmt.Printf("Verts = %v, Index = %d\n", e.Verts, pI)
	if e.Verts[0] == pI || e.Verts[1] == pI {
		return true
	}
	return
}

func (tm *TriMeshG) NewEdge(verts [2]int, isImmovableO ...bool) (e *Edge) {
	var (
		isImmovable bool
	)
	if len(isImmovableO) != 0 {
		isImmovable = isImmovableO[0]
	}
	// TODO: If this is necessary, use a data structure to make it not O(N^2)
	for _, ee := range tm.Edges {
		if ee.containsIndex(verts[0]) && ee.containsIndex(verts[1]) {
			// fmt.Printf("reusing edge: %s, inputIsImmovable: %v\n", ee.Print(), isImmovable)
			return ee
		}
	}
	e = newEdge(verts, isImmovable)
	tm.Edges = append(tm.Edges, e)
	return
}

func newEdge(verts [2]int, isImmovableO ...bool) (e *Edge) {
	if verts[0] == verts[1] {
		err := fmt.Errorf("degenerate edge, points are the same: [%d,%d]", verts[0], verts[1])
		panic(err)
	}
	e = &Edge{
		Verts: verts,
	}
	if len(isImmovableO) != 0 {
		e.IsImmovable = isImmovableO[0]
	}
	return
}

func (e *Edge) Print() string {
	var (
		label = "Movable"
	)
	if e.IsImmovable {
		label = "Fixed"
	}
	return fmt.Sprintf("[%d,%d]%s ", e.Verts[0], e.Verts[1], label)
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

func (e *Edge) DeleteTri(tri *Tri) {
	delete(e.Tris, tri)
}

func (e *Edge) AddTri(tri *Tri) {
	if e.Tris == nil {
		e.Tris = make(map[*Tri]struct{})
	}
	e.Tris[tri] = struct{}{}
}

func (e *Edge) getOpposingVertices() (pts []int, tris []*Tri) {
	// Get vertices opposing this edge from all tris with this edge, could be 1, 2 or 0 points (unconnected edge)
	for tri := range e.Tris {
		verts, _ := tri.GetVertices()
		for _, ptI := range verts {
			if ptI != e.Verts[0] && ptI != e.Verts[1] {
				pts = append(pts, ptI)
				tris = append(tris, tri)
				break
			}
		}
	}
	return
}

type TriMeshG struct {
	Tris     []*Tri
	Points   []PointG
	Edges    []*Edge
	TriGraph *TriGraphNode // Graph used to determine which triangle a new point is inside
}

type TriGraphNode struct { // Tracks nested triangles during formation by using "point inside triangle" approach
	Triangle *Tri
	Children []*TriGraphNode
}

func (tgn *TriGraphNode) GetDotOutput() (output []string) {
	var (
		descend func(tgn *TriGraphNode)
		visited = make(map[string]struct{})
	)
	output = append(output, "digraph Trigraph {")
	descend = func(ttgn *TriGraphNode) {
		myName := ttgn.Triangle.GetName()
		if _, ok := visited[myName]; ok {
			return
		}
		visited[myName] = struct{}{}
		for _, tt := range ttgn.Children {
			entry := fmt.Sprintf("\t_%s -> _%s", myName, tt.Triangle.GetName())
			output = append(output, entry)
		}
		for _, tt := range ttgn.Children {
			descend(tt)
		}
	}
	descend(tgn)
	output = append(output, "}")
	return
}

func (tgn *TriGraphNode) PrintDotOutput() {
	for _, line := range tgn.GetDotOutput() {
		fmt.Printf("%s\n", line)
	}
}

func (tgn *TriGraphNode) Print(labelO ...string) (str string) {
	var (
		tri   = tgn.Triangle
		label string
	)
	if len(labelO) != 0 {
		label = labelO[0]
	}
	return fmt.Sprintf("%s Node Edges: %s, %s, %s",
		label,
		tri.Edges[0].Print(),
		tri.Edges[1].Print(),
		tri.Edges[2].Print())
}

func (tgn *TriGraphNode) PrintAll() {
	var (
		descend func(t *TriGraphNode)
	)
	descend = func(t *TriGraphNode) {
		var label string
		if len(t.Children) == 0 {
			label = "**leaf** node"
		} else {
			label = "interior node"
		}
		fmt.Println(t.Print(label))
		for _, tt := range t.Children {
			descend(tt)
		}
	}
	descend(tgn)
}

type Tri struct {
	Edges     []*Edge
	Reversals [3]int
	TGN       *TriGraphNode
}

func (tri *Tri) GetName() (name string) {
	verts, fixed := tri.GetVertices()
	fixedName := func(f bool) string {
		if f {
			return "F"
		} else {
			return "M"
		}
	}
	name = fmt.Sprintf("%d%s_%d%s_%d%s", verts[0], fixedName(fixed[0]), verts[1], fixedName(fixed[1]), verts[2], fixedName(fixed[2]))
	return
}

func (tm *TriMeshG) NewTri(edges ...*Edge) (tri *Tri) {
	if !(len(edges) >= 2) {
		panic("need three connected edges or two connectable edges to make a triangle")
	}
	if len(edges) == 2 {
		// Construct needed edge
		// Find the two points not connected on the two input edges
		ptMap := make(map[int]int)
		for _, e := range edges {
			for _, vert := range e.Verts {
				ptMap[vert]++
			}
		}
		var verts [2]int
		var ii int
		var connected bool
		for key, count := range ptMap {
			if count == 1 {
				verts[ii] = key
				ii++
			} else if count == 2 {
				connected = true
			}
		}
		if !connected {
			panic("input edges not connected")
		}
		if ii != 2 {
			panic("unable to find two unconnected ends in the two input edges")
		}
		edges = append(edges, tm.NewEdge(verts))
	}
	// Check whether triangle has the correct number of unique vertices
	numVerts := len(getUniqueVertices(edges[0], edges[1], edges[2]))
	if numVerts != 3 {
		tm.TriGraph.PrintDotOutput()
		panic(fmt.Errorf("wrong number of vertices, need 3, have %d\nEdges: %s, %s, %s\n",
			numVerts, edges[0].Print(), edges[1].Print(), edges[2].Print()))
	}
	tri = &Tri{}
	for _, e := range edges {
		tri.AddEdge(e)
	}
	// Orient the edges CCW
	// fmt.Println(tm.PrintTri(tri))
	tm.orientEdges(tri)
	return
}

func (tm *TriMeshG) orientEdges(tri *Tri) {
	var (
		edges = tri.Edges
		pts   = tm.Points
	)
	// First connect each edge, reversing when needed to get the correct matching vertex
	orientFirstTwoEdges := func(e1, e2 *Edge) (r1, r2 int) {
		// Orient edges so that they connect the right of e1 to the left of e2
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
	signF := func(p [3]PointG) float64 {
		return (p[0].X[0]-p[2].X[0])*(p[1].X[1]-p[2].X[1]) - (p[1].X[0]-p[2].X[0])*(p[0].X[1]-p[2].X[1])
	}
	var verts [3]PointG
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

func (tri *Tri) GetVertices() (verts [3]int, fixed [3]bool) {
	if len(tri.Edges) != 3 {
		panic("not enough edges")
	}
	verts[0] = tri.Edges[0].Verts[0+tri.Reversals[0]]
	verts[1] = tri.Edges[0].Verts[1-tri.Reversals[0]]
	verts[2] = tri.Edges[1].Verts[1-tri.Reversals[1]]
	for _, e := range tri.Edges {
		if e.containsIndex(verts[0]) && e.containsIndex(verts[1]) {
			fixed[0] = e.IsImmovable
		}
		if e.containsIndex(verts[1]) && e.containsIndex(verts[2]) {
			fixed[1] = e.IsImmovable
		}
		if e.containsIndex(verts[2]) && e.containsIndex(verts[0]) {
			fixed[2] = e.IsImmovable
		}
	}
	return
}

func (tri *Tri) AddEdge(e *Edge) {
	e.AddTri(tri)
	tri.Edges = append(tri.Edges, e)
	return
}

func NewTriMesh(X, Y []float64) (tris *TriMeshG) {
	pts := make([]PointG, len(X))
	for i, x := range X {
		y := Y[i]
		pts[i].X[0] = x
		pts[i].X[1] = y
	}
	tris = &TriMeshG{
		Tris:   nil,
		Points: pts,
	}
	return
}

func (tm *TriMeshG) PrintTri(tri *Tri, labelO ...string) string {
	if tri == nil {
		return "nil triangle"
	}
	var (
		pts   = tm.Points
		v, _  = tri.GetVertices()
		label = "triangle"
	)
	if len(labelO) != 0 {
		label = labelO[0]
	}
	return fmt.Sprintf("%s: Edges: %s, %s, %s - v1: %8.5f, v2: %8.5f, v3: %8.5f\n",
		label,
		tri.Edges[0].Print(),
		tri.Edges[1].Print(),
		tri.Edges[2].Print(),
		pts[v[0]], pts[v[1]], pts[v[2]])
}

func (tm *TriMeshG) AddBoundingTriangle(tri *Tri) {
	tm.orientEdges(tri)
	tm.TriGraph = &TriGraphNode{Triangle: tri}
	tri.TGN = tm.TriGraph
}

func (tm *TriMeshG) triContainsPoint(tri *Tri, pt PointG, traceO ...bool) (contains bool) {
	var (
		verts, _ = tri.GetVertices()
		pts      = tm.Points
		// v1, v2, v3 = pts[verts[0]], pts[verts[1]], pts[verts[2]]
		vPts  = []PointG{pts[verts[0]], pts[verts[1]], pts[verts[2]]}
		trace = len(traceO) != 0 && traceO[0]
	)
	// fmt.Println("pt = ", pt, tm.PrintTri(tri, "inside TCP"))
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
		if trace {
			fmt.Printf("point %v outside bounding box\n", pt)
		}
		return false
	}
	// Edge test - is point on a triangle edge?
	if tm.whichEdgeIsPointOn(pt.X[0], pt.X[1], tri) != -1 {
		// PointG is on an edge of this tri
		if trace {
			fmt.Printf("point %v on an edge\n", pt)
		}
		return true
	}
	// Interior point test
	signF := func(p1, p2, p3 PointG) float64 {
		return (p1.X[0]-p3.X[0])*(p2.X[1]-p3.X[1]) - (p2.X[0]-p3.X[0])*(p1.X[1]-p3.X[1])
	}
	b1 := math.Signbit(signF(pt, vPts[0], vPts[1]))
	b2 := math.Signbit(signF(pt, vPts[1], vPts[2]))
	b3 := math.Signbit(signF(pt, vPts[2], vPts[0]))
	interior := (b1 == b2) && (b2 == b3)
	if interior && trace {
		fmt.Printf("point %v is interior\n", pt)
	}
	return interior
}

func (tm *TriMeshG) getLeafTri(tgn *TriGraphNode, pt PointG, traceO ...bool) (triLeaf *Tri, leafNode *TriGraphNode) {
	var (
		trace = len(traceO) != 0 && traceO[0]
	)
	if len(tgn.Children) == 0 {
		if trace {
			fmt.Println(tgn.Print("leaf"))
		}
		if tm.triContainsPoint(tgn.Triangle, pt, trace) {
			triLeaf = tgn.Triangle
			leafNode = tgn
		} else {
			panic("unable to confirm point is inside triangle on leaf node of graph")
		}
	} else {
		if trace {
			fmt.Printf("%s", tm.PrintTri(tgn.Triangle, "interior"))
			fmt.Printf("#Children = %d\n", len(tgn.Children))
		}
		for _, tgnDown := range tgn.Children {
			tri := tgnDown.Triangle
			if tm.triContainsPoint(tri, pt, trace) {
				triLeaf, leafNode = tm.getLeafTri(tgnDown, pt, trace)
				return
			}
		}
	}
	return
}

func (tm *TriMeshG) isPointOnEdge(X, Y float64, e *Edge) (onEdge bool) {
	var (
		pts        = tm.Points
		e1x, e1y   = pts[e.Verts[0]].X[0], pts[e.Verts[0]].X[1]
		e2x, e2y   = pts[e.Verts[1]].X[0], pts[e.Verts[1]].X[1]
		xmin, ymin = math.Min(e1x, e2x), math.Min(e1y, e2y)
		xmax, ymax = math.Max(e1x, e2x), math.Max(e1y, e2y)
		tol        = 1.e-09
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
	if math.Abs(leftLength+rightLength-eLength) < tol {
		return true
	}
	return
}

func (tm *TriMeshG) whichEdgeIsPointOn(X, Y float64, tri *Tri) (edgeNumber int) {
	// Returns -1 if point is not on triangle edge
	edgeNumber = -1
	for i, e := range tri.Edges {
		if tm.isPointOnEdge(X, Y, e) {
			return i
		}
	}
	return
}

func (tm *TriMeshG) AddPoint(X, Y float64, traceO ...bool) {
	var (
		trace                = len(traceO) != 0 && traceO[0]
		legalize1, legalize2 = true, true
	)
	pr := PointG{X: [2]float64{X, Y}}
	tm.Points = append(tm.Points, pr)
	pR := len(tm.Points) - 1
	// Find a triangle containing the point
	baseTri, leafNode := tm.getLeafTri(tm.TriGraph, pr)
	if baseTri == nil {
		tm.getLeafTri(tm.TriGraph, pr, true) // call again with tracing
		// tm.TriGraph.PrintAll()
		err := fmt.Errorf(
			"unable to add point to triangulation, point %v is outside %v",
			pr, tm.PrintTri(tm.TriGraph.Triangle, "bounding triangle"))
		panic(err)
	}
	eNumber := tm.whichEdgeIsPointOn(X, Y, baseTri)
	if eNumber == -1 { // PointG is inside base tri, make 3 new triangles by connecting vertices of base tri to pr
		// Remove base tri from edges, it will be replaced with new tris
		for _, e := range baseTri.Edges {
			e.DeleteTri(baseTri)
		}
		v, _ := baseTri.GetVertices()
		e1 := tm.NewEdge([2]int{pR, v[0]})
		e2 := tm.NewEdge([2]int{pR, v[1]})
		e3 := tm.NewEdge([2]int{pR, v[2]})
		/*
			fmt.Printf("About to split triangle into three\n%sPoint: %s\n%s\n%s\n%s\n",
				tm.PrintTri(baseTri, "baseTri"),
				tm.Points[pR].Print(),
				e1.Print(),
				e2.Print(),
				e3.Print(),
			)
		*/
		tri := tm.NewTri(e1, e2)
		tm.addTriToGraph(tri, leafNode)
		tri = tm.NewTri(e2, e3)
		tm.addTriToGraph(tri, leafNode)
		tri = tm.NewTri(e3, e1)
		tm.addTriToGraph(tri, leafNode)
		// Legalize edges opposing pR
		if legalize1 {
			tm.legalizeEdge(baseTri.Edges[0], pR)
			tm.legalizeEdge(baseTri.Edges[1], pR)
			tm.legalizeEdge(baseTri.Edges[2], pR)
		}
	} else {
		edge := baseTri.Edges[eNumber]
		oppoVerts, oppoTris := edge.getOpposingVertices() // Must happen before removing tris from edge
		if trace {
			fmt.Printf("X, Y = [%8.5f,%8.5f], Edge = %s, Oppoverts = %v\n",
				X, Y, edge.Print(), oppoVerts)
		}
		// Remove base tris from central edge, they will be replaced with 4 new tris
		for tri := range edge.Tris {
			for _, ee := range tri.Edges {
				ee.DeleteTri(tri)
			}
		}
		var newTris [4]*Tri
		var numTris int
		for i, pLorK := range oppoVerts {
			oTri := oppoTris[i]
			if trace {
				fmt.Printf("pR, pLorK = %d, %d\n", pR, pLorK)
			}
			e1 := tm.NewEdge([2]int{pR, pLorK})
			var eForNewTri []*Edge
			for _, ee := range oTri.Edges {
				if ee.containsIndex(pLorK) { // We want only the edges that share pLorK
					eForNewTri = append(eForNewTri, ee)
				}
			}
			if len(eForNewTri) != 2 {
				err := fmt.Errorf("not enough edges, should be 2, is: %d", len(eForNewTri))
				panic(err)
			}
			// Split edge IJ into two new edges, I-R and J-R for two new triangles
			var e2 [2]*Edge
			e2[0] = tm.NewEdge([2]int{edge.Verts[0], pR}, edge.IsImmovable)
			e2[1] = tm.NewEdge([2]int{edge.Verts[1], pR}, edge.IsImmovable)
			for ii, vI := range edge.Verts { // Form two new triangles from pR->pLorK->pIorJ
				var ee *Edge
				if eForNewTri[0].containsIndex(vI) {
					ee = eForNewTri[0]
				} else {
					ee = eForNewTri[1]
				}
				if trace {
					fmt.Printf("e1, e2, e3 = %s, %s, %s\n", e1.Print(), e2[ii].Print(), ee.Print())
				}
				newTris[numTris] = tm.NewTri(e2[ii], e1, ee)
				if trace {
					fmt.Println(tm.PrintTri(newTris[numTris]), "new tri")
				}
				tm.addTriToGraph(newTris[numTris], leafNode)
				numTris++
			}
		}
		if legalize2 {
			// Legalize the outer boundary edges of baseTri
			// fmt.Printf("Working on %d new tris...\n", numTris)
			for j := 0; j < numTris; j++ {
				for _, ee := range newTris[j].Edges {
					if !ee.containsIndex(pR) { // We only want the edge opposite of pR
						// fmt.Printf("legalizing: %s\n", ee.Print())
						tm.legalizeEdge(ee, pR)
					}
				}
			}
		}
	}
}

func (tm *TriMeshG) getEdgeTriWithoutPoint(e *Edge, ptI int) (oppoTri *Tri) {
	// Get the triangle on the other side of the edge, relative to ptI
	for tri := range e.Tris {
		var found bool
		verts, _ := tri.GetVertices()
		for _, vertI := range verts {
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

func isIllegalEdge(prX, prY, piX, piY, pjX, pjY, pkX, pkY float64) bool {
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

func (tm *TriMeshG) legalizeEdge(e *Edge, testPtI int) {
	var (
		pts = tm.Points
		tri = tm.getEdgeTriWithoutPoint(e, testPtI)
	)
	if tri == nil || e.IsImmovable { // There is no opposing tri, edge may be a boundary
		return
	}
	prX, prY := pts[testPtI].X[0], pts[testPtI].X[1]
	v, _ := tri.GetVertices()
	p1x, p2x, p3x := pts[v[0]].X[0], pts[v[1]].X[0], pts[v[2]].X[0]
	p1y, p2y, p3y := pts[v[0]].X[1], pts[v[1]].X[1], pts[v[2]].X[1]
	if isIllegalEdge(prX, prY, p1x, p1y, p2x, p2y, p3x, p3y) {
		// fmt.Printf("illegal edge, flipping\n")
		// flip edge, update leaves of trigraph
		tm.flipEdge(e)
		// Call recursively for the other two edges in the (formerly) opposing tri
		for _, eee := range tri.Edges {
			if eee != e {
				tm.legalizeEdge(eee, testPtI)
			}
		}
	}
	return
}

func (tm *TriMeshG) flipEdge(e *Edge) {
	// Reformulate the pair of triangles adjacent to edge into two new triangles connecting the opposing vertices
	if len(e.Tris) != 2 || e.IsImmovable { // Not able to flip edge
		fmt.Printf("unable to flip edge, #tris = %d, isImmovable = %v\n", len(e.Tris), e.IsImmovable)
		return
	}
	vv := [2][3]int{}
	tris := [2]*Tri{}
	var ii int
	for tri := range e.Tris {
		vv[ii], _ = tri.GetVertices()
		tris[ii] = tri
		ii++
	}
	eMap := make(map[*Edge]struct{}) // a bucket of edges for use in forming two new triangles
	for tri := range e.Tris {
		for _, ee := range tri.Edges {
			eMap[ee] = struct{}{}
			delete(ee.Tris, tri) // delete triangles on illegal edge from surrounding edges
		}
	}
	delete(eMap, e) // Remove "illegal" edge from bucket prior to triangle formation
	// Get opposing points
	getPtsExclEdge := func(verts [3]int) (op1 int) {
		for _, val := range verts {
			if !e.containsIndex(val) {
				return val
			}
		}
		panic("unable to find opposing point")
	}
	findConnectedEdge := func(ptI int) (ee *Edge) {
		for ee = range eMap {
			if ee.containsIndex(ptI) {
				delete(eMap, ee) // remove edge from bucket
				// fmt.Printf("Number of tris (%d) on edge %s\n", len(ee.Tris), ee.Print())
				return ee
			}
		}
		panic("unable to find connected edge")
	}
	// Form new edge from points opposing "illegal" edge
	eNew := tm.NewEdge([2]int{getPtsExclEdge(vv[0]), getPtsExclEdge(vv[1])})

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

	// tris[0].TGN.Triangle = triNew2
	// tris[1].TGN.Triangle = triNew1
	tm.addTriToGraph(triNew1, tris[0].TGN)
	tm.addTriToGraph(triNew1, tris[1].TGN)
	tm.addTriToGraph(triNew2, tris[0].TGN)
	tm.addTriToGraph(triNew2, tris[1].TGN)
}

func (tm *TriMeshG) extractFinishedTris() {
	var (
		extractTris func(tgn *TriGraphNode)
		triMap      = make(map[*Tri]struct{})
	)
	// Recurse the Trigraph to extract leaves and load finished tri slice
	extractTris = func(tgn *TriGraphNode) {
		if len(tgn.Children) == 0 {
			triMap[tgn.Triangle] = struct{}{}
			return
		} else {
			for _, node := range tgn.Children {
				extractTris(node)
			}
		}
	}
	extractTris(tm.TriGraph)
	for tri := range triMap {
		tm.Tris = append(tm.Tris, tri)
	}
}

func (tm *TriMeshG) ToGraphMesh() (trisOut TriMesh) {
	if len(tm.Tris) == 0 {
		tm.extractFinishedTris()
		if len(tm.Tris) == 0 {
			return
		}
	}
	pts := make([]Point, len(tm.Points))
	for i, pt := range tm.Points {
		pts[i].X[0] = float32(pt.X[0])
		pts[i].X[1] = float32(pt.X[1])
	}
	tris := make([]Triangle, len(tm.Tris))
	Attributes := make([][]float32, len(tm.Tris))
	for i, tri := range tm.Tris {
		pts, fixed := tri.GetVertices()
		for ptI := 0; ptI < 3; ptI++ {
			tris[i].Nodes[ptI] = int32(pts[ptI])
			value := float32(0.5)
			if fixed[ptI] {
				value = float32(1.0)
			}
			Attributes[i] = append(Attributes[i], value)
		}
	}
	trisOut = TriMesh{
		BaseGeometryClass: BaseGeometryClass{
			Geometry: pts,
		},
		Triangles:  tris,
		Attributes: Attributes,
	}
	return
}

func (tm *TriMeshG) addTriToGraph(tri *Tri, leaf *TriGraphNode) {
	if leaf != nil { // Root node
		if tri.TGN == nil {
			tri.TGN = &TriGraphNode{Triangle: tri}
		}
		leaf.Children = append(leaf.Children, tri.TGN)
	}
}
