package geometry2D

import (
	"fmt"
	"math"
)

type BoundingBox struct {
	XMin [2]float32
	XMax [2]float32
}

func NewBoundingBox(Geometry []Point) (Box *BoundingBox) {
	if len(Geometry) == 0 {
		return nil
	}
	Box = new(BoundingBox)
	Box.XMin[0], Box.XMin[1] = Geometry[0].X[0], Geometry[0].X[1]
	Box.XMax[0], Box.XMax[1] = Geometry[0].X[0], Geometry[0].X[1]
	for _, point := range Geometry {
		for i := 0; i < 2; i++ {
			if point.X[i] < Box.XMin[i] {
				Box.XMin[i] = point.X[i]
			}
			if point.X[i] > Box.XMax[i] {
				Box.XMax[i] = point.X[i]
			}
		}
	}
	return Box
}
func BoundingBoxFromGoogleMapsZoom(vCenter Point, canvasWidth int, zoomLevel float64) (viewBox *BoundingBox) {
	/*
			Calculate longitude width from Google Maps zoomLevel
		    From: https://developers.google.com/maps/documentation/javascript/coordinates
	*/
	lonWidth := float64(canvasWidth/256) * 360. / math.Pow(2, math.Max(0, zoomLevel))
	halfDelta := 0.5 * float32(lonWidth)

	lowLeft := vCenter.Minus(&Point{X: [2]float32{
		halfDelta, halfDelta,
	}})
	upRight := vCenter.Plus(&Point{X: [2]float32{
		halfDelta, halfDelta,
	}})

	viewBox = NewBoundingBox([]Point{*lowLeft, *upRight})

	centroid := viewBox.Centroid()
	trans := vCenter.Minus(centroid)

	return viewBox.Translate(trans.X)
}

func (bb *BoundingBox) Centroid() (centroid *Point) {
	return &Point{X: [2]float32{
		0.5 * (bb.XMax[0] + bb.XMin[0]),
		0.5 * (bb.XMax[1] + bb.XMin[1]),
	}}
}
func (bb *BoundingBox) Scale(scale float32) (bbOut *BoundingBox) {
	bbOut = new(BoundingBox)
	for i := 0; i < 2; i++ {
		xRange := bb.XMax[i] - bb.XMin[i]
		centroid := bb.XMin[i] + 0.5*xRange
		bbOut.XMin[i] = scale*(bb.XMin[i]-centroid) + centroid
		bbOut.XMax[i] = scale*(bb.XMax[i]-centroid) + centroid
	}
	return bbOut
}
func (bb *BoundingBox) ScaleX(scale float32) (bbOut *BoundingBox) {
	bbOut = new(BoundingBox)
	*bbOut = *bb
	xRange := bb.XMax[0] - bb.XMin[0]
	centroid := bb.XMin[0] + 0.5*xRange
	bbOut.XMin[0] = scale*(bb.XMin[0]-centroid) + centroid
	bbOut.XMax[0] = scale*(bb.XMax[0]-centroid) + centroid
	return bbOut
}
func (bb *BoundingBox) ScaleY(scale float32) (bbOut *BoundingBox) {
	bbOut = new(BoundingBox)
	*bbOut = *bb
	xRange := bb.XMax[1] - bb.XMin[1]
	centroid := bb.XMin[1] + 0.5*xRange
	bbOut.XMin[1] = scale*(bb.XMin[1]-centroid) + centroid
	bbOut.XMax[1] = scale*(bb.XMax[1]-centroid) + centroid
	return bbOut
}

func (bb *BoundingBox) Translate(panX [2]float32) (bbOut *BoundingBox) {
	bbOut = new(BoundingBox)
	for i := 0; i < 2; i++ {
		bbOut.XMin[i] = bb.XMin[i] + panX[i]
		bbOut.XMax[i] = bb.XMax[i] + panX[i]
	}
	return bbOut
}
func (bb *BoundingBox) MoveToOrigin() (bbOut *BoundingBox) {
	bbOut = new(BoundingBox)
	for i := 0; i < 2; i++ {
		bbOut.XMin[i] = 0
		bbOut.XMax[i] = bb.XMax[i] - bb.XMin[i]
	}
	return bbOut
}
func (bb *BoundingBox) Grow(newBB *BoundingBox) {
	for i := 0; i < 2; i++ {
		bb.XMin[i] = float32(math.Min(float64(bb.XMin[i]), float64(newBB.XMin[i])))
		bb.XMax[i] = float32(math.Max(float64(bb.XMax[i]), float64(newBB.XMax[i])))
	}
}
func (bb *BoundingBox) Divide(denom *BoundingBox) (scaleX []float64) {
	scaleX = make([]float64, 2)
	for i := 0; i < 2; i++ {
		denomRange := denom.XMax[i] - denom.XMin[i]
		bbRange := bb.XMax[i] - bb.XMin[i]
		scaleX[i] = float64(bbRange / denomRange)
	}
	return scaleX
}
func (bb *BoundingBox) Outline() (pLine *PolyLine) {
	pt1 := Point{X: [2]float32{
		bb.XMin[0], bb.XMin[1],
	}}
	pt2 := Point{X: [2]float32{
		bb.XMax[0], bb.XMin[1],
	}}
	pt3 := Point{X: [2]float32{
		bb.XMax[0], bb.XMax[1],
	}}
	pt4 := Point{X: [2]float32{
		bb.XMin[0], bb.XMax[1],
	}}
	return NewPolyLine([]Point{pt1, pt2, pt3, pt4, pt1})
}
func (bb *BoundingBox) PointInside(point *Point) (within bool) {
	for ii := 0; ii < 2; ii++ {
		if point.X[ii] > bb.XMax[ii] || point.X[ii] < bb.XMin[ii] {
			return false
		}
	}
	return true
}

type GeometryInterface interface {
	GetGeometry() []Point
	SetGeometry([]Point) // Allows for external transformations
	GetBoundingBox() *BoundingBox
	Area() float64
	Centroid() *Point
}

type BaseGeometryClass struct {
	Box      *BoundingBox
	Geometry []Point
}

func (bg *BaseGeometryClass) GetGeometry() (geom []Point) {
	return bg.Geometry
}
func (bg *BaseGeometryClass) SetGeometry(geom []Point) {
	bg.Geometry = geom
}
func (bg *BaseGeometryClass) GetBoundingBox() (bb *BoundingBox) {
	return bg.Box
}
func (bg *BaseGeometryClass) Area() (area float64)  { return 0 }
func (bg *BaseGeometryClass) Centroid() (ct *Point) { return bg.Box.Centroid() }
func (bg *BaseGeometryClass) TransformLambertAzimuthal(radius float64, center *Point) {
	/*
		From: www.epsg.org/Portals/0/373-07-2.pdf
		Projects the spherical coordinates geometry (in degrees of lon,lat) to
		a circular representation of a sphere of radius R
		Areas computed in the resulting projection are correct

		We use the center coordinates of the transformation for the FE and FN ("False Easting"
		and "False Northing" coordinates, which are arbitrary
	*/
	lambda0 := float64(center.X[0]) * math.Pi / 180
	phi1 := float64(center.X[1]) * math.Pi / 180
	cosPhi1 := math.Cos(phi1)
	inverseCosPhi1 := 1 / cosPhi1
	radScale := math.Pi / 180
	for i := range bg.Geometry {
		lambda := float64(bg.Geometry[i].X[0]) * radScale
		phi := float64(bg.Geometry[i].X[1]) * radScale
		bg.Geometry[i].X[0] = float32(radius*cosPhi1*(lambda-lambda0)) + center.X[0]
		bg.Geometry[i].X[1] = float32(radius*math.Sin(phi)*inverseCosPhi1) + center.X[1]
	}
}
func (bg *BaseGeometryClass) InverseTransformLambertAzimuthal(radius float64, center *Point) {
	/*
		Inverse transform - note that we use the center coordinate as both the [FE,FN] and [lambda,phi]0
	*/
	phi1 := float64(center.X[1]) * math.Pi / 180
	cosPhi1 := math.Cos(phi1)
	inverseR := 1 / radius
	inverseRCosPhi1 := 1 / (radius * cosPhi1)
	degScale := 180 / math.Pi
	for i := range bg.Geometry {
		E := float64(bg.Geometry[i].X[0])
		N := float64(bg.Geometry[i].X[1])
		bg.Geometry[i].X[0] = center.X[0] + float32(degScale*inverseRCosPhi1*(E-float64(center.X[0])))
		bg.Geometry[i].X[1] = float32(degScale * math.Asin(cosPhi1*inverseR*(N-float64(center.X[1]))))
	}
}

type Point struct {
	X [2]float32
}

func NewPoint(x, y float32) *Point {
	a := new(Point)
	a.X = [2]float32{x, y}
	return a
}

func (pt *Point) GetGeometry() (geom []Point) {
	return []Point{*pt}
}
func (pt *Point) SetGeometry(geom []Point) {
	*pt = geom[0]
}
func (pt *Point) GetBoundingBox() (box *BoundingBox) {
	return nil
}
func (pt *Point) Area() (area float64) { return 0 }
func (pt *Point) Centroid() (centroid *Point) {
	return pt
}

func (pt *Point) Minus(rhs *Point) (res *Point) {
	return &Point{X: [2]float32{
		pt.X[0] - rhs.X[0],
		pt.X[1] - rhs.X[1],
	}}
}
func (pt *Point) Plus(rhs *Point) (res *Point) {
	return &Point{X: [2]float32{
		pt.X[0] + rhs.X[0],
		pt.X[1] + rhs.X[1],
	}}
}
func (pt *Point) Equal(rhs Point) bool {
	return pt.X[0] == rhs.X[0] && pt.X[1] == rhs.X[1]
}
func (pt *Point) Scale(scale Point) {
	pt.X[0] *= scale.X[0]
	pt.X[1] *= scale.X[1]
}

type Line struct {
	Box      *BoundingBox
	geometry [2]Point
}

func NewLine(pt1, pt2 Point) (line *Line) {
	x1, y1, x2, y2 :=
		float64(pt1.X[0]),
		float64(pt1.X[1]),
		float64(pt2.X[0]),
		float64(pt2.X[1])
	pl := new(Line)
	pl.Box = new(BoundingBox)
	pl.Box.XMin[0] = float32(math.Min(x1, x2))
	pl.Box.XMax[0] = float32(math.Max(x1, x2))
	pl.Box.XMin[1] = float32(math.Min(y1, y2))
	pl.Box.XMax[1] = float32(math.Max(y1, y2))
	pl.geometry = [2]Point{pt1, pt2}
	return pl
}

func (ln *Line) GetGeometry() (geom []Point) {
	return []Point{ln.geometry[0], ln.geometry[1]}
}
func (ln *Line) SetGeometry(geom []Point) {
	ln.geometry[0] = geom[0]
	ln.geometry[1] = geom[1]
}
func (ln *Line) GetBoundingBox() (box *BoundingBox) {
	return ln.Box
}
func (ln *Line) Area() (area float64) { return 0 }
func (ln *Line) Centroid() (centroid *Point) {
	return ln.Box.Centroid()
}

type Triangle struct {
	Nodes [3]int32
}

type QuadMesh struct {
	BaseGeometryClass
	Dimensions [2]int
	Attributes [][]float32
}

type TriMesh struct {
	BaseGeometryClass
	Triangles  []Triangle
	Attributes [][]float32
}

func (tm *TriMesh) Area() (area float64) {
	/*
		From: https://www.mathopenref.com/coordtrianglearea.html
	*/
	for _, tri := range tm.Triangles {
		x1 := float64(tm.Geometry[tri.Nodes[0]].X[0])
		y1 := float64(tm.Geometry[tri.Nodes[0]].X[1])
		x2 := float64(tm.Geometry[tri.Nodes[1]].X[0])
		y2 := float64(tm.Geometry[tri.Nodes[1]].X[1])
		x3 := float64(tm.Geometry[tri.Nodes[2]].X[0])
		y3 := float64(tm.Geometry[tri.Nodes[2]].X[1])
		area += 0.5 * (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
	}
	return area
}
func (tm *TriMesh) Centroid() (centroid *Point) {
	var area, a, xc, yc float64
	for _, tri := range tm.Triangles {
		x1 := float64(tm.Geometry[tri.Nodes[0]].X[0])
		y1 := float64(tm.Geometry[tri.Nodes[0]].X[1])
		x2 := float64(tm.Geometry[tri.Nodes[1]].X[0])
		y2 := float64(tm.Geometry[tri.Nodes[1]].X[1])
		x3 := float64(tm.Geometry[tri.Nodes[2]].X[0])
		y3 := float64(tm.Geometry[tri.Nodes[2]].X[1])
		a = 0.5 * (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
		xc += a * (x1 + x2 + x3) / 3
		yc += a * (y1 + y2 + y3) / 3
		area += a
	}
	centroid = &Point{X: [2]float32{
		float32(xc / area), float32(yc / area),
	}}
	return centroid
}

func (tm *TriMesh) AddQuadMesh(mesh *QuadMesh) error {
	iDim, jDim := mesh.Dimensions[0], mesh.Dimensions[1]
	Len := len(mesh.Geometry)
	if iDim*jDim != Len {
		return fmt.Errorf("dimension mismatch %d x %d is not equal to %d",
			iDim, jDim, Len)
	}

	tm.Geometry = mesh.Geometry
	tm.Attributes = mesh.Attributes

	/*
		Create the triangle mesh
	*/
	NodeNum := func(i, j, iDim int) (nodeNum int32) {
		return int32(i + j*iDim)
	}

	for j := 0; j < jDim-1; j++ {
		for i := 0; i < iDim-1; i++ {
			nn1 := NodeNum(i, j, iDim)
			nn2 := NodeNum(i+1, j, iDim)
			nn3 := NodeNum(i, j+1, iDim)
			nn4 := NodeNum(i+1, j+1, iDim)
			tm.Triangles = append(
				tm.Triangles,
				Triangle{[3]int32{nn1, nn2, nn3}},
			)
			tm.Triangles = append(
				tm.Triangles,
				Triangle{[3]int32{nn3, nn4, nn2}},
			)
		}
	}

	return nil
}

type PolyLine struct {
	BaseGeometryClass
}

func NewPolyLine(geom []Point) (pl *PolyLine) {
	p_pl := new(PolyLine)
	p_pl.Box = NewBoundingBox(geom)
	p_pl.Geometry = geom
	return p_pl
}

type Polygon struct {
	BaseGeometryClass
}

func NewPolygon(geom []Point) (poly *Polygon) {
	/*
		Close off the polygon if needed
	*/
	if !geom[len(geom)-1].Equal(geom[0]) {
		geom = append(geom, geom[0])
	}
	pPoly := new(Polygon)
	pPoly.Box = NewBoundingBox(geom)
	pPoly.Geometry = geom
	return pPoly
}
func NewNgon(centroid Point, radius float64, n int) (poly *Polygon) {
	nF := float64(n)
	angleInc := 2 * math.Pi / nF
	var geom []Point
	for i := 0; i < n; i++ {
		angle := 2*math.Pi - float64(i)*angleInc // Generate in counterclockwise order for positive normal
		geom = append(geom,
			*centroid.Plus(&Point{X: [2]float32{
				float32(math.Sin(angle) * radius),
				float32(math.Cos(angle) * radius),
			}}))
	}
	return NewPolygon(geom)
}
func NewNgonGivenArea(centroid Point, area float64, n int) (poly *Polygon) {
	area = math.Max(-area, area)
	nF := float64(n)
	angle := 2 * math.Pi / nF
	radius := math.Sqrt(2 * area / (nF * math.Sin(angle)))
	return NewNgon(centroid, radius, n)
}

func (pg *Polygon) Centroid() (centroid *Point) {
	/*
		From: https://en.wikipedia.org/wiki/Centroid#Centroid_of_a_polygon
	*/
	centroid = &Point{X: [2]float32{0, 0}}
	area := pg.Area()
	ct := [2]float64{0, 0}
	for i := 0; i < len(pg.Geometry)-1; i++ {
		pt0 := pg.Geometry[i]
		pt1 := pg.Geometry[i+1]
		x0, y0 := float64(pt0.X[0]), float64(pt0.X[1])
		x1, y1 := float64(pt1.X[0]), float64(pt1.X[1])
		metric := x0*y1 - y0*x1
		ct[0] += (x0 + x1) * metric
		ct[1] += (y0 + y1) * metric
	}
	for i := 0; i < 2; i++ {
		centroid.X[i] = float32(ct[i] / (6 * area))
	}
	return centroid
}
func (pg *Polygon) Area() (area float64) {
	/*
		Algorithm: Green's theorem in the plane
	*/
	var a64 float64
	for i := 0; i < len(pg.Geometry)-1; i++ {
		pt0 := pg.Geometry[i]
		pt1 := pg.Geometry[i+1]
		x0, y0 := float64(pt0.X[0]), float64(pt0.X[1])
		x1, y1 := float64(pt1.X[0]), float64(pt1.X[1])
		a64 += x0*y1 - x1*y0
	}
	return 0.5 * a64
}

func (pg *Polygon) PointInside(point Point) (inside bool) {
	if !pg.Box.PointInside(&point) {
		return false
	}
	/*
		Algorithm:
		Winding Number from http://geomalgorithms.com/a03-_inclusion.html#wn_PnPoly()
		if wn = 0, the point is outside
	*/

	/*
		isLeft(): tests if a point is Left|On|Right of an infinite line.
		Input:  three points P0, P1, and P2
		Return:
			>0 for P2 left of the line through P0 and P1
			=0 for P2  on the line
			<0 for P2  right of the line
		See: Algorithm 1 "Area of Triangles and Polygons"
	*/
	isLeft := func(P0, P1, P2 Point) float32 {
		return (P1.X[0]-P0.X[0])*(P2.X[1]-P0.X[1]) -
			(P2.X[0]-P0.X[0])*(P1.X[1]-P0.X[1])
	}

	var wn int
	for i := 0; i < len(pg.Geometry)-1; i++ {
		pt0 := pg.Geometry[i]
		pt1 := pg.Geometry[i+1]
		if pt0.X[1] <= point.X[1] {
			if pt1.X[1] > point.X[1] {
				if isLeft(pt0, pt1, point) > 0 {
					wn++
				}
			}
		} else {
			if pt1.X[1] <= point.X[1] {
				if isLeft(pt0, pt1, point) < 0 {
					wn--
				}
			}
		}
	}
	return wn != 0
}

func LineLineIntersection(line1, line2 *Line) *Point {
	/*
		Check intersection of lines in two phases - calculate denominator, then full
		From: https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
	*/
	//epsilon := float32(math.SmallestNonzeroFloat32 * 1000000)
	x1, x2 := float64(line1.geometry[0].X[0]), float64(line1.geometry[1].X[0])
	y1, y2 := float64(line1.geometry[0].X[1]), float64(line1.geometry[1].X[1])
	x3, x4 := float64(line2.geometry[0].X[0]), float64(line2.geometry[1].X[0])
	y3, y4 := float64(line2.geometry[0].X[1]), float64(line2.geometry[1].X[1])

	denom := (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
	/*
		if float32(math.Abs(float64(denom))) < epsilon {
			return nil
		}
	*/
	m1 := x1*y2 - y1*x2
	m2 := x3*y4 - y3*x4
	Xnum := m1*(x3-x4) - (x1-x2)*m2
	Ynum := m1*(y3-y4) - (y1-y2)*m2
	return &Point{X: [2]float32{
		float32(Xnum / denom),
		float32(Ynum / denom),
	}}
}
