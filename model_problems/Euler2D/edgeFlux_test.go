package Euler2D

import (
	"fmt"
	"image/color"
	"math"
	"testing"

	"github.com/notargets/gocfd/model_problems/Euler2D/isentropic_vortex"

	"github.com/notargets/avs/chart2d"
	"github.com/notargets/avs/geometry"
	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/gocfd/DG2D"
)

type TestField uint8

const (
	NORMALSHOCKTESTM2 = TestField(iota)
	NORMALSHOCKTESTM5 = TestField(iota)
	NORMALSHOCKTESTM12
	FIXEDVORTEXTEST
	MOVINGVORTEXTEST
	RADIAL1TEST
	RADIAL2TEST
	RADIAL3TEST
	RADIAL4TEST
)

func (tf TestField) String() string {
	switch tf {
	case NORMALSHOCKTESTM2:
		return "NORMALSHOCKTESTM2"
	case NORMALSHOCKTESTM5:
		return "NORMALSHOCKTESTM5"
	case NORMALSHOCKTESTM12:
		return "NORMALSHOCKTESTM12"
	case FIXEDVORTEXTEST:
		return "FIXEDVORTEXTEST"
	case MOVINGVORTEXTEST:
		return "MOVINGVORTEXTEST"
	case RADIAL1TEST:
		return "RADIAL1TEST"
	case RADIAL2TEST:
		return "RADIAL2TEST"
	case RADIAL3TEST:
		return "RADIAL3TEST"
	case RADIAL4TEST:
		return "RADIAL4TEST"
	}
	return ""
}

func TestEdgeFlux(t *testing.T) {
	var (
		N = 3
	)
	if !testing.Verbose() {
		return
	}
	dfr := DG2D.NewDFR2D(N, false,
		"../../DG2D/test_data/test_10tris_centered.neu")
	gmRT := dfr.TriangulateRTElement()

	Np_RTBoundary := 3 * (1 + dfr.FluxElement.NpEdge)

	// fmt.Println("TriVerts = ", gmRT.TriVerts)

	// Index of edge segment pairs
	var targetPairs [][2]int
	targetPairs = make([][2]int, Np_RTBoundary-1)
	for i := 0; i < Np_RTBoundary-1; i++ { // Edges of boundary
		targetPairs[i] = [2]int{i, i + 1}
	}

	// Filter out "corner" triangles, those composed of only edge vertices
	var triList [][3]int
	for _, tri := range gmRT.TriVerts {
		edgeBound := int64(Np_RTBoundary)
		if !(tri[0] < edgeBound && tri[1] < edgeBound && tri[2] < edgeBound) {
			// ! Corner tri
			triList = append(triList, [3]int{int(tri[0]), int(tri[1]), int(tri[2])})
		}
	}

	// Build an index from vertex to the list of triangle indices that include that vertex.
	vertexToTris := make(map[int][]int)
	for idx, tri := range triList {
		for _, v := range tri {
			vertexToTris[v] = append(vertexToTris[v], idx)
		}
	}

	// For each target pair, intersect the lists of triangle indices for each vertex.
	triMap := make([][][3]int, len(targetPairs))
	for idx, pair := range targetPairs {
		a, b := pair[0], pair[1]
		trisA := vertexToTris[a]
		trisB := vertexToTris[b]

		// Build a set from one of the lists (choose the smaller list for efficiency).
		set := make(map[int]struct{})
		if len(trisA) < len(trisB) {
			for _, triIdx := range trisA {
				set[triIdx] = struct{}{}
			}
			for _, triIdx := range trisB {
				if _, found := set[triIdx]; found {
					triMap[idx] = append(triMap[idx], triList[triIdx])
				}
			}
		} else {
			for _, triIdx := range trisB {
				set[triIdx] = struct{}{}
			}
			for _, triIdx := range trisA {
				if _, found := set[triIdx]; found {
					triMap[idx] = append(triMap[idx], triList[triIdx])
				}
			}
		}
	}

	// Print result
	var edgeTris [][3]int
	for _, tris := range triMap {
		// fmt.Printf("Target Pair %v matches:\n", targetPairs[idx])
		for _, tri := range tris {
			fmt.Printf("  [%d:%d:%d]\n", tri[0], tri[1], tri[2])
			edgeTris = append(edgeTris, tri)
		}
	}
	fmt.Printf("Number of edge tris: %d\n", len(edgeTris))
	lines := make(map[color.RGBA][]float32)
	for i := 0; i < Np_RTBoundary; i++ {
		ip := i + 1
		if i == Np_RTBoundary-1 {
			ip = 0
		}
		lines[utils2.RED] = append(lines[utils2.RED],
			gmRT.XY[2*i], gmRT.XY[2*i+1],
			gmRT.XY[2*ip], gmRT.XY[2*ip+1])
	}
	for _, tri := range edgeTris {
		i1, i2, i3 := tri[0], tri[1], tri[2]
		lines[utils2.WHITE] = append(lines[utils2.WHITE],
			gmRT.XY[2*i1], gmRT.XY[2*i1+1],
			gmRT.XY[2*i2], gmRT.XY[2*i2+1],
			gmRT.XY[2*i2], gmRT.XY[2*i2+1],
			gmRT.XY[2*i3], gmRT.XY[2*i3+1],
			gmRT.XY[2*i3], gmRT.XY[2*i3+1],
			gmRT.XY[2*i1], gmRT.XY[2*i1+1])
	}
	var xyCross []float32
	for i := 0; i < dfr.SolutionElement.Np; i++ {
		x, y := dfr.SolutionElement.R.DataP[i], dfr.SolutionElement.S.DataP[i]
		xyCross = append(xyCross, float32(x), float32(y))
	}
	addCrossHairs(xyCross, utils2.GREEN, lines)
	plotLines(lines)
}

func addCrossHairs(xy []float32, col color.RGBA, lines map[color.RGBA][]float32) {
	var (
		lenXY = len(xy) / 2
		size  = float32(0.02)
	)
	for i := 0; i < lenXY; i++ {
		lines[col] = append(lines[col],
			xy[2*i]-size, xy[2*i+1],
			xy[2*i]+size, xy[2*i+1],
			xy[2*i], xy[2*i+1]-size,
			xy[2*i], xy[2*i+1]+size)
	}
}

func plotLines(lines map[color.RGBA][]float32) {
	var (
		xMin, xMax = float32(math.MaxFloat32), -float32(math.MaxFloat32)
		yMin, yMax = float32(math.MaxFloat32), -float32(math.MaxFloat32)
	)
	for _, line := range lines {
		xMin, xMax, yMin, yMax = getMinMax(line, xMin, xMax, yMin, yMax)
	}
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	for col, line := range lines {
		ch.AddLine(line, col)
	}
	for {
	}
}
func plotMesh(gm geometry.TriMesh) {
	var (
		xMin, xMax = float32(math.MaxFloat32), -float32(math.MaxFloat32)
		yMin, yMax = float32(math.MaxFloat32), -float32(math.MaxFloat32)
	)
	xMin, xMax, yMin, yMax = getMinMax(gm.XY, xMin, xMax, yMin, yMax)
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	ch.AddTriMesh(gm)
	for {
	}
}

func getMinMax(XY []float32, xi, xa, yi, ya float32) (xMin, xMax, yMin, yMax float32) {
	var (
		x, y  float32
		lenXY = len(XY) / 2
	)
	for i := 0; i < lenXY; i++ {
		x, y = XY[i*2+0], XY[i*2+1]
		if i == 0 {
			xMin = xi
			xMax = xa
			yMin = yi
			yMax = ya
		} else {
			if x < xMin {
				xMin = x
			}
			if x > xMax {
				xMax = x
			}
			if y < yMin {
				yMin = y
			}
			if y > yMax {
				yMax = y
			}
		}
	}
	return
}

func getFieldMinMax(field []float64) (fMin, fMax float64) {
	for i, f := range field {
		if i == 0 {
			fMin = f
			fMax = f
		}
		if f < fMin {
			fMin = f
		}
		if f > fMax {
			fMax = f
		}
	}
	return
}

func plotField(field []float64, gm geometry.TriMesh, FMin, FMax float64,
	xMM ...float64) {
	var xMin, xMax, yMin, yMax float32

	if len(xMM) == 4 {
		xMin, xMax = float32(xMM[0]), float32(xMM[1])
		yMin, yMax = float32(xMM[2]), float32(xMM[3])
	} else {
		xMin, xMax = -1, 1
		yMin, yMax = -1, 1
	}
	ch := chart2d.NewChart2D(xMin, xMax, yMin, yMax,
		1024, 1024, utils2.WHITE, utils2.BLACK)
	// Create a vector field including the three vertices
	var pField []float32
	var fMin, fMax float32
	fMin, fMax = math.MaxFloat32, -math.MaxFloat32
	pField = make([]float32, len(field))
	for i, f := range field {
		f32 := float32(f)
		if fMin > f32 {
			fMin = f32
		}
		if fMax < f32 {
			fMax = f32
		}
		pField[i] = float32(f)
	}
	vs := geometry.VertexScalar{
		TMesh:       &gm,
		FieldValues: pField,
	}
	fmt.Printf("Interpolated fMin: %f, fMax: %f\n", fMin, fMax)
	ch.AddShadedVertexScalar(&vs, float32(FMin), float32(FMax))
	ch.AddTriMesh(gm)
	line := []float32{0, -5, 0, 5, -5, 0, 5, 0}
	ch.AddLine(line, utils2.RED)
	for {
	}
}

func setTestField(X, Y []float64, tf TestField) (field []float64) {
	switch tf {
	case NORMALSHOCKTESTM12:
		field = setShockConditions(X, 1.2, 0)
	case NORMALSHOCKTESTM2:
		field = setShockConditions(X, 2, 0)
	case NORMALSHOCKTESTM5:
		field = setShockConditions(X, 5, 0)
	case FIXEDVORTEXTEST:
		iv := isentropic_vortex.NewIVortex(5, 0, 0, 1.4, 0)
		field = setIsoVortexConditions(X, Y, iv)
	case MOVINGVORTEXTEST:
		iv := isentropic_vortex.NewIVortex(5, 0, 0, 1.4)
		field = setIsoVortexConditions(X, Y, iv)
	case RADIAL1TEST:
		field = setRadial(X, Y, 1)
	case RADIAL2TEST:
		field = setRadial(X, Y, 2)
	case RADIAL3TEST:
		field = setRadial(X, Y, 3)
	case RADIAL4TEST:
		field = setRadial(X, Y, 4)
	}
	return
}

func setRadial(X, Y []float64, order float64) (field []float64) {
	var (
		Np = len(X)
	)
	field = make([]float64, Np*4)
	var x, y float64
	for i := 0; i < Np; i++ {
		x, y = X[i], Y[i]
		for n := 0; n < 4; n++ {
			field[i+n*Np] = math.Pow(x, order) + math.Pow(y, order)
		}
	}
	return
}

func setShockConditions(X []float64, Mach, Alpha float64) (field []float64) {
	var (
		U1, U2 = ShockConditions(Mach, Alpha)
		Np     = len(X)
	)
	// fmt.Printf("Mach: %.2f\n", Mach)
	// fmt.Printf("U1: ")
	// for n := 0; n < 4; n++ {
	// 	fmt.Printf(" %.2f ", U1[n])
	// }
	// fmt.Printf("\n")
	// fmt.Printf("U2: ")
	// for n := 0; n < 4; n++ {
	// 	fmt.Printf(" %.2f ", U2[n])
	// }
	// fmt.Printf("\n")
	field = make([]float64, Np*4)
	var x float64
	for i := 0; i < Np; i++ {
		x = X[i]
		for n := 0; n < 4; n++ {
			if x < 0 {
				field[i+n*Np] = U1[n]
			} else {
				field[i+n*Np] = U2[n]
			}
		}
	}
	return
}

func setIsoVortexConditions(X, Y []float64,
	iv *isentropic_vortex.IVortex) (field []float64) {
	var (
		Np = len(X)
	)
	field = make([]float64, Np*4)
	var x, y float64
	for i := 0; i < Np; i++ {
		x, y = X[i], Y[i]
		rho, rhoU, rhoV, rhoE := iv.GetStateC(0, x, y)
		field[i+0*Np] = rho
		field[i+1*Np] = rhoU
		field[i+2*Np] = rhoV
		field[i+3*Np] = rhoE
	}
	return
}

// ShockConditions computes post-shock properties given pre-shock Mach number (M1) and shock angle (alpha)
// where alpha = 0 corresponds to a normal shock.
func ShockConditions(M1, alpha float64) (U1, U2 [4]float64) {
	// Alpha = angle of incoming flow relative to the shock front (0 = normal shock).
	var (
		gamma = 1.4
	)

	// Convert angles to radians
	alphaRad := alpha * math.Pi / 180.0

	// Pre-shock conditions (non-dimensional)
	rho1 := 1.0
	p1 := 1.0
	u1 := M1
	v1 := 0.0 // Flow is initially aligned along x (free stream)

	// Resolve incoming velocity into normal/tangential components relative to the shock
	u1n := u1*math.Cos(alphaRad) + v1*math.Sin(alphaRad)  // Normal velocity component
	u1t := -u1*math.Sin(alphaRad) + v1*math.Cos(alphaRad) // Tangential component (unchanged)

	// Normal Mach number
	M1n := u1n / math.Sqrt(gamma*p1/rho1)

	// If subsonic normal component, no shock forms
	if M1n < 1.0 {
		En1 := p1/(gamma-1) + 0.5*rho1*(u1*u1+v1*v1)
		U1 = [4]float64{rho1, rho1 * u1, rho1 * v1, En1}
		U2 = U1
		return
	}

	// Rankine-Hugoniot relations (normal direction)
	rhoRatio := ((gamma + 1) * M1n * M1n) / ((gamma-1)*M1n*M1n + 2)
	pRatio := 1 + (2*gamma/(gamma+1))*(M1n*M1n-1)
	u2n_to_u1n := ((gamma-1)*M1n*M1n + 2) / ((gamma + 1) * M1n * M1n)

	// Post-shock normal values
	rho2 := rho1 * rhoRatio
	p2 := p1 * pRatio
	u2n := u1n * u2n_to_u1n

	// Tangential velocity is unchanged across the shock (slip condition for inviscid flow)
	u2t := u1t

	// Now transform back to x,y coordinates (undo rotation)
	u2 := u2n*math.Cos(alphaRad) - u2t*math.Sin(alphaRad)
	v2 := u2n*math.Sin(alphaRad) + u2t*math.Cos(alphaRad)

	// Total energy (consistent with post-shock pressure and velocity)
	En1 := p1/(gamma-1) + 0.5*rho1*(u1*u1+v1*v1)
	En2 := p2/(gamma-1) + 0.5*rho2*(u2*u2+v2*v2)

	// Pack conservative variables into U1 and U2
	U1 = [4]float64{rho1, rho1 * u1, rho1 * v1, En1}
	U2 = [4]float64{rho2, rho2 * u2, rho2 * v2, En2}

	return
}
