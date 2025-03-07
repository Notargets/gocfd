package DG2D

import (
	"fmt"
	"math"
	"testing"
	"time"

	"github.com/notargets/avs/chart2d"
	"github.com/notargets/avs/screen"
	utils2 "github.com/notargets/avs/utils"
	"github.com/notargets/gocfd/types"

	"github.com/notargets/gocfd/utils"
	"github.com/stretchr/testify/assert"
)

func TestDFR2D(t *testing.T) {
	// Basic test of interpolation matrix
	{
		N := 1
		dfr := NewDFR2D(N, false)
		el := dfr.SolutionElement
		fluxEl := dfr.FluxElement
		s := make([]float64, el.Np)
		for i := 0; i < el.Np; i++ {
			s[i] = float64(2 * i)
		}
		sM := utils.NewMatrix(el.Np, 1, s)
		// For each nodal location, interpolate a value (should equal the nodal function value)
		// Build an interpolating polynomial matrix using the nodal geometry
		interpM := el.JB2D.GetInterpMatrix(fluxEl.R, fluxEl.S)
		values := interpM.Mul(sM)
		// Verify the interpolated vals match the input solution values from the same [R,S]
		assert.InDeltaSlicef(t, s, values.DataP[0:3], 0.0000001, "")
		assert.InDeltaSlicef(t, s, values.DataP[3:6], 0.0000001, "")
		// After 2*NpInt points, the values have unknown expected interpolated values
	}
	// Test accuracy of interpolation
	{
		Nmax := 7
		for N := 1; N <= Nmax; N++ {
			dfr := NewDFR2D(N, false)
			el := dfr.SolutionElement
			fluxEl := dfr.FluxElement
			// Construct a 2D polynomial at the flux element geo locations, the first NpInt of which match the interior
			sFlux := make([]float64, fluxEl.Np)
			for i := 0; i < fluxEl.Np; i++ {
				sFlux[i] = utils.POW(fluxEl.R.DataP[i], N) + utils.POW(fluxEl.S.DataP[i], N)
			}
			// Copy the polynomial values from the first NpInt positions in the flux solution
			s := make([]float64, el.Np)
			for i := 0; i < fluxEl.NpInt; i++ {
				s[i] = sFlux[i]
			}
			sM := utils.NewMatrix(el.Np, 1, s)
			// Build an interpolating polynomial matrix using the nodal geometry
			interpM := el.JB2D.GetInterpMatrix(fluxEl.R, fluxEl.S)
			values := interpM.Mul(sM)
			// Verify the interpolated values match the input polynomial
			assert.InDeltaSlicef(t, sFlux, values.DataP, 0.00001, "")
		}
	}
	// Test point distribution
	{
		N := 1
		dfr := NewDFR2D(N, false)
		el := dfr.SolutionElement
		rt := dfr.FluxElement
		assert.InDeltaSlicef(t, rt.GetInternalLocations(rt.R.DataP),
			el.R.DataP, 0.000001, "")
		assert.InDeltaSlicef(t, rt.GetInternalLocations(rt.S.DataP),
			el.S.DataP, 0.000001, "")
		// t.Logf("Edge R = %v\n", rt.GetEdgeLocations(rt.R))
		// t.Logf("Edge S = %v\n", rt.GetEdgeLocations(rt.S))
	}
	// Test interpolation from solution points to flux points
	{
		N := 1
		dfr := NewDFR2D(N, false)
		el := dfr.SolutionElement
		rt := dfr.FluxElement
		assert.Equal(t, el.Np, rt.NpInt)
		solution := make([]float64, rt.NpInt)
		// Load some values into the solution space
		for i := 0; i < rt.NpInt; i++ {
			solution[i] = float64(i) / float64(rt.NpInt-1)
		}
		// Interpolate from interior to flux points
		sV := utils.NewMatrix(rt.NpInt, 1, solution)
		_ = dfr.FluxEdgeInterp.Mul(sV)
		// t.Logf("%s\n", fluxInterp.Print("fluxInterp"))
		// t.Logf("%s\n", sV.Print("sV"))
	}
	// Test triangulation
	{
		transpose := func(J []float64) (JT []float64) {
			JT = []float64{J[0], J[2], J[1], J[3]}
			return
		}
		multiply := func(J []float64, det float64, f [2]float64) (fm [2]float64) {
			fm[0], fm[1] = J[0]*f[0]+J[1]*f[1], J[2]*f[0]+J[3]*f[1]
			fm[0] /= det
			fm[1] /= det
			return
		}
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
		N := 1
		dfr := NewDFR2D(N, false, "test_data/test_tris_5.neu")
		/*
				Coordinates in test triangle case:
				1     0.000000000e+000    0.000000000e+000
				2     1.000000000e+000    0.000000000e+000
				3     0.500000000e+001    1.000000000e+001
				4    -1.000000000e+001    1.000000000e+001
			    Two tris:
				1, 2, 3
				1, 3, 4
				Five Edges:
				0-2: 		Shared between Tri 0 and Tri 1
				0-1, 2-3: 	BC Inflow
				1-2, 0-3: 	Unconnected
		*/
		trn := dfr.Tris
		{ // Test edge construction
			{ // 0-2: 		Shared between Tri 0 and Tri 1
				en := types.NewEdgeKey([2]int{0, 2}) // Edge formed by vertices 0 and 2
				e := trn.Edges[en]
				assert.Equal(t, uint8(2), e.NumConnectedTris)
				assert.Equal(t, uint32(0), e.ConnectedTris[0])
				assert.Equal(t, InternalEdgeNumber(2), e.ConnectedTriEdgeNumber[First])
				assert.Equal(t, uint32(1), e.ConnectedTris[1])
				assert.Equal(t, InternalEdgeNumber(0), e.ConnectedTriEdgeNumber[Second])
			}
			{ // 0-1, 2-3: 	BC Inflow
				en := types.NewEdgeKey([2]int{0, 1})
				e := trn.Edges[en]
				assert.Equal(t, types.BC_In, e.BCType)
				en = types.NewEdgeKey([2]int{2, 3})
				e = trn.Edges[en]
				assert.Equal(t, types.BC_In, e.BCType)
			}
			assert.Equal(t, dfr.Tris.EtoE[0], [3]int{-1, -1, 1})
			assert.Equal(t, dfr.Tris.EtoE[1], [3]int{0, -1, -1})
		}
		// Test Piola transform and jacobian
		for en, e := range dfr.Tris.Edges {
			assert.True(t, e.NumConnectedTris > 0)
			for conn := 0; conn < int(e.NumConnectedTris); conn++ {
				k := e.ConnectedTris[conn]
				revDir := bool(e.ConnectedTriDirection[conn])
				J, Jinv, Jdet := dfr.GetJacobian(int(k))
				Jt := transpose(J)       // Transpose is applied to normals, which are cross products of vectors
				JtInv := transpose(Jinv) // Transpose is applied to normals, which are cross products of vectors
				x1, x2 := GetEdgeCoordinates(en, revDir, dfr.VX, dfr.VY)
				dx := [2]float64{x2[0] - x1[0], x2[1] - x1[1]}
				normal := [2]float64{-dx[1], dx[0]}

				edgeNumber := e.ConnectedTriEdgeNumber[conn]
				{ // Check that a full edge scaled normal transforms into a full edge scaled transformed normal
					nxT := multiply(Jt, Jdet, normal)
					switch edgeNumber {
					case First: // A whole 1st edge is 2 units long in unit tri space
						assert.Equal(t, [2]float64{0, -2}, nxT)
					case Second: // A whole 2nd edge is 2*Sqrt(2) units long in unit tri space
						assert.Equal(t, [2]float64{2, 2}, nxT)
					case Third: // A whole 3rd edge is 2 units long in unit tri space
						assert.Equal(t, [2]float64{-2, 0}, nxT)
					}
					nxTT := multiply(JtInv, 1./Jdet, nxT)
					assert.InDeltaf(t, normal[0], nxTT[0], 0.00001, "")
					assert.InDeltaf(t, normal[1], nxTT[1], 0.00001, "")
				}
				{ // Check scaling factor ||n||, used in transforming face normals
					normal = normalize(normal)
					normal[0] *= e.IInII[conn] // scale ||n|| the vector coords on the edge prior to transform
					normal[1] *= e.IInII[conn] // scale ||n|| the vector coords on the edge prior to transform
					nxT := multiply(Jt, Jdet, normal)
					oosr2 := 1. / math.Sqrt(2)
					switch edgeNumber { // The transformed and scaled normal should be unit for each edge direction
					case First:
						assert.InDeltaf(t, 0, nxT[0], 0.000001, "")
						assert.InDeltaf(t, -1, nxT[1], 0.000001, "")
					case Second:
						assert.InDeltaf(t, oosr2, nxT[0], 0.000001, "")
						assert.InDeltaf(t, oosr2, nxT[1], 0.000001, "")
					case Third:
						assert.InDeltaf(t, -1, nxT[0], 0.000001, "")
						assert.InDeltaf(t, 0, nxT[1], 0.000001, "")
					}
				}
			}
		}
	}
}

func _TestPlotEquiTri(t *testing.T) {
	dfr := CreateEquiTriMesh(1, 47)
	if testing.Verbose() {
		PlotDFRElements(dfr)
	}
}

func TestDFRP3(t *testing.T) {
	dfr := CreateEquiTriMesh(3, 0)
	dfr.VX.Transpose().Print("VX")
	dfr.VY.Transpose().Print("VY")
	dfr.FluxX.Transpose().Print("FluxX")
	dfr.FluxY.Transpose().Print("FluxY")
	// gm := CreateGraphMesh(dfr)
	// _ = gm
}

func CreateEquiTriMesh(N int, angle float64, dfrO ...*DFR2D) (dfr *DFR2D) {
	var (
		scale   = float64(5)
		yHeight = scale * math.Sin(math.Pi*60./180.)
	)
	if len(dfrO) > 0 {
		dfr = dfrO[0]
	} else {
		dfr = NewDFR2D(N, false)
	}
	dfr.K = 1
	vx, vy := []float64{-scale * 0.5, scale * 0.5, 0.},
		[]float64{-yHeight / 3., -yHeight / 3, 2. * yHeight / 3.}
	dfr.VX, dfr.VY = utils.NewVector(3), utils.NewVector(3)
	for i := 0; i < dfr.VX.Len(); i++ {
		dfr.VX.DataP[i] = vx[i]*math.Cos(2.*math.Pi*angle/360.) +
			vy[i]*math.Sin(2.*math.Pi*angle/360.)
		dfr.VY.DataP[i] = -vx[i]*math.Sin(2.*math.Pi*angle/360.) +
			vy[i]*math.Cos(2.*math.Pi*angle/360.)
	}
	EToV := utils.NewMatrix(dfr.K, 3, []float64{0, 1, 2})
	dfr.ProcessGeometry(dfr.VX, dfr.VY, EToV, nil)
	return
}

func PlotDFRElements(dfr *DFR2D) {
	var Line []float32
	addVertex := func(x, y, width float64) {
		wo2 := float32(width / 2.)
		xl, yl := float32(x), float32(y)
		Line = append(Line, xl-wo2, yl, xl+wo2, yl)
		Line = append(Line, xl, yl-wo2, xl, yl+wo2)
	}
	addLineSegment := func(x1, y1, x2, y2 float64) {
		x1l, y1l := float32(x1), float32(y1)
		x2l, y2l := float32(x2), float32(y2)
		Line = append(Line, x1l, y1l, x2l, y2l)
	}

	ch := chart2d.NewChart2D(-1, 1, -1, 1, 1024, 1024, utils2.WHITE,
		utils2.BLACK, 0.8)
	Np := dfr.FluxElement.Np
	for k := 0; k < dfr.K; k++ {
		for i := 0; i < Np; i++ {
			addVertex(dfr.FluxX.At(i, k), dfr.FluxY.At(i, k), 0.025)
		}
	}
	for k := 0; k < dfr.K; k++ {
		verts := dfr.Tris.GetTriVerts(uint32(k))
		for ii := 0; ii < 3; ii++ {
			x1, y1 := dfr.VX.DataP[verts[ii]], dfr.VY.DataP[verts[ii]]
			ip := ii + 1
			if ii == 2 {
				ip = 0
			}
			x2, y2 := dfr.VX.DataP[verts[ip]], dfr.VY.DataP[verts[ip]]
			addLineSegment(x1, y1, x2, y2)
		}
	}
	ch.AddLine(Line, utils2.WHITE)
	ch.NewWindow("Unit Triangle", 0.9, screen.AUTO)
	Line = []float32{}
	addLineSegment(-1, -1, 1, -1)
	addLineSegment(1, -1, -1, 1)
	addLineSegment(-1, 1, -1, -1)
	ch.AddLine(Line, utils2.WHITE)
	time.Sleep(30 * time.Second)
}

func DivergenceCheck(t *testing.T, dfr *DFR2D) {
	// Test divergence
	checkSolution := func(dfr *DFR2D, Order int) (Fx, Fy, divCheck utils.Matrix) {
		var (
			rt  = dfr.FluxElement
			K   = dfr.K
			Np  = rt.Np
			ccf = float64(Order)
		)
		divCheck = utils.NewMatrix(K, Np)
		Fx, Fy = utils.NewMatrix(K, Np), utils.NewMatrix(K, Np)
		for k := 0; k < K; k++ {
			var (
				xD, yD   = dfr.FluxX.Col(k).DataP[0:Np], dfr.FluxY.Col(k).DataP[0:Np]
				fxD, fyD = Fx.DataP, Fy.DataP
				dcD      = divCheck.DataP
			)
			for n := 0; n < Np; n++ {
				ind := n + k*Np
				x, y := xD[n], yD[n]
				fxD[ind], fyD[ind] = utils.POW(x, Order), utils.POW(y, Order)
				dcD[ind] = ccf * (utils.POW(x, Order-1) + utils.POW(y, Order-1))
			}
		}
		return
	}
	{ // Check divergence for polynomial vector fields of order < N against analytical solution
		// dfr := NewDFR2D(N, false, "test_data/test_tris_5.neu")
		N := dfr.N
		rt := dfr.FluxElement
		// rt.Div.Print("Div")
		// rt.DivInt.Print("DivInt")
		for cOrder := 1; cOrder <= N; cOrder++ { // Run a test on polynomial flux vector fields up to Nth order
			t.Logf("checking RT order[%d]...", cOrder+1)
			Fx, Fy, divCheck := checkSolution(dfr, cOrder)
			// project the flux onto the RT basis directly
			Fp := dfr.ProjectFluxOntoRTSpace(Fx, Fy)
			for k := 0; k < dfr.K; k++ {
				var (
					Fpk  = Fp.Row(k).ToMatrix()
					Jdet = dfr.Jdet.Row(k).DataP[0]
				)
				divM := rt.Div.Mul(Fpk).Scale(1. / Jdet)
				// divM.Print("divM")
				// t.Logf("divCheck.Row(k) = %v\n", divCheck.Row(k).DataP)
				assert.InDeltaSlicef(t, divCheck.Row(k).DataP, divM.DataP, 0.00001, "")
			}
			// Now project flux onto untransformed normals using ||n|| scaling factor, divergence should be the same
			SetNormalFluxOnEdges(dfr, Fx, Fy, &Fp)
			for k := 0; k < dfr.K; k++ {
				var (
					Fpk  = Fp.Row(k).ToMatrix()
					Jdet = dfr.Jdet.Row(k).DataP[0]
				)
				divM := rt.Div.Mul(Fpk).Scale(1. / Jdet)
				assert.InDeltaSlicef(t, divCheck.Row(k).DataP, divM.DataP, 0.00001, "")
			}
			t.Logf("passed.\n")
		}
	}
}

func TestDivergence(t *testing.T) {
	for _, ang := range []float64{0, 15, 25, 45, 90, 180, 210, 270, 310} {
		fmt.Printf("Checking angle:%f ==================\n", ang)
		DivergenceCheck(t, CreateEquiTriMesh(2, ang))
	}
	DivergenceCheck(t, NewDFR2D(7, false, "test_data/test_tris_5.neu"))
}

func TestGradient(t *testing.T) {
	// Test gradient on Flux element points, derived from a solution field interpolated from solution pts to Flux pts
	dfr := NewDFR2D(3, false, "test_data/test_tris_6.neu")
	var (
		Kmax   = dfr.K
		NpInt  = dfr.SolutionElement.Np
		NpFlux = dfr.FluxElement.Np
		X, Y   = dfr.FluxX, dfr.FluxY
		Dr     = dfr.FluxDr
		Ds     = dfr.FluxDs
		Jinv   = dfr.Jinv
		// Jdet   = dfr.Jdet
		// Scalar fields, linear, quadratic, cubic
		XY1, XY2, XY3 = utils.NewMatrix(NpInt, Kmax), utils.NewMatrix(NpInt, Kmax), utils.NewMatrix(NpInt, Kmax)
		// Directional derivatives of scalar fields, linear, quadratic, cubic
		DX1, DX2, DX3 = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
		DY1, DY2, DY3 = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
	)
	assert.Equal(t, 2, Kmax)
	for k := 0; k < Kmax; k++ {
		for i := 0; i < NpFlux; i++ {
			ind := k + i*Kmax
			x, y := X.DataP[ind], Y.DataP[ind]
			DX1.DataP[ind] = 1         // Derivative of Linear field
			DX2.DataP[ind] = 2 * x     // Derivative of Quadratic field
			DX3.DataP[ind] = 3 * x * x // Derivative of Cubic field
			DY1.DataP[ind] = 1         // Derivative of Linear field
			DY2.DataP[ind] = 2 * y     // Derivative of Quadratic field
			DY3.DataP[ind] = 3 * y * y // Derivative of Cubic field
			if i < NpInt {             // Solution points [R,S] coords are identical to RT [R,S] up to Nintx2
				XY1.DataP[ind] = x + y         // Linear field
				XY2.DataP[ind] = x*x + y*y     // Quadratic field
				XY3.DataP[ind] = x*x*x + y*y*y // Cubic field
			}
		}
	}
	/*
		Test gradient in physical coordinates after interpolating from solution to flux pts and taking derivatives
		using the Lagrange element derivative operators
	*/
	{
		Qr, Qs := Dr.Mul(XY1), Ds.Mul(XY1) // In a single multiplication, interpolate from solution to flux pts and Dr/Ds
		Qx, Qy := utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)

		Qr2, Qs2 := Dr.Mul(XY2), Ds.Mul(XY2)
		Qx2, Qy2 := utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)

		Qr3, Qs3 := Dr.Mul(XY3), Ds.Mul(XY3)
		Qx3, Qy3 := utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
		transform := func(ddr, dds utils.Matrix, JinvD []float64, index int) (ddx, ddy float64) {
			var (
				ddrD, ddsD = ddr.DataP, dds.DataP
			)
			rx, ry, sx, sy := JinvD[0], JinvD[1], JinvD[2], JinvD[3]
			ddx = ddrD[index]*rx + ddsD[index]*sx
			ddy = ddrD[index]*ry + ddsD[index]*sy
			return
		}
		for k := 0; k < Kmax; k++ {
			var (
				JinvD = Jinv.DataP[4*k : 4*(k+1)]
			)
			for i := 0; i < NpFlux; i++ {
				ind := k + i*Kmax
				Qx.DataP[ind], Qy.DataP[ind] = transform(Qr, Qs, JinvD, ind)
				Qx2.DataP[ind], Qy2.DataP[ind] = transform(Qr2, Qs2, JinvD, ind)
				Qx3.DataP[ind], Qy3.DataP[ind] = transform(Qr3, Qs3, JinvD, ind)
			}
		}
		// Compare known analytical gradient to calculation
		assert.InDeltaSlicef(t, Qx.DataP, DX1.DataP, 0.000001, "err msg %s")
		assert.InDeltaSlicef(t, Qy.DataP, DY1.DataP, 0.000001, "err msg %s")
		assert.InDeltaSlicef(t, Qx2.DataP, DX2.DataP, 0.000001, "err msg %s")
		assert.InDeltaSlicef(t, Qy2.DataP, DY2.DataP, 0.000001, "err msg %s")
		assert.InDeltaSlicef(t, Qx3.DataP, DX3.DataP, 0.000001, "err msg %s")
		assert.InDeltaSlicef(t, Qy3.DataP, DY3.DataP, 0.000001, "err msg %s")
	}
	// Test Gradient derived from Raviart Thomas element, from solution field interpolated from solution pts
	{
		var (
			DXmd, DYmd   = dfr.DXMetric.DataP, dfr.DYMetric.DataP
			DOFX, DOFY   = utils.NewMatrix(NpFlux, Kmax), utils.NewMatrix(NpFlux, Kmax)
			DOFXd, DOFYd = DOFX.DataP, DOFY.DataP
		)
		// Create DX and DY for each field type in a loop and check against analytic values
		// var Uint, Uedge utils.Matrix
		fI := dfr.FluxEdgeInterp.Mul
		Uint := []utils.Matrix{XY1, XY2, XY3}
		Uedge := []utils.Matrix{fI(Uint[0]), fI(Uint[1]), fI(Uint[2])}
		DXCheck, DYCheck := []utils.Matrix{DX1, DX2, DX3}, []utils.Matrix{DY1, DY2, DY3}
		{
			for n := range Uint {
				// for n := 0; n < 3; n++ {
				var Un float64
				for k := 0; k < Kmax; k++ {
					for i := 0; i < NpFlux; i++ {
						ind := k + i*Kmax
						switch {
						case i < NpInt: // The first NpInt points are the solution element nodes
							Un = Uint[n].DataP[ind]
						case i >= NpInt && i < 2*NpInt: // The second NpInt points are duplicates of the first NpInt values
							Un = Uint[n].DataP[ind-NpInt*Kmax]
						case i >= 2*NpInt:
							Un = Uedge[n].DataP[ind-2*NpInt*Kmax] // The last 3*NpEdge points are the edges in [0-1,1-2,2-0] order
						}
						DOFXd[ind] = DXmd[ind] * Un
						DOFYd[ind] = DYmd[ind] * Un
					}
				}
				DX := dfr.FluxElement.Div.Mul(DOFX) // X Derivative, divergence x RT_DOF is X derivative for this DOF
				DY := dfr.FluxElement.Div.Mul(DOFY) // Y Derivative, divergence x RT_DOF is Y derivative for this DOF
				t.Logf("Order[%d] check ...", n+1)
				assert.Equal(t, len(DX.DataP), len(DXCheck[n].DataP))
				assert.Equal(t, len(DY.DataP), len(DYCheck[n].DataP))
				assert.InDeltaSlicef(t, DX.DataP, DXCheck[n].DataP, 0.000001, "err msg %s")
				assert.InDeltaSlicef(t, DY.DataP, DYCheck[n].DataP, 0.000001, "err msg %s")
				t.Logf("... validates\n")
			}
		}
	}
}

func SetNormalFluxOnEdges(dfr *DFR2D, Fx, Fy utils.Matrix, Fp *utils.Matrix) {
	var (
		Np       = dfr.FluxElement.Np
		fxD, fyD = Fx.DataP, Fy.DataP
		fpD      = Fp.DataP
		Nint     = dfr.FluxElement.NpInt
		Nedge    = dfr.FluxElement.NpEdge
	)
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
	for en, e := range dfr.Tris.Edges {
		for conn := 0; conn < int(e.NumConnectedTris); conn++ {
			var (
				k = int(e.ConnectedTris[conn])
			)
			revDir := bool(e.ConnectedTriDirection[conn])
			x1, x2 := GetEdgeCoordinates(en, revDir, dfr.VX, dfr.VY)
			dx := [2]float64{x2[0] - x1[0], x2[1] - x1[1]}
			normal := normalize([2]float64{-dx[1], dx[0]})
			normal[0] *= e.IInII[conn]
			normal[1] *= e.IInII[conn]
			edgeNumber := int(e.ConnectedTriEdgeNumber[conn])
			for n := 2*Nint + edgeNumber*Nedge; n < 2*Nint+(edgeNumber+1)*Nedge; n++ {
				ind := n + k*Np
				fpD[ind] = normal[0]*fxD[ind] + normal[1]*fyD[ind]
			}
		}
	}
	return
}
