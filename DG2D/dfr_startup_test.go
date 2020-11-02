package DG2D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"
	"github.com/stretchr/testify/assert"
)

func TestDFR2D(t *testing.T) {
	{ // Test Interpolation
		N := 2
		dfr := NewDFR2D(N)
		el := dfr.SolutionElement
		s := make([]float64, el.Np)
		for i := 0; i < el.Np; i++ {
			s[i] = float64(2 * i)
		}
		// For each nodal location, interpolate a value (should equal the nodal function value)
		// Build an interpolating polynomial matrix using the nodal geometry
		interpM := el.Simplex2DInterpolatingPolyMatrix(el.R, el.S)
		values := interpM.Mul(utils.NewMatrix(el.Np, 1, s))
		assert.True(t, nearVec(s, values.Data(), 0.0000001))
	}
	{ // Test point distribution
		N := 1
		dfr := NewDFR2D(N)
		el := dfr.SolutionElement
		rt := dfr.FluxElement
		assert.True(t, nearVec(rt.GetInternalLocations(rt.R), el.R.Data(), 0.000001))
		assert.True(t, nearVec(rt.GetInternalLocations(rt.S), el.S.Data(), 0.000001))
		//fmt.Printf("Edge R = %v\n", rt.GetEdgeLocations(rt.R))
		//fmt.Printf("Edge S = %v\n", rt.GetEdgeLocations(rt.S))
	}
	{ // Test interpolation from solution points to flux points
		N := 1
		dfr := NewDFR2D(N)
		el := dfr.SolutionElement
		rt := dfr.FluxElement
		assert.Equal(t, el.Np, rt.Nint)
		solution := make([]float64, rt.Nint)
		// Load some values into the solution space
		for i := 0; i < rt.Nint; i++ {
			solution[i] = float64(i) / float64(rt.Nint-1)
		}
		// Interpolate from interior to flux points
		sV := utils.NewMatrix(rt.Nint, 1, solution)
		_ = dfr.FluxInterpMatrix.Mul(sV)
		//fmt.Printf("%s\n", fluxInterp.Print("fluxInterp"))
		//fmt.Printf("%s\n", sV.Print("sV"))
	}
	{ // Test packed int for edge labeling
		en := NewEdgeNumber([2]int{1, 0})
		assert.Equal(t, EdgeNumber(1<<32), en)
		assert.Equal(t, [2]int{0, 1}, en.GetVertices(false))

		en = NewEdgeNumber([2]int{0, 1})
		assert.Equal(t, EdgeNumber(1<<32), en)
		assert.Equal(t, [2]int{0, 1}, en.GetVertices(false))

		en = NewEdgeNumber([2]int{0, 10})
		assert.Equal(t, EdgeNumber(10*(1<<32)), en)
		assert.Equal(t, [2]int{0, 10}, en.GetVertices(false))

		en = NewEdgeNumber([2]int{100, 0})
		assert.Equal(t, EdgeNumber(100*(1<<32)), en)
		assert.Equal(t, [2]int{0, 100}, en.GetVertices(false))

		en = NewEdgeNumber([2]int{100, 1})
		assert.Equal(t, EdgeNumber(100*(1<<32)+1), en)
		assert.Equal(t, [2]int{1, 100}, en.GetVertices(false))

		en = NewEdgeNumber([2]int{100, 100001})
		assert.Equal(t, EdgeNumber(100001*(1<<32)+100), en)
		assert.Equal(t, [2]int{100, 100001}, en.GetVertices(false))

		// Test maximum/minimum indices
		en = NewEdgeNumber([2]int{1, 1<<32 - 1})
		assert.Equal(t, EdgeNumber((1<<32-1)<<32+1), en)
		assert.Equal(t, [2]int{1, 1<<32 - 1}, en.GetVertices(false))

		en = NewEdgeNumber([2]int{1<<32 - 1, 1<<32 - 1})
		assert.Equal(t, EdgeNumber(1<<64-1), en)
		assert.Equal(t, [2]int{1<<32 - 1, 1<<32 - 1}, en.GetVertices(false))

		en = NewEdgeNumber([2]int{1<<32 - 1, 1})
		assert.Equal(t, EdgeNumber((1<<32-1)<<32+1), en)
		assert.Equal(t, [2]int{1, 1<<32 - 1}, en.GetVertices(false))
	}
	{ // Test triangulation
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
		dfr := NewDFR2D(N, "test_tris_5.neu")
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
				en := NewEdgeNumber([2]int{0, 2}) // Edge formed by vertices 0 and 2
				e := trn.Edges[en]
				assert.Equal(t, uint8(2), e.NumConnectedTris)
				assert.Equal(t, uint32(0), e.ConnectedTris[0])
				assert.Equal(t, InternalEdgeNumber(2), e.ConnectedTriEdgeNumber[First])
				assert.Equal(t, uint32(1), e.ConnectedTris[1])
				assert.Equal(t, InternalEdgeNumber(0), e.ConnectedTriEdgeNumber[Second])
			}
			{ // 0-1, 2-3: 	BC Inflow
				en := NewEdgeNumber([2]int{0, 1})
				e := trn.Edges[en]
				assert.Equal(t, BC_In, e.BCType)
				en = NewEdgeNumber([2]int{2, 3})
				e = trn.Edges[en]
				assert.Equal(t, BC_In, e.BCType)
			}
		}
		// Test Piola transform and jacobian
		for en, e := range dfr.Tris.Edges {
			assert.True(t, e.NumConnectedTris > 0)
			for conn := 0; conn < int(e.NumConnectedTris); conn++ {
				k := e.ConnectedTris[conn]
				revDir := bool(e.ConnectedTriDirection[conn])
				Jt := transpose(dfr.J.Row(int(k)).Data()[0:4])       // Transpose is applied to normals, which are cross products of vectors
				JtInv := transpose(dfr.Jinv.Row(int(k)).Data()[0:4]) // Transpose is applied to normals, which are cross products of vectors
				Jdet := math.Abs(dfr.Jdet.Row(int(k)).Data()[0])
				x1, x2 := dfr.Tris.GetEdgeCoordinates(en, revDir, dfr.VX, dfr.VY)
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
					assert.True(t, near(normal[0], nxTT[0], 0.00001))
					assert.True(t, near(normal[1], nxTT[1], 0.00001))
				}
				{ // Check scaling factor ||n||, used in transforming face normals
					normal = normalize(normal)
					normal[0] *= e.IInII[conn] // scale ||n|| the vector coords on the edge prior to transform
					normal[1] *= e.IInII[conn] // scale ||n|| the vector coords on the edge prior to transform
					nxT := multiply(Jt, Jdet, normal)
					oosr2 := 1. / math.Sqrt(2)
					switch edgeNumber { // The transformed and scaled normal should be unit for each edge direction
					case First:
						assert.True(t, near(0, nxT[0], 0.000001))
						assert.True(t, near(-1, nxT[1], 0.000001))
					case Second:
						assert.True(t, near(oosr2, nxT[0], 0.000001))
						assert.True(t, near(oosr2, nxT[1], 0.000001))
					case Third:
						assert.True(t, near(-1, nxT[0], 0.000001))
						assert.True(t, near(0, nxT[1], 0.000001))
					}
				}
			}
		}
	}
	{ // Test divergence
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
					xD, yD   = dfr.FluxX.Col(k).Data()[0:Np], dfr.FluxY.Col(k).Data()[0:Np]
					fxD, fyD = Fx.Data(), Fy.Data()
					dcD      = divCheck.Data()
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
		{ // Check Divergence for polynomial vector fields of order < N against analytical solution
			N := 7
			dfr := NewDFR2D(N, "test_tris_5.neu")
			rt := dfr.FluxElement
			for cOrder := 1; cOrder <= N; cOrder++ {
				//fmt.Printf("Check Order = %d, \n", cOrder)
				Fx, Fy, divCheck := checkSolution(dfr, cOrder)
				Fp := dfr.ProjectFluxOntoRTSpace(Fx, Fy)
				for k := 0; k < dfr.K; k++ {
					var (
						Fpk  = Fp.Row(k).ToMatrix()
						Jdet = dfr.Jdet.Row(k).Data()[0]
					)
					divM := rt.Div.Mul(Fpk).Scale(1. / Jdet)
					assert.True(t, nearVec(divCheck.Row(k).Data(), divM.Data(), 0.00001))
				}
			}
		}
		{ // Use transformed edge normal projected flux in place of flux and check divergence again
			N := 7
			dfr := NewDFR2D(N, "test_tris_5.neu")
			rt := dfr.FluxElement
			for cOrder := 1; cOrder <= N; cOrder++ {
				//fmt.Printf("Check Order = %d, \n", cOrder)
				Fx, Fy, divCheck := checkSolution(dfr, cOrder)
				Fp := dfr.ProjectFluxOntoRTSpace(Fx, Fy)
				SetNormalFluxOnEdges(dfr, Fx, Fy, &Fp)
				for k := 0; k < dfr.K; k++ {
					var (
						Fpk  = Fp.Row(k).ToMatrix()
						Jdet = dfr.Jdet.Row(k).Data()[0]
					)
					divM := rt.Div.Mul(Fpk).Scale(1. / Jdet)
					assert.True(t, nearVec(divCheck.Row(k).Data(), divM.Data(), 0.00001))
					//_, _ = divM, divCheck
				}
			}
		}
	}
}

func SetNormalFluxOnEdges(dfr *DFR2D, Fx, Fy utils.Matrix, Fp *utils.Matrix) {
	var (
		Np       = dfr.FluxElement.Np
		fxD, fyD = Fx.Data(), Fy.Data()
		fpD      = Fp.Data()
		Nint     = dfr.FluxElement.Nint
		Nedge    = dfr.FluxElement.Nedge
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
			x1, x2 := dfr.Tris.GetEdgeCoordinates(en, revDir, dfr.VX, dfr.VY)
			dx := [2]float64{x2[0] - x1[0], x2[1] - x1[1]}
			normal := normalize([2]float64{-dx[1], dx[0]})
			normal[0] *= e.IInII[conn]
			normal[1] *= e.IInII[conn]
			var edgeStart, edgeEnd int
			edgeNumber := e.ConnectedTriEdgeNumber[conn]
			switch edgeNumber {
			case Second: // The hypotenuse is stored first
				edgeStart, edgeEnd = 2*Nint, 2*Nint+Nedge
			case Third:
				edgeStart, edgeEnd = 2*Nint+Nedge, 2*Nint+2*Nedge
			case First:
				edgeStart, edgeEnd = 2*Nint+2*Nedge, 2*Nint+3*Nedge
			}
			//fmt.Printf("k, edgeNum, edgeStart, edgeEnd, IInII = %d, %s, %d, %d, %8.5f\n",
			//	k, edgeNumber.String(), edgeStart, edgeEnd, e.IInII[conn])
			for n := edgeStart; n < edgeEnd; n++ {
				ind := n + k*Np
				fpD[ind] = normal[0]*fxD[ind] + normal[1]*fyD[ind]
			}
		}
	}
	return
}
