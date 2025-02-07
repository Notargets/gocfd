package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"

	"github.com/stretchr/testify/assert"
)

func TestRTElementRT1Interpolation(t *testing.T) {
	// Verify the interpolation of a constant vector field onto the element
	P := 1
	rt := NewRTElement(P)
	rt.V.Print("V")
	rt.VInv.Print("VInv")
	f1, f2 := utils.NewVectorConstant(rt.Np, -1),
		utils.NewVectorConstant(rt.Np, 10)
	rt.ProjectFunctionOntoDOF(f1.DataP, f2.DataP)

	C := rt.VInv.Mul(rt.Projection)
	for i := 0; i < rt.Np; i++ {
		rt.Phi[i].Coefficient = C.At(i, 0)
	}

	// For each polynomial evaluation at (r,s)i
	f_rt_dot := make([]float64, rt.Np)
	for i := 0; i < rt.Np; i++ {
		r_i, s_i := rt.R.AtVec(i), rt.S.AtVec(i)
		b_i := rt.Phi[i].BasisVector.Eval(r_i, s_i)
		// Sum of the basis polynomials over j, each dotted with basis vector_i
		for j := 0; j < rt.Np; j++ {
			f_rt_dot[i] += rt.Phi[j].Dot(r_i, s_i, b_i)
		}
		r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		fmt.Printf("f_rt[%f,%f]=%f, f_proj=%f\n",
			r, s, f_rt_dot[i], rt.Projection.At(i, 0))
	}
	assert.InDeltaSlicef(t, rt.Projection.DataP, f_rt_dot, 0.000001,
		"Interpolation Check")
}

func TestRTElementErvinRT1(t *testing.T) {
	var (
		P      = 1
		Np     = (P + 1) * (P + 3)
		NpInt  = (P) * (P + 1) / 2
		NpEdge = P + 1
		g1     = 0.5 - math.Sqrt(3)/6
		g2     = 0.5 + math.Sqrt(3)/6
	)
	scalarMult := func(p float64, v [2]float64) (v2 [2]float64) {
		v2 = [2]float64{p * v[0], p * v[1]}
		return
	}
	conv := func(r float64) float64 { return (r + 1.) / 2. }
	l1 := func(t_rs float64) (val float64) {
		tt := conv(t_rs)
		val = (tt - g2) / (g1 - g2)
		return
	}
	l2 := func(t_rs float64) (val float64) {
		tt := conv(t_rs)
		val = (tt - g1) / (g2 - g1)
		return
	}
	e1 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64,
		div float64) {
		// Bottom edge
		var (
			xi, eta = conv(r), conv(s)
		)
		if len(derivO) > 0 {
			div = 0
			return
		}
		v[0] = xi - 0.5
		v[1] = eta - 1
		return
	}
	e2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64,
		div float64) {
		// Hypotenuse
		var (
			xi, eta = conv(r), conv(s)
			sr2     = math.Sqrt2
		)
		if len(derivO) > 0 {
			div = 0
			return
		}
		v[0] = sr2 * xi
		v[1] = sr2 * eta
		return
	}
	e3 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64,
		div float64) {
		// Left edge
		var (
			xi, eta = conv(r), conv(s)
		)
		if len(derivO) > 0 {
			div = 0
			return
		}
		v[0] = xi - 1
		v[1] = eta - 0.5
		return
	}
	e4 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64,
		div float64) {
		var (
			xi, eta = conv(r), conv(s)
		)
		if len(derivO) != 0 {
			div = (3.*s + 1.) / 4.
			return
		}
		v[0] = eta * xi
		v[1] = eta * (eta - 1)
		return
	}
	e5 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64,
		div float64) {
		var (
			xi, eta = conv(r), conv(s)
		)
		if len(derivO) != 0 {
			div = (3.*r + 1.) / 4.
			return
		}
		v[0] = xi * (xi - 1)
		v[1] = xi * eta
		return
	}
	phiInt1 := func(r, s float64, derivO ...DerivativeDirection) (
		v [2]float64, div float64) {
		v, div = e4(r, s)
		return
	}
	phiInt2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e5(r, s)
		return
	}
	phiEdge1_1 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e1(r, s)
		// TODO: calculate div contribution from l1
		v = scalarMult(l1(r), v)
		return
	}
	phiEdge1_2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e1(r, s)
		// TODO: calculate div contribution from l2
		v = scalarMult(l2(r), v)
		return
	}
	phiEdge2_1 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e2(r, s)
		// TODO: calculate div contribution from l1
		v = scalarMult(l1(s), v)
		return
	}
	phiEdge2_2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e2(r, s)
		// TODO: calculate div contribution from l2
		v = scalarMult(l2(s), v)
		return
	}
	phiEdge3_1 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e3(r, s)
		// TODO: calculate div contribution from l2
		v = scalarMult(l2(s), v)
		return
	}
	phiEdge3_2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e3(r, s)
		// TODO: calculate div contribution from l1
		v = scalarMult(l1(s), v)
		return
	}
	dot := func(v1, v2 [2]float64) (val float64) {
		val = v1[0]*v2[0] + v1[1]*v2[1]
		return
	}

	type phi func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64)

	phi_j := []phi{phiInt1, phiInt2, phiEdge1_1, phiEdge1_2, phiEdge2_1,
		phiEdge2_2, phiEdge3_1, phiEdge3_2}

	RInt, SInt := NodesEpsilon(P - 1)
	R, S := utils.NewVector(Np), utils.NewVector(Np)
	for i := 0; i < NpInt; i++ {
		R.DataP[i] = RInt.DataP[i]
		S.DataP[i] = SInt.DataP[i]
		R.DataP[i+NpInt] = RInt.DataP[i]
		S.DataP[i+NpInt] = SInt.DataP[i]
	}
	rconv := func(xi float64) (r float64) {
		r = 2*xi - 1
		return
	}
	i := 2 * NpInt
	R.DataP[i], S.DataP[i] = rconv(g1), rconv(0)
	i++
	R.DataP[i], S.DataP[i] = rconv(g2), rconv(0)
	i++
	R.DataP[i], S.DataP[i] = -rconv(g1), rconv(g1)
	i++
	R.DataP[i], S.DataP[i] = -rconv(g2), rconv(g2)
	i++
	R.DataP[i], S.DataP[i] = rconv(0), rconv(g2)
	i++
	R.DataP[i], S.DataP[i] = rconv(0), rconv(g1)

	edgeNum := func(j int) (eNum RTFunctionNumber) {
		switch {
		case j < NpInt:
			eNum = E4
		case j >= NpInt && j < 2*NpInt:
			eNum = E5
		case j >= 2*NpInt && j < 2*NpInt+NpEdge:
			eNum = E1
		case j >= 2*NpInt+NpEdge && j < 2*NpInt+2*NpEdge:
			eNum = E2
		case j >= 2*NpInt+2*NpEdge && j < 2*NpInt+3*NpEdge:
			eNum = E3
		}
		return
	}
	V := utils.NewMatrix(Np, Np)
	for i = 0; i < Np; i++ {
		r_i, s_i := R.DataP[i], S.DataP[i]
		var v_i [2]float64
		switch edgeNum(i) {
		case E4:
			v_i, _ = e4(r_i, s_i)
		case E5:
			v_i, _ = e5(r_i, s_i)
		case E1:
			v_i, _ = e1(r_i, s_i)
		case E2:
			v_i, _ = e2(r_i, s_i)
		case E3:
			v_i, _ = e3(r_i, s_i)
		}
		// fmt.Printf("v_%d[%f,%f] = [%f,%f]\n", i, r_i, s_i, v_i[0], v_i[1])
		for j := 0; j < Np; j++ {
			v_j, _ := phi_j[j](r_i, s_i)
			V.Set(i, j, dot(v_j, v_i))
		}
	}
	rt := NewRTElement(P)
	assert.InDeltaSlicef(t, V.DataP, rt.V.DataP, 0.000001,
		"Basis Matrix Comparison - Manual Test vs. Element")
}

func _TestRTElementPerformanceRT2(t *testing.T) {
	// We test RT2 in isolation because RT:
	// - uses only the analytic interior basis functions E4 and E5 times the
	//   analytic polynomial multiplier functions in Ervin
	var (
		dt DivTest
	)
	dt = PolyField{}

	fmt.Println("Begin Divergence Test")
	P := 2
	rt := NewRTElement(P)
	rt.Div.Print("Div")

	// First check a constant field - divergence should be zero
	// TODO: Interior points 1 and 2 of [0,1,2]E4 and 1 and 2 of [0,1,2]E5 show
	// TODO: non zero divergence
	Np := rt.Np
	divFcalc := make([]float64, Np)
	s1, s2 := make([]float64, Np), make([]float64, Np)
	for i := 0; i < Np; i++ {
		r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		f1, f2 := dt.F(r, s, 0)
		// f1, f2 := dt.F(r, s, P-1)
		s1[i], s2[i] = f1, f2
		dF := dt.divF(r, s, 0)
		// dF := dt.divF(r, s, P-1)
		divFcalc[i] = dF
		// fmt.Printf("F[%f,%f]=[%f,%f], divF[%f,%f]=%f\n", r, s, f1, f2, r, s, dF)
	}
	dFcalc := utils.NewMatrix(Np, 1, divFcalc)
	dFcalc.Transpose().Print("Reference Div")
	rt.ProjectFunctionOntoDOF(s1, s2)
	dB := rt.Projection
	rt.Div.Mul(dB).Transpose().Print("Calculated Divergence")
}

func TestRTElementDivergenceRT1(t *testing.T) {
	// We test RT1 first in isolation because RT1:
	// - uses only the analytic interior basis functions E4 and E5
	// - uses the Lagrange 1D polynomial on edges
	// - is the simplest construction to test divergence
	var (
		dt DivTest
	)
	dt = PolyField{}

	fmt.Println("Begin Divergence Test")
	P := 1
	rt := NewRTElement(P)
	rt.Div.Print("Div")

	Np := rt.Np
	divFcalc := make([]float64, Np)
	s1, s2 := make([]float64, Np), make([]float64, Np)
	for i := 0; i < Np; i++ {
		r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		f1, f2 := dt.F(r, s, P-1)
		s1[i], s2[i] = f1, f2
		dF := dt.divF(r, s, P-1)
		divFcalc[i] = dF
	}
	dFcalc := utils.NewMatrix(Np, 1, divFcalc)
	dFcalc.Transpose().Print("Reference Div")
	rt.ProjectFunctionOntoDOF(s1, s2)
	dB := rt.Projection
	rt.Div.Mul(dB).Transpose().Print("Calculated Divergence")
}

func TestRTElementVerifyErvinRT1(t *testing.T) {
	// Check the basis vectors for RT1 against Ervin's RT1 construction
	// Altered to check against edge normal vectors instead of E1,E2,E3 in Ervin
	P := 1
	rt := NewRTElement(P)
	BasisVector := [2]utils.Matrix{utils.NewMatrix(rt.Np, 1),
		utils.NewMatrix(rt.Np, 1)}
	var v [2]float64
	fmt.Printf("RT%d Element Basis\n", P)
	fmt.Printf("Np:%d; NpInt:%d; NpEdge:%d\n", rt.Np, rt.NpInt, rt.NpEdge)
	for j := 0; j < rt.Np; j++ {
		r, s := rt.R.AtVec(j), rt.S.AtVec(j)
		v = rt.Phi[j].BasisVector.Eval(r, s)
		BasisVector[0].Set(j, 0, v[0])
		BasisVector[1].Set(j, 0, v[1])
	}
	g1, g2 := 0.5-math.Sqrt(3)/6, 0.5+math.Sqrt(3)/6
	l1 := func(t float64) (ll float64) {
		ll = (t - g2) / (g1 - g2)
		return
	}
	l2 := func(t float64) (ll float64) {
		ll = (t - g1) / (g2 - g1)
		return
	}
	e1 := func(xi, eta float64) (v [2]float64) {
		var (
			s2 = math.Sqrt(2)
		)
		v = [2]float64{s2 * xi, s2 * eta} // Ervin
		return
	}
	e2 := func(xi, eta float64) (v [2]float64) {
		v = [2]float64{xi - 1, eta - 0.5} // Ervin
		return
	}
	e3 := func(xi, eta float64) (v [2]float64) {
		v = [2]float64{xi - 0.5, eta - 1} // Ervin
		return
	}
	e4 := func(xi, eta float64) (v [2]float64) {
		v = [2]float64{eta * xi, eta * (eta - 1)}
		return
	}
	e5 := func(xi, eta float64) (v [2]float64) {
		v = [2]float64{xi * (xi - 1), xi * eta}
		return
	}
	edge1Basis := func(xi, eta float64, j int) (v [2]float64) {
		// Hypotenuse
		var (
			l float64
		)
		e := e1(xi, eta)
		switch j {
		case 0:
			l = l1(eta)
		case 1:
			l = l2(eta)
		}
		v = [2]float64{l * e[0], l * e[1]}
		return
	}
	edge2Basis := func(xi, eta float64, j int) (v [2]float64) {
		// Left Edge
		var (
			l float64
		)
		e := e2(xi, eta)
		switch j {
		case 0:
			l = l2(eta)
		case 1:
			l = l1(eta)
		}
		v = [2]float64{l * e[0], l * e[1]}
		return
	}
	edge3Basis := func(xi, eta float64, j int) (v [2]float64) {
		// Bottom Edge
		var (
			l float64
		)
		e := e3(xi, eta)
		switch j {
		case 0:
			l = l1(xi)
		case 1:
			l = l2(xi)
		}
		v = [2]float64{l * e[0], l * e[1]}
		return
	}
	convTo := func(r float64) (x float64) {
		// xi = (r+1)/2, eta = (s+1)/2
		// r = 2*xi -1, s = 2*eta -1
		x = (r + 1.) / 2.
		return
	}
	i := 0
	r, s := rt.R.AtVec(i), rt.S.AtVec(i)
	v1 := e4(convTo(r), convTo(s))
	i = 1
	r, s = rt.R.AtVec(i), rt.S.AtVec(i)
	v2 := e5(convTo(r), convTo(s))
	i = 2
	r, s = rt.R.AtVec(i), rt.S.AtVec(i)
	v3 := edge3Basis(convTo(r), convTo(s), 0)
	i = 3
	r, s = rt.R.AtVec(i), rt.S.AtVec(i)
	v4 := edge3Basis(convTo(r), convTo(s), 1)
	i = 4
	r, s = rt.R.AtVec(i), rt.S.AtVec(i)
	v5 := edge1Basis(convTo(r), convTo(s), 0)
	i = 5
	r, s = rt.R.AtVec(i), rt.S.AtVec(i)
	v6 := edge1Basis(convTo(r), convTo(s), 1)
	i = 6
	r, s = rt.R.AtVec(i), rt.S.AtVec(i)
	v7 := edge2Basis(convTo(r), convTo(s), 0)
	i = 7
	r, s = rt.R.AtVec(i), rt.S.AtVec(i)
	v8 := edge2Basis(convTo(r), convTo(s), 1)

	i = 0
	BV0E := utils.NewVector(8, []float64{
		v1[i], v2[i], v3[i], v4[i], v5[i], v6[i], v7[i], v8[i],
	})
	i = 1
	BV1E := utils.NewVector(8, []float64{
		v1[i], v2[i], v3[i], v4[i], v5[i], v6[i], v7[i], v8[i],
	})
	// BV0E.Print("BV0E")
	// BV1E.Print("BV1E")
	for i := 0; i < 8; i++ {
		r, s = rt.R.AtVec(i), rt.S.AtVec(i)
		fmt.Printf("%d;[r,s]=[%5.2f,%5.2f] "+
			"[xi,eta]=[%5.2f,%5.2f] h:[%5.2f,%5.2f] Ervin:[%5.2f,%5.2f]\n",
			i+1, r, s, convTo(r), convTo(s),
			BasisVector[0].At(i, 0), BasisVector[1].At(i, 0),
			BV0E.AtVec(i), BV1E.AtVec(i))
		b_i := rt.Phi[i].BasisVector.Eval(r, s)
		assert.InDeltaf(t, b_i[0], BV0E.AtVec(i), 0.0001, "")
		assert.InDeltaf(t, b_i[1], BV1E.AtVec(i), 0.0001, "")
	}
}

func TestRTElement(t *testing.T) {
	{
		// Check term-wise orthogonal 2D polynomial basis
		N := 2
		R, S := NodesEpsilon(N - 1)
		JB2D := NewJacobiBasis2D(N-1, R, S, 0, 0)
		ii, jj := 1, 1
		p := JB2D.Simplex2DP(R, S, ii, jj)
		ddr, dds := JB2D.GradSimplex2DP(R, S, ii, jj)
		Np := R.Len()
		pCheck, ddrCheck, ddsCheck := make([]float64, Np), make([]float64, Np), make([]float64, Np)
		for i, rVal := range R.DataP {
			sVal := S.DataP[i]
			ddrCheck[i] = JB2D.PolynomialTermDr(rVal, sVal, ii, jj)
			ddsCheck[i] = JB2D.PolynomialTermDs(rVal, sVal, ii, jj)
			pCheck[i] = JB2D.PolynomialTerm(rVal, sVal, ii, jj)
		}
		assert.True(t, nearVec(pCheck, p, 0.000001))
		assert.True(t, nearVec(ddrCheck, ddr, 0.000001))
		assert.True(t, nearVec(ddsCheck, dds, 0.000001))
	}
	errorCheck := func(N int, div, divCheck []float64) (minInt, maxInt, minEdge, maxEdge float64) {
		var (
			Npm    = len(div)
			errors = make([]float64, Npm)
		)
		for i := 0; i < Npm; i++ {
			// var ddr, dds float64
			errors[i] = div[i] - divCheck[i]
		}
		minInt, maxInt = errors[0], errors[0]
		Nint := N * (N + 1) / 2
		minEdge, maxEdge = errors[Nint], errors[Nint]
		for i := 0; i < Nint; i++ {
			errAbs := math.Abs(errors[i])
			if minInt > errAbs {
				minInt = errAbs
			}
			if maxInt < errAbs {
				maxInt = errAbs
			}
		}
		for i := Nint; i < Npm; i++ {
			errAbs := math.Abs(errors[i])
			if minEdge > errAbs {
				minEdge = errAbs
			}
			if maxEdge < errAbs {
				maxEdge = errAbs
			}
		}
		fmt.Printf("Order = %d, ", N)
		fmt.Printf("Min, Max Int Err = %8.5f, %8.5f, Min, Max Edge Err = %8.5f, %8.5f\n", minInt, maxInt, minEdge, maxEdge)
		return
	}
	checkSolution := func(rt *RTElement, Order int) (s1, s2, divCheck []float64) {
		var (
			Np = rt.Np
		)
		s1, s2 = make([]float64, Np), make([]float64, Np)
		divCheck = make([]float64, Np)
		var ss1, ss2 float64
		for i := 0; i < Np; i++ {
			r := rt.R.DataP[i]
			s := rt.S.DataP[i]
			ccf := float64(Order)
			s1[i] = utils.POW(r, Order)
			s2[i] = utils.POW(s, Order)
			ss1, ss2 = ccf*utils.POW(r, Order-1), ccf*utils.POW(s, Order-1)
			divCheck[i] = ss1 + ss2
		}
		return
	}
	// against analytical solution
	// Nend := 8
	// for N := 1; N < Nend; N++ {
	N := 1
	rt := NewRTElement(N)
	for cOrder := 0; cOrder < N; cOrder++ {
		fmt.Printf("Check Order = %d, ", cOrder)
		// [s1,s2] values for each location in {R,S}
		s1, s2, divCheck := checkSolution(rt, cOrder)
		rt.ProjectFunctionOntoDOF(s1, s2)
		divM := rt.Div.Mul(rt.Projection)
		// fmt.Println(divM.Print("divM"))
		minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, divM.DataP, divCheck)
		assert.True(t, near(minerrInt, 0.0, 0.00001))
		assert.True(t, near(maxerrInt, 0.0, 0.00001))
		assert.True(t, near(minerrEdge, 0.0, 0.00001))
		assert.True(t, near(maxerrEdge, 0.0, 0.00001))
	}
	// }
	plot := false
	if plot {
		N := 2
		rt := NewRTElement(N)
		s1, s2 := make([]float64, rt.R.Len()), make([]float64, rt.R.Len())
		for i := range rt.R.DataP {
			s1[i] = 1
			s2[i] = 1
		}
		if plot {
			chart := PlotTestTri(true)
			points := utils.ArraysToPoints(rt.R.DataP, rt.S.DataP)
			f := utils.ArraysTo2Vector(s1, s2, 0.1)
			_ = chart.AddVectors("test function", points, f, chart2d.Solid, utils.GetColor(utils.Green))
			utils.SleepFor(500000)
		}
	}
}

func PlotTestTri(plotGeom bool) (chart *chart2d.Chart2D) {
	var (
		points  []graphics2D.Point
		trimesh graphics2D.TriMesh
		K       = 1
	)

	points = make([]graphics2D.Point, 3)
	points[0].X[0], points[0].X[1] = -1, -1
	points[1].X[0], points[1].X[1] = 1, -1
	points[2].X[0], points[2].X[1] = -1, 1

	trimesh.Triangles = make([]graphics2D.Triangle, K)
	colorMap := utils2.NewColorMap(0, 1, 1)
	trimesh.Triangles[0].Nodes[0] = 0
	trimesh.Triangles[0].Nodes[1] = 1
	trimesh.Triangles[0].Nodes[2] = 2
	trimesh.Geometry = points
	box := graphics2D.NewBoundingBox(trimesh.GetGeometry())
	chart = chart2d.NewChart2D(1024, 1024, box.XMin[0], box.XMax[0], box.XMin[1], box.XMax[1])
	chart.AddColorMap(colorMap)
	go chart.Plot()

	if plotGeom {
		if err := chart.AddTriMesh("TriMesh", trimesh,
			chart2d.CrossGlyph, 0.1, chart2d.Solid,
			utils.GetColor(utils.Black)); err != nil {
			panic("unable to add graph series")
		}
	}
	return
}

func checkIfUnitMatrix(t *testing.T, A utils.Matrix) (isDiag bool) {
	var (
		Np, _ = A.Dims()
	)
	for j := 0; j < Np; j++ {
		for i := 0; i < Np; i++ {
			if i == j {
				assert.InDeltaf(t, 1., A.At(i, j), 0.00001, "")
			} else {
				assert.InDeltaf(t, 0., A.At(i, j), 0.00001, "")
			}
		}
	}
	return
}

type DivTest interface {
	F(r, s float64, P int) (f1, f2 float64)
	divF(r, s float64, P int) (div float64)
}

type SinCosField struct{}

func (scf SinCosField) F(r, s float64, P int) (f1, f2 float64) {
	var (
		Pi = math.Pi
	)
	conv := func(r float64) (xi float64) {
		xi = Pi * (r + 1)
		return
	}
	f1, f2 = math.Sin(conv(r)), math.Cos(conv(s))
	return
}

func (scf SinCosField) divF(r, s float64, P int) (div float64) {
	var (
		Pi = math.Pi
	)
	conv := func(r float64) (xi float64) {
		xi = Pi * (r + 1)
		return
	}
	div = (math.Cos(conv(r)) - math.Sin(conv(s)))
	div += Pi * (math.Sin(conv(r)) + math.Cos(conv(s)))
	return
}

type PolyField struct{}

func (lpf PolyField) F(r, s float64, P int) (f1, f2 float64) {
	var (
		Pi = math.Pi
		p  = float64(P)
	)
	conv := func(r float64) (xi float64) {
		xi = Pi * (r + 1)
		return
	}
	f1, f2 = math.Pow(conv(r), p), math.Pow(conv(s), p)
	return
}

func (lpf PolyField) divF(r, s float64, P int) (div float64) {
	var (
		Pi = math.Pi
		p  = float64(P)
	)
	conv := func(r float64) (xi float64) {
		xi = Pi * (r + 1)
		return
	}
	// val1 = (Pi*(r+1))^p
	// dval1/dr = Pi*p*((Pi*(r+1))^(p-1))
	// val2 = (Pi*(s+1))^p
	// dval2/ds = Pi*p*((Pi*(s+1))^(p-1))
	// div = p*Pi*((Pi*(r+1))^(p-1) + (Pi*(s+1))^(p-1))
	if P > 0 {
		div = p * Pi * (math.Pow(conv(r), p-1) + math.Pow(conv(s), p-1))
	} else {
		div = 0
	}
	return
}
