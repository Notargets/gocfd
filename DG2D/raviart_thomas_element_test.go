package DG2D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/notargets/avs/chart2d"

	"github.com/stretchr/testify/assert"
)

func TestRTElementDivergence2(t *testing.T) {
	// We test RT1 first in isolation because RT1:
	// - uses only the analytic interior basis functions E4 and E5
	// - uses the Lagrange 1D polynomial on edges
	// - is the simplest construction to test divergence
	var (
		dt VectorTestField
	)
	dt = SinCosVectorField{}

	t.Log("Begin Divergence Test")
	PStart := 1
	PEnd := 2
	for P := PStart; P <= PEnd; P++ {
		t.Logf("---------------------------------------------\n")
		t.Logf("Checking Divergence for RT%d\n", P)
		t.Logf("---------------------------------------------\n")
		rt := NewRTElement(P)
		Np := rt.Np
		divFcalc := make([]float64, Np)
		s1, s2 := make([]float64, Np), make([]float64, Np)
		// for PField := 0; PField <= (P - 1); PField++ {
		t.Logf("\nReference Vector Field Sin/Cos\n")
		t.Logf("-------------------------------\n")
		for i := 0; i < Np; i++ {
			r, s := rt.R.AtVec(i), rt.S.AtVec(i)
			f1, f2 := dt.F(r, s, 0)
			s1[i], s2[i] = f1, f2
			dF := dt.Divergence(r, s, 0)
			divFcalc[i] = dF
		}
		dFReference := utils.NewMatrix(Np, 1, divFcalc)
		if testing.Verbose() {
			dFReference.Transpose().Print("Reference Div")
		}
		rt.ProjectFunctionOntoDOF(s1, s2)
		dB := rt.Projection
		calcDiv := rt.Div.Mul(dB)
		if testing.Verbose() {
			calcDiv.Transpose().Print("Calculated Divergence")
		}
		var err float64
		for i := 0; i < Np; i++ {
			err += math.Pow(calcDiv.At(i, 0)-dFReference.At(i, 0), 2)
		}
		rms := math.Sqrt(err / float64(Np))
		t.Logf("RMS Err = %f\n", rms)
		// assert.InDeltaSlice(t, dFReference.DataP, calcDiv.DataP, 0.0001)
	}
}

func TestRTElementRTInterpolation(t *testing.T) {
	// Verify the interpolation of a constant vector field onto the element
	PStart := 1
	PEnd := 6
	for P := PStart; P <= PEnd; P++ {
		var (
			dt VectorTestField
		)
		dt = PolyVectorField{}

		rt := NewRTElement(P)
		if testing.Verbose() {
			rt.V.Print("V")
			rt.VInv.Print("VInv")
		}
		s1, s2 := make([]float64, rt.Np), make([]float64, rt.Np)
		for PField := 0; PField <= P; PField++ {
			for i := 0; i < rt.Np; i++ {
				r, s := rt.R.AtVec(i), rt.S.AtVec(i)
				f1, f2 := dt.F(r, s, PField)
				s1[i], s2[i] = f1, f2
			}
			rt.ProjectFunctionOntoDOF(s1, s2)

			C := rt.VInv.Mul(rt.Projection)

			// For each polynomial evaluation at (r,s)i
			f_rt_dot := make([]float64, rt.Np)
			for i := 0; i < rt.Np; i++ {
				r_i, s_i := rt.R.AtVec(i), rt.S.AtVec(i)
				b_i := rt.Phi[i].BasisVector.Eval(r_i, s_i)
				// Sum of the basis polynomials over j, each dotted with basis vector_i
				for j := 0; j < rt.Np; j++ {
					f_rt_dot[i] += rt.Phi[j].Dot(r_i, s_i, b_i) * C.At(j, 0)
				}
				if PField >= P+1 {
					r, s := rt.R.AtVec(i), rt.S.AtVec(i)
					t.Logf("f_rt[%f,%f]=%f, f_proj=%f\n",
						r, s, f_rt_dot[i], rt.Projection.At(i, 0))
				}
			}
			assert.InDeltaSlicef(t, rt.Projection.DataP, f_rt_dot, 0.000001,
				"Interpolation Check")
		}
	}
}

func TestRTElementDivergence(t *testing.T) {
	// We test RT1 first in isolation because RT1:
	// - uses only the analytic interior basis functions E4 and E5
	// - uses the Lagrange 1D polynomial on edges
	// - is the simplest construction to test divergence
	var (
		dt VectorTestField
	)
	dt = PolyVectorField{}

	t.Log("Begin Divergence Test")
	// P := 1
	PStart := 1
	PEnd := 3
	for P := PStart; P <= PEnd; P++ {
		PFieldStart := 0
		PFieldEnd := P
		t.Logf("---------------------------------------------\n")
		t.Logf("Checking Divergence for RT%d\n", P)
		t.Logf("---------------------------------------------\n")
		rt := NewRTElement(P)
		Np := rt.Np
		divFcalc := make([]float64, Np)
		s1, s2 := make([]float64, Np), make([]float64, Np)
		// for PField := 0; PField <= (P - 1); PField++ {
		for PField := PFieldStart; PField <= PFieldEnd; PField++ {
			t.Logf("\nReference Vector Field Order:%d\n", PField)
			t.Logf("-------------------------------\n")
			for i := 0; i < Np; i++ {
				r, s := rt.R.AtVec(i), rt.S.AtVec(i)
				f1, f2 := dt.F(r, s, PField)
				s1[i], s2[i] = f1, f2
				dF := dt.Divergence(r, s, PField)
				divFcalc[i] = dF
			}
			dFReference := utils.NewMatrix(Np, 1, divFcalc)
			if testing.Verbose() {
				dFReference.Transpose().Print("Reference Div")
			}
			rt.ProjectFunctionOntoDOF(s1, s2)
			dB := rt.Projection
			// dB.Transpose().Print("F Projection")
			// rt.VInv.Mul(dB).Print("Coefficients")
			calcDiv := rt.Div.Mul(dB)
			if testing.Verbose() {
				calcDiv.Transpose().Print("Calculated Divergence")
				// calcCoeffs := rt.V.Mul(dB)
				// dB.Transpose().Print("Projected Field")
				// calcCoeffs.Transpose().Print("Calculated Coeffs")
			}
			assert.InDeltaSlice(t, dFReference.DataP, calcDiv.DataP, 0.0001)
		}
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
		assert.InDeltaSlicef(t, pCheck, p, 0.000001, "")
		assert.InDeltaSlicef(t, ddrCheck, ddr, 0.000001, "")
		assert.InDeltaSlicef(t, ddsCheck, dds, 0.000001, "")
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
		t.Logf("Order = %d, ", N)
		t.Logf("Min, Max Int Err = %8.5f, %8.5f, Min, Max Edge Err = %8.5f, %8.5f\n", minInt, maxInt, minEdge, maxEdge)
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
		t.Logf("Check Order = %d, ", cOrder)
		// [s1,s2] values for each location in {R,S}
		s1, s2, divCheck := checkSolution(rt, cOrder)
		rt.ProjectFunctionOntoDOF(s1, s2)
		divM := rt.Div.Mul(rt.Projection)
		// fmt.Println(divM.Print("divM"))
		minerrInt, maxerrInt, minerrEdge, maxerrEdge := errorCheck(N, divM.DataP, divCheck)
		assert.InDeltaf(t, minerrInt, 0.0, 0.00001, "")
		assert.InDeltaf(t, maxerrInt, 0.0, 0.00001, "")
		assert.InDeltaf(t, minerrEdge, 0.0, 0.00001, "")
		assert.InDeltaf(t, maxerrEdge, 0.0, 0.00001, "")
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
