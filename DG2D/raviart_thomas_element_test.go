package DG2D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/utils"

	"github.com/stretchr/testify/assert"
)

func TestRTElement(t *testing.T) {
	// for _, rtb := range []RTBasisType{SimplexRTBasis} {
	// for _, rtb := range []RTBasisType{ErvinBasisRT, SimplexRTBasis} {
	for _, rtb := range []RTBasisType{ErvinBasisRT} {
		var PMax int
		switch rtb {
		case ErvinBasisRT:
			PMax = 2
		case SimplexRTBasis:
			PMax = 2
		}
		InterpolationFunctions_Test(t, PMax)
		// t.Logf("Testing RT Interpolation for %v\n", rtb.String())
		// RTInterpolation_Test(t, rtb, PMax)

		// t.Logf("Testing RT divergence on SinCos Fields for %v\n",
		// 	rtb.String())
		// RTDivergenceSinCos_Test(t, rtb, PMax)
	}
}

func InterpolationFunctions_Test(t *testing.T, PMax int) {
	var (
		dt  VectorTestField
		tol = 0.0000001
	)
	dt = PolyVectorField{}

	t.Log("Begin Interpolation Function Test")
	// P := 1
	PStart := 1
	PEnd := PMax
	for P := PStart; P <= PEnd; P++ {
		PFieldStart := 0
		PFieldEnd := P
		t.Logf("---------------------------------------------\n")
		t.Logf("Checking Interpolation Functions for RT%d\n", P)
		t.Logf("---------------------------------------------\n")
		rt := NewRTElement(P, ErvinBasisRT)
		Np := rt.Np
		A := utils.NewMatrix(Np, Np)
		f1, f2 := make([]float64, Np), make([]float64, Np)
		divF := make([]float64, Np)
		for PField := PFieldStart; PField <= PFieldEnd; PField++ {
			t.Logf("\nReference Vector Field Order:%d\n", PField)
			t.Logf("-------------------------------\n")
			for i := 0; i < Np; i++ {
				r, s := rt.R.AtVec(i), rt.S.AtVec(i)
				f1[i], f2[i] = dt.F(r, s, PField)
				divF[i] = dt.Divergence(r, s, PField)
			}
			rt.ProjectFunctionOntoDOF(f1, f2)
			for j := 0; j < Np; j++ {
				for i := 0; i < Np; i++ {
					r, s := rt.R.AtVec(i), rt.S.AtVec(i)
					b_i := rt.DOFVectors[i].Eval(0, 0)
					dot := rt.Psi[j].dot(r, s, b_i)
					A.Set(i, j, dot)
					tol := 0.0000001
					if i == j {
						assert.InDeltaf(t, 1., dot, tol, "")
					} else {
						assert.InDeltaf(t, 0., dot, tol, "")
					}
				}
			}
			for i := 0; i < Np; i++ {
				r, s := rt.R.AtVec(i), rt.S.AtVec(i)
				var div float64
				for j := 0; j < Np; j++ {
					div += rt.Projection.At(j, 0) * rt.Psi[j].Divergence(r, s)
				}
				// t.Logf("Div[%f,%f] = %f(Ref), %f(Calc)\n",
				// 	r, s, divF[i], div)
				assert.InDeltaf(t, divF[i], div, tol, "")
			}
		}
		// if testing.Verbose() {
		// 	A.Print("Should Be Identity")
		// }
		assert.True(t, isIdentityMatrix(A, tol))
	}
}

func RTDivergenceSinCos_Test(t *testing.T, BasisType RTBasisType, PMax int) {
	var (
		dt VectorTestField
	)
	dt = SinCosVectorField{}

	t.Log("Begin divergence Test")
	PStart := 1
	PEnd := PMax
	for P := PStart; P <= PEnd; P++ {
		t.Logf("---------------------------------------------\n")
		t.Logf("Checking divergence for RT%d\n", P)
		t.Logf("---------------------------------------------\n")
		rt := NewRTElement(P, BasisType)
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
			calcDiv.Transpose().Print("Calculated divergence")
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

func RTInterpolation_Test(t *testing.T, BasisType RTBasisType, PMax int) {
	// Verify the interpolation of a constant vector field onto the element
	PStart := 1
	PEnd := PMax
	for P := PStart; P <= PEnd; P++ {
		t.Logf("---------------------------------------------\n")
		t.Logf("Checking Interpolation for RT%d\n", P)
		t.Logf("---------------------------------------------\n")
		var (
			dt VectorTestField
		)
		dt = PolyVectorField{}

		rt := NewRTElement(P, BasisType)
		if testing.Verbose() {
			rt.V.Print("V")
			rt.VInv.Print("VInv")
		}
		s1, s2 := make([]float64, rt.Np), make([]float64, rt.Np)
		for PField := 0; PField <= P; PField++ {
			t.Logf("\nReference Vector Field Order:%d\n", PField)
			t.Logf("-------------------------------\n")
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
				b_i := rt.DOFVectors[i].Eval(0, 0)
				// Sum of the basis polynomials over j, each dotted with basis vector_i
				for j := 0; j < rt.Np; j++ {
					f_rt_dot[i] += rt.Psi[j].dot(r_i, s_i, b_i) * C.At(j, 0)
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
