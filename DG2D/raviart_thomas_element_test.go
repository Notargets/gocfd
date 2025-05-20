package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

func TestRTElementConvergence(t *testing.T) {
	DivergenceConvergence_Test(t, SimplexRTBasis, 7)
}

func TestRTElement(t *testing.T) {
	for _, rtb := range []RTBasisType{SimplexRTBasis} {
		// for _, rtb := range []RTBasisType{ErvinBasisRT, SimplexRTBasis} {
		// for _, rtb := range []RTBasisType{ErvinBasisRT} {
		var PMin, PMax int
		switch rtb {
		case ErvinBasisRT:
			PMin = 1
			PMax = 2
		case SimplexRTBasis:
			PMin = 1
			PMax = 7
		}
		t.Logf("===================> %s\n", rtb.String())
		DivergencePolynomialField_Test(t, rtb, PMin, PMax)
	}
}

func DivergencePolynomialField_Test(t *testing.T, BasisType RTBasisType, PMin, PMax int) {
	var (
		dt []VectorTestField
	)
	dt = []VectorTestField{
		PolyVectorField{},
		PolyVectorField2{},
		PolyVectorField3{},
	}
	for _, field := range dt {
		t.Logf("Begin Divergence Test for [%s]", field.String())
		// P := 1
		PStart := PMin
		PEnd := PMax
		for P := PStart; P <= PEnd; P++ {
			t.Logf("---------------------------------------------\n")
			t.Logf("Checking Divergence for RT%d\n", P)
			t.Logf("---------------------------------------------\n")
			rt := NewRTElement(P, BasisType, WSJ)
			CheckDivergence(t, rt, field, 0, P)
		}
	}
}

func DivergenceConvergence_Test(t *testing.T, BasisType RTBasisType, PMax int) {
	PStart := 1
	PEnd := PMax
	rmsCheck := make([]float64, PEnd-PStart+1)
	for P := PStart; P <= PEnd; P++ {
		t.Logf("---------------------------------------------\n")
		t.Logf("Checking RMS convergence for RT%d\n", P)
		t.Logf("---------------------------------------------\n")
		rt := NewRTElement(P, BasisType, WSJ)
		checkField := SinCosVectorField{}
		rmsCheck[P-1] = CheckDivergenceRMS(t, rt, checkField)
		fmt.Printf("RMS Error on %s: %f\n", checkField.String(), rmsCheck[P-1])
	}
	for P := PStart; P <= PEnd; P++ {
		fmt.Printf("Log10 RMS Error(%d) =%f\n", P, math.Log10(rmsCheck[P-1]))
		// Shows quadratic convergence in Log10 RMS error scaling in P
	}
	// Let'S do a least squares fit of the Log10 of the error to discover the
	// polynomial convergence order
	// l10RMS := make([]float64, len(rmsCheck))
	// for i := 0; i < len(rmsCheck); i++ {
	// 	l10RMS[i] = math.Log10(rmsCheck[i])
	// }
	// Results := utils.NewVector(PEnd, l10RMS)
	Results := utils.NewVector(PEnd, rmsCheck)
	jp1d := NewJacobiBasis1D(PEnd-1, Results, 0, 0)
	V := jp1d.Vandermonde1D()
	V.Print("V")
	VInter, _ := V.Transpose().Mul(V).Inverse()
	ls := VInter.Mul(V.Transpose()).Mul(Results.ToMatrix())
	ls.Transpose().Print("Polynomial Coefficients")
	// Quadratic dominant convergence is expected
	assert.True(t, ls.DataP[1] > 0.8)
	for i := 0; i < PEnd; i++ {
		if i != 1 {
			assert.True(t, ls.DataP[i] < 0.01)
		}
	}
	assert.True(t, ls.DataP[1] > 0.8)
}

func CheckDivergenceRMS(t *testing.T, rt *RTElement, dt VectorTestField) (rmsError float64) {
	var (
		Np = rt.Np
	)
	// A := utils.NewMatrix(Np, Np)
	f1, f2 := make([]float64, Np), make([]float64, Np)
	DivRef := make([]float64, Np)
	FProj := utils.NewMatrix(Np, 1)
	t.Logf("\nReference Vector Field Infinite Order: %s\n", dt.String())
	t.Logf("-------------------------------\n")
	for i := 0; i < Np; i++ {
		r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		f1[i], f2[i] = dt.F(r, s, 0)
		DivRef[i] = dt.Divergence(r, s, 0)
	}
	rt.ProjectFunctionOntoDOF(f1, f2, FProj.DataP)
	DivCalc := rt.Div.Mul(FProj)
	for i := 0; i < Np; i++ {
		baseErr := DivRef[i] - DivCalc.DataP[i]
		rmsError += baseErr * baseErr
	}
	rmsError = math.Sqrt(rmsError / float64(Np))
	return
}

func CheckDivergence(t *testing.T, rt *RTElement, dt VectorTestField,
	PFieldStart, PFieldEnd int) {
	var (
		Np   = rt.Np
		tolM = 0.00001
	)
	// A := utils.NewMatrix(Np, Np)
	f1, f2 := make([]float64, Np), make([]float64, Np)
	DivRef := make([]float64, Np)
	FProj := utils.NewMatrix(Np, 1)
	for PField := PFieldStart; PField <= PFieldEnd; PField++ {
		t.Logf("\nReference Vector Field Order:%d\n", PField)
		t.Logf("-------------------------------\n")
		for i := 0; i < Np; i++ {
			r, s := rt.R.AtVec(i), rt.S.AtVec(i)
			f1[i], f2[i] = dt.F(r, s, PField)
			DivRef[i] = dt.Divergence(r, s, PField)
		}
		rt.ProjectFunctionOntoDOF(f1, f2, FProj.DataP)
		DivCalc := rt.Div.Mul(FProj)
		if testing.Verbose() {
			DivCalc.Transpose().Print("Div calc")
			fmt.Printf("Ref Div = \n[ ")
			for i := 0; i < Np; i++ {
				fmt.Printf("%6.5f ", DivRef[i])
			}
			fmt.Printf("]\n")
		}
		// Calculate appropriate tolerance for check
		var maxF float64
		for i := 0; i < Np; i++ {
			maxF = math.Max(math.Abs(DivRef[i]), maxF)
		}
		tol := tolM * maxF
		if math.Abs(maxF) < tolM {
			tol = tolM
		}
		assert.InDeltaSlicef(t, DivRef, DivCalc.DataP, tol, "")
	}
}
