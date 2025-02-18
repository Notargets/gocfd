package DG2D

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

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
			PMax = 6
		}
		t.Logf("===================> %s\n", rtb.String())
		DivergencePolynomialField_Test(t, rtb, PMin, PMax)

		// t.Logf("Testing RT divergence on SinCos Fields for %v\n",
		// 	rtb.String())
		// RTDivergenceSinCos_Test(t, rtb, PMax)
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
			PFieldStart := 0
			PFieldEnd := P
			t.Logf("---------------------------------------------\n")
			t.Logf("Checking Divergence for RT%d\n", P)
			t.Logf("---------------------------------------------\n")
			rt := NewRTElement(P, BasisType)
			CheckDivergence(t, rt, field, PFieldStart, PFieldEnd)
		}
	}
}

func CheckDivergence(t *testing.T, rt *RTElement, dt VectorTestField,
	PFieldStart, PFieldEnd int) {
	var (
		Np  = rt.Np
		tol = 0.000001
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
		rt.ProjectFunctionOntoDOF(f1, f2, FProj)
		DivCalc := rt.Div.Mul(FProj)
		if testing.Verbose() {
			DivCalc.Transpose().Print("Div Calc")
			fmt.Printf("Ref Div = \n[ ")
			for i := 0; i < Np; i++ {
				fmt.Printf("%6.5f ", DivRef[i])
			}
			fmt.Printf("]\n")
		}
		assert.InDeltaSlicef(t, DivRef, DivCalc.DataP, tol, "")
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
		f1, f2 := make([]float64, Np), make([]float64, Np)
		// for PField := 0; PField <= (P - 1); PField++ {
		t.Logf("\nReference Vector Field Sin/Cos\n")
		t.Logf("-------------------------------\n")
		for i := 0; i < Np; i++ {
			r, s := rt.R.AtVec(i), rt.S.AtVec(i)
			f1[i], f2[i] = dt.F(r, s, 0)
			dF := dt.Divergence(r, s, 0)
			divFcalc[i] = dF
		}
		dFReference := utils.NewMatrix(Np, 1, divFcalc)
		if testing.Verbose() {
			dFReference.Transpose().Print("Reference Div")
		}
		FProj := utils.NewMatrix(Np, 1)
		rt.ProjectFunctionOntoDOF(f1, f2, FProj)
		calcDiv := rt.Div.Mul(FProj)
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
