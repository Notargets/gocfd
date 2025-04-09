package DG2D

import (
	"math"
	"testing"

	"github.com/notargets/gocfd/DG1D"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

func TestInterpolation1D(t *testing.T) {
	// We have density on 6 interior points:
	// Interior point index:       0       1       2       3       4       5
	// Density values:        1.1691  1.0000  1.1691  1.0000  1.1691  1.0000
	// Interior point R:     -0.8168  0.6337 -0.8168 -0.1081 -0.7838 -0.1081
	// Interior point S:     -0.8168 -0.8168  0.6337 -0.7838 -0.1081 -0.1081
	// In a test case, we have the projected Edge2 (Left) projected distribution
	// of points on the edge from the RT3 interior as:
	// Interior point index:       0       1       3       4
	// Edge mapped locations:-0.8168  0.6337 -0.1081 -0.7838	      Edge 0 R
	// Interior point index:       1       2       3       4       5
	// Edge mapped locations: 0.7253 -0.7253  0.3378 -0.3378  0.0000  Edge 1 R
	// Edge mapped locations:-0.7253  0.7253 -0.3378  0.3378  0.0000  Edge 1 S
	// Interior point index:       0       2       3       4
	// Edge mapped locations:-0.8168  0.6337 -0.7838 -0.1081	      Edge 2 S
}

func TestJacobiBasis2D_IndividualTerms(t *testing.T) {
	tol := 0.000001
	P := 2
	R, S := NodesEpsilon(P)
	jb2d := NewJacobiBasis2D(P, R, S, 0, 0)
	if testing.Verbose() {
		jb2d.V.Print("JacobiBasis2D V")
	}

	// IJ := make([][2]int, jb2d.Np)
	// var sk int
	// for j := 0; j <= jb2d.P; j++ {
	// 	for i := 0; i <= (jb2d.P - j); i++ {
	// 		IJ[sk] = [2]int{j, i}
	// 		sk++
	// 	}
	// }
	// assert.Equal(t, sk, jb2d.Np)

	A := utils.NewMatrix(jb2d.Np, jb2d.Np)
	B := utils.NewMatrix(jb2d.Np, jb2d.Np)
	for j := 0; j < jb2d.Np; j++ {
		ij := jb2d.IJ[j]
		for i := 0; i < jb2d.Np; i++ {
			r, s := R.AtVec(i), S.AtVec(i)
			// fmt.Println("IJ = ", j, ij)
			A.Set(i, j, jb2d.PolynomialTerm(r, s, ij[0], ij[1]))
			B.Set(i, j, jb2d.GetPolynomialAtJ(r, s, j))
		}
	}
	if testing.Verbose() {
		A.Print("A")
		B.Print("B")
	}
	assert.InDeltaSlicef(t, jb2d.V.DataP, A.DataP, tol, "")
	assert.InDeltaSlicef(t, jb2d.V.DataP, B.DataP, tol, "")
}

func TestJacobiBasis2D_Gradient(t *testing.T) {
	var (
		scalarPoly = PolyScalarField{}
		tol        = 0.000001
	)
	PStart := 1
	PEnd := 6
	for P := PStart; P <= PEnd; P++ {
		t.Logf("----------------------------------------\n")
		t.Logf("Testing Order %d\n", P)
		t.Logf("----------------------------------------\n")
		R, S := NodesEpsilon(P)
		jb2d := NewJacobiBasis2D(P, R, S, 0, 0)
		PTestStart := 0
		PTestEnd := P
		// PTestEnd = 1
		for testP := PTestStart; testP <= PTestEnd; testP++ {
			t.Logf("Scalar Test Field Order %d\n", testP)
			t.Logf("----------------------------------------\n")
			Np := jb2d.Np
			field := make([]float64, Np)
			grad := make([][2]float64, Np)
			grad0 := make([]float64, Np)
			grad1 := make([]float64, Np)
			for i := 0; i < Np; i++ {
				r, s := R.AtVec(i), S.AtVec(i)
				field[i] = scalarPoly.P(r, s, testP)
				grad[i] = scalarPoly.Gradient(r, s, testP)
				grad0[i] = grad[i][0]
				grad1[i] = grad[i][1]
				// t.Logf("P[%f,%f]=%f\n", R, S, field[i])
			}
			VField := utils.NewMatrix(Np, 1, field)
			Grad0 := utils.NewVector(Np, grad0)
			Grad1 := utils.NewVector(Np, grad1)
			calcGrad0 := jb2d.Vr.Mul(jb2d.Vinv).Mul(VField)
			calcGrad1 := jb2d.Vs.Mul(jb2d.Vinv).Mul(VField)
			if testing.Verbose() {
				VField.Transpose().Print("VField")
				calcGrad0.Transpose().Print("VField Dr")
				Grad0.Transpose().Print("VTestField Grad0")
				calcGrad1.Transpose().Print("VField Ds")
				Grad1.Transpose().Print("VTestField Grad1")
			}
			assert.InDeltaSlicef(t, Grad0.DataP, calcGrad0.DataP, tol, "")
			assert.InDeltaSlicef(t, Grad1.DataP, calcGrad1.DataP, tol, "")
		}
	}

}

func TestJacobiBasis2D_GetOrthogonalPolynomialAtJ(t *testing.T) {
	tol := 0.000001
	PStart := 1
	PEnd := 6
	for P := PStart; P <= PEnd; P++ {
		R, S := NodesEpsilon(P)
		jb2d := NewJacobiBasis2D(P, R, S, 0, 0)
		A := utils.NewMatrix(jb2d.Np, jb2d.Np)
		for j := 0; j < jb2d.Np; j++ {
			for i := 0; i < jb2d.Np; i++ {
				r, s := R.AtVec(i), S.AtVec(i)
				A.Set(i, j, jb2d.GetOrthogonalPolynomialAtJ(r, s, j))
			}
		}
		if testing.Verbose() {
			A.Print("A")
		}
		assert.True(t, isIdentityMatrix(A, tol))
	}
}

func TestJacobiBasis1D_GetOrthogonalPolynomialAtJ(t *testing.T) {
	tol := 0.000001
	PStart := 0
	PEnd := 6
	for P := PStart; P <= PEnd; P++ {
		R := utils.NewVector(P+1, DG1D.LegendreZeros(P))
		jb1d := NewJacobiBasis1D(P, R, 0, 0)
		A := utils.NewMatrix(jb1d.Np, jb1d.Np)
		DR := utils.NewMatrix(jb1d.Np, jb1d.Np)
		for j := 0; j < jb1d.Np; j++ {
			for i := 0; i < jb1d.Np; i++ {
				r := R.AtVec(i)
				A.Set(i, j, jb1d.GetOrthogonalPolynomialAtJ(r, j))
				DR.Set(i, j, jb1d.GetOrthogonalPolynomialAtJ(r, j, Dr))
			}
		}
		if testing.Verbose() {
			A.Print("A")
			DR.Print("DR")
		}
		assert.True(t, isIdentityMatrix(A, tol))
		F_Ref := Poly1DForTests(R)
		DR_Ref := Poly1DForTests(R, Dr)
		DR_Calc := DR.Mul(F_Ref)
		if testing.Verbose() {
			F_Ref.Transpose().Print("F_Ref")
			DR_Ref.Transpose().Print("DR_Ref")
			DR_Calc.Transpose().Print("DR_Calc")
		}
		assert.InDeltaSlicef(t, DR_Ref.DataP, DR_Calc.DataP, tol, "")
	}
}

func TestLagrangePoly1D(t *testing.T) {
	tol := 0.000001
	PStart := 0
	PEnd := 6
	for P := PStart; P <= PEnd; P++ {
		Np := P + 1
		R := utils.NewVector(P+1, DG1D.LegendreZeros(P))
		A := utils.NewMatrix(Np, Np)
		DR := utils.NewMatrix(Np, Np)
		for j := 0; j < Np; j++ {
			for i := 0; i < Np; i++ {
				r := R.AtVec(i)
				A.Set(i, j, Lagrange1DPoly(r, R.DataP, j))
				DR.Set(i, j, Lagrange1DPoly(r, R.DataP, j, Dr))
			}
		}
		if testing.Verbose() {
			A.Print("A")
			DR.Print("DR")
		}
		assert.True(t, isIdentityMatrix(A, tol))
		F_Ref := Poly1DForTests(R)
		DR_Ref := Poly1DForTests(R, Dr)
		DR_Calc := DR.Mul(F_Ref)
		if testing.Verbose() {
			F_Ref.Transpose().Print("F_Ref")
			DR_Ref.Transpose().Print("DR_Ref")
			DR_Calc.Transpose().Print("DR_Calc")
		}
		assert.InDeltaSlicef(t, DR_Ref.DataP, DR_Calc.DataP, tol, "")
	}
}

func Poly1DForTests(R utils.Vector, derivO ...DerivativeDirection) (
	valM utils.Matrix) {
	var (
		Np = R.Len()
		P  = Np - 1
	)
	val := make([]float64, Np)
	for i := 0; i < Np; i++ {
		if len(derivO) == 0 {
			val[i] = math.Pow(R.AtVec(i), float64(P))
		} else {
			if P == 0 {
				val[i] = 0.
			} else {
				val[i] = float64(P) * math.Pow(R.AtVec(i), float64(P-1))
			}
		}
	}
	valM = utils.NewMatrix(Np, 1, val)
	return
}

func isIdentityMatrix(A utils.Matrix, epsilon float64) bool {
	n, _ := A.Dims()
	sum := 0.0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			diff := A.At(i, j)
			if i == j {
				diff -= 1 // A[i][i] - 1
			}
			sum += diff * diff
		}
	}
	return math.Sqrt(sum) < epsilon
}
