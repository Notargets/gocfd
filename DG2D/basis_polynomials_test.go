package DG2D

import (
	"testing"

	"github.com/stretchr/testify/assert"

	"github.com/notargets/gocfd/utils"
)

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
