package isentropic_vortex

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestIVortex(t *testing.T) {
	{ // Test base state
		iv := NewIVortex(5, 5, 0, 1.4)
		rho, rhoU, rhoV, E := iv.GetStateC(0, 5, 0)
		// Compare to matlab script results (from isentropic_vortex.m)
		assert.True(t, nearVec([]float64{rho, rhoU, rhoV, E}, []float64{0.361673, 0.361673, 0.000000, 0.782817}, 0.00001))
		//fmt.Printf("Q = %v\n", []float64{rho, rhoU, rhoV, E})
		// Compare to matlab script results (from isentropic_vortex.m)
		Fx, Fy := iv.GetFlux(0, 5, 0)
		assert.True(t, near4Vec(Fx, [4]float64{0.361673, 0.602465, 0.000000, 1.023609}, 0.00001))
		assert.True(t, near4Vec(Fy, [4]float64{0.000000, 0.000000, 0.240792, 0.000000}, 0.00001))
		Div := iv.GetDivergence(0, 5, 0)
		assert.True(t, near4Vec(Div, [4]float64{0.000000, 0.000000, 0.782349, 0.000000}, 0.00001))
	}
	{ // Test values
		iv := NewIVortex(5, 5, 0, 1.4)
		rho, rhoU, rhoV, E := iv.GetStateC(3.99793, 10.00000, -0.14684)
		fmt.Printf("rho, rhoU, rhoV, E = %8.5f,%8.5f,%8.5f,%8.5f\n", rho, rhoU, rhoV, E)
		rho, rhoU, rhoV, E = iv.GetStateC(0, 5.00000, 0.0)
		fmt.Printf("rho, rhoU, rhoV, E = %8.5f,%8.5f,%8.5f,%8.5f\n", rho, rhoU, rhoV, E)
		rho, u, v, p := iv.GetState(4., 10.00000, 0.0)
		fmt.Printf("rho, u, v, p = %8.5f,%8.5f,%8.5f,%8.5f\n", rho, u, v, p)
	}
}

func near4Vec(a, b [4]float64, tol float64) (l bool) {
	A, B := make([]float64, 4), make([]float64, 4)
	for i := 0; i < 4; i++ {
		A[i] = a[i]
		B[i] = b[i]
	}
	return nearVec(A, B, tol)
}

func nearVec(a, b []float64, tol float64) (l bool) {
	for i, val := range a {
		if !near(b[i], val, tol) {
			fmt.Printf("Diff = %v, Left[%d] = %v, Right[%d] = %v\n", math.Abs(val-b[i]), i, val, i, b[i])
			return false
		}
	}
	return true
}

func nearVecScalar(a []float64, b float64, tol float64) (l bool) {
	for i, val := range a {
		if !near(b, val, tol) {
			fmt.Printf("Diff = %v, Left[%d] = %v, Right[%d] = %v\n", math.Abs(val-b), i, val, i, b)
			return false
		}
	}
	return true
}

func near(a, b float64, tolI ...float64) (l bool) {
	var (
		tol float64
	)
	if len(tolI) == 0 {
		tol = 1.e-08
	} else {
		tol = tolI[0]
	}
	bound := math.Max(tol, tol*math.Abs(a))
	if math.Abs(a-b) <= bound {
		l = true
	}
	return
}
