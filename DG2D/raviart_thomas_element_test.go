package DG2D

import (
	"fmt"
	"math"
	"strconv"
	"testing"
	"time"

	"github.com/notargets/gocfd/DG1D"
	"github.com/notargets/gocfd/utils"

	utils2 "github.com/notargets/avs/utils"

	"github.com/notargets/avs/chart2d"
	graphics2D "github.com/notargets/avs/geometry"

	"github.com/stretchr/testify/assert"
)

func TestRTElementErvinRT1(t *testing.T) {
	var (
		P      = 1
		Np     = (P + 1) * (P + 3)
		NpInt  = (P) * (P + 1) / 2
		NpEdge = P + 1
		g1     = 0.5 - math.Sqrt(3)/6
		g2     = 0.5 + math.Sqrt(3)/6
	)
	conv := func(r float64) (xi float64) {
		xi = (r + 1) / 2
		return
	}
	scalarMult := func(p float64, v [2]float64) (v2 [2]float64) {
		v2 = [2]float64{p * v[0], p * v[1]}
		return
	}
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
		// xi, eta = conv(r), conv(s)
		)
		if len(derivO) > 0 {
			div = 0
			return
		}
		// v[0] = xi
		// v[1] = eta - 1
		v[0] = 0
		v[1] = -1
		return
	}
	e2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64,
		div float64) {
		// Hypotenuse
		var (
			// xi, eta = conv(r), conv(s)
			sr2 = math.Sqrt2
		)
		if len(derivO) > 0 {
			div = 0
			return
		}
		// v[0] = sr2 * xi
		// v[1] = sr2 * eta
		v[0] = 0.5 * sr2
		v[1] = 0.5 * sr2
		return
	}
	e3 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64,
		div float64) {
		// Left edge
		var (
		// xi, eta = conv(r), conv(s)
		)
		if len(derivO) > 0 {
			div = 0
			return
		}
		// v[0] = xi - 1
		// v[1] = eta
		v[0] = -1
		v[1] = 0
		return
	}
	e4 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64,
		div float64) {
		var (
			xi, eta = conv(r), conv(s)
		)
		if len(derivO) != 0 {
			// v0 = (r+1)/2 * (s+1)/2
			// dv0/dr = (1/2) * (s+1)/2 = (s+1)/4
			// v1 = (s+1)/2 * ((s-1)/2) = (1/4) * (s^2 - 1)
			// dv1/ds = 2 * s / 4 = s/2
			// div = (s+1)/4 + s/2 = (s+1)/4 + 2s/4 = (3s+1)/4
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
			// v0 = (r+1)/2 * ((r-1)/2) = (1/4) * (r^2 - 1)
			// dv0/dr = 2 * r / 4 = r/2
			// v1 = (r+1)/2 * (s+1)/2
			// dv1/ds = 1/2 * (r+1)/2 = (r+1)/4
			// div = (r+1)/4 + r/2 = (r+1)/4 + 2r/4 = (3*r+1)/4
			div = (3.*r + 1.) / 4.
			return
		}
		v[0] = xi * (xi - 1)
		v[1] = xi * eta
		return
	}
	psiInt1 := func(r, s float64, derivO ...DerivativeDirection) (
		v [2]float64, div float64) {
		v, div = e4(r, s)
		return
	}
	psiInt2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e5(r, s)
		return
	}
	psiEdge1_1 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e1(r, s)
		// TODO: calculate div contribution from l1
		v = scalarMult(l1(r), v)
		return
	}
	psiEdge1_2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e1(r, s)
		// TODO: calculate div contribution from l2
		v = scalarMult(l2(r), v)
		return
	}
	psiEdge2_1 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e2(r, s)
		// TODO: calculate div contribution from l1
		v = scalarMult(l1(s), v)
		return
	}
	psiEdge2_2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e2(r, s)
		// TODO: calculate div contribution from l2
		v = scalarMult(l2(s), v)
		return
	}
	psiEdge3_1 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e3(r, s)
		// TODO: calculate div contribution from l2
		v = scalarMult(l2(s), v)
		return
	}
	psiEdge3_2 := func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64) {
		v, div = e3(r, s)
		// TODO: calculate div contribution from l1
		v = scalarMult(l1(s), v)
		return
	}
	dot := func(v1, v2 [2]float64) (val float64) {
		val = v1[0]*v2[0] + v1[1]*v2[1]
		return
	}

	type psi func(r, s float64, derivO ...DerivativeDirection) (v [2]float64, div float64)

	psi_j := []psi{psiInt1, psiInt2, psiEdge1_1, psiEdge1_2, psiEdge2_1,
		psiEdge2_2, psiEdge3_1, psiEdge3_2}

	RInt, SInt := NodesEpsilon(P)
	R, S := utils.NewVector(Np), utils.NewVector(Np)
	fmt.Printf("NpInt = %d\n", NpInt)
	for i := 0; i < 2*NpInt; i++ {
		R.DataP[i] = RInt.DataP[i]
		S.DataP[i] = SInt.DataP[i]
	}
	rconv := func(xi float64) (r float64) {
		r = 2*xi - 1
		return
	}
	i := 2 * NpInt
	// fmt.Printf("g1,g2 = %f, %f\n", g1, g2)
	// fmt.Printf("rconv(g1,g2) = %f, %f\n", rconv(g1), rconv(g2))
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
	// R.Transpose().Print("R")
	// S.Transpose().Print("S")

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
		default:
			fmt.Printf("j = %d\n", j)
			panic("wrong j")
		}
		return
	}
	V := utils.NewMatrix(Np, Np)
	for i = 0; i < Np; i++ {
		r_i, s_i := R.DataP[i], S.DataP[i]
		v_i, _ := psi_j[i](r_i, s_i)
		// fmt.Printf("v_%d[%f,%f] = [%f,%f]\n", i, r_i, s_i, v_i[0], v_i[1])
		for j := 0; j < Np; j++ {
			// Don't evaluate edge basis on other edges
			switch edgeNum(j) {
			case E1, E2, E3:
				switch edgeNum(i) {
				case E1, E2, E3:
					if edgeNum(i) != edgeNum(j) {
						continue
					}
				}
			}
			v_j, _ := psi_j[j](r_i, s_i)
			V.Set(i, j, dot(v_j, v_i))
		}
	}
	// i = 2 * NpInt
	// r_1, s_1 := R.DataP[i], S.DataP[i] // Edge 1 Node 1
	// j := 2 * NpInt                     // Edge 1 Psi_1_1
	// v0 := psi_j[j](r_1, s_1)
	// j = 2*NpInt + NpEdge // Edge 2 Psi_2_1
	// v1 := psi_j[j](r_1, s_1)
	// a_3_5 := dot(v0, v1)
	// fmt.Printf("v0 = [%f,%f], v1 = [%f,%f], a_3_5 = %f\n",
	// 	v0[0], v0[1], v1[0], v1[1], a_3_5)

	V.Print("V")
	V.InverseWithCheck().Print("Vinv")

	// V := utils.NewMatrix(Np, Np, []float64{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
	// 	16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
	// 	33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
	// 	51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63})
	// V.Print("V")
}

func TestRTElementPerformanceRT2(t *testing.T) {
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
	rt.BasisVector[0].Transpose().Print("Basis Vector 0")
	rt.BasisVector[1].Transpose().Print("Basis Vector 1")
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

func TestRTElementPerformanceRT1(t *testing.T) {
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
	rt.BasisVector[0].Transpose().Print("Basis Vector 0")
	rt.BasisVector[1].Transpose().Print("Basis Vector 1")
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
		// fmt.Printf("F[%f,%f]=[%f,%f], divF[%f,%f]=%f\n", r, s, f1, f2, r, s, dF)
	}
	dFcalc := utils.NewMatrix(Np, 1, divFcalc)
	dFcalc.Transpose().Print("Reference Div")
	rt.ProjectFunctionOntoDOF(s1, s2)
	dB := rt.Projection
	rt.Div.Mul(dB).Transpose().Print("Calculated Divergence")
}

func TestRTElementConstruction3(t *testing.T) {
	// Check the Lagrange Polynomial basis to verify the Lagrange property
	for P := 1; P < 7; P++ {
		R, S := NodesEpsilon(P)
		lp2d := NewLagrangePolynomialBasis2D(P, R, S)
		A := utils.NewMatrix(lp2d.Np, lp2d.Np)
		// Evaluate the j-th lagrange polynomial at (r,s)i
		// It should evaluate to 1 at each (r,s)i=j and 0 at (r,s)i!=j
		// In other words, this should be the diagonal matrix
		for j := 0; j < lp2d.Np; j++ {
			for i := 0; i < lp2d.Np; i++ {
				r, s := R.AtVec(i), S.AtVec(i)
				// fmt.Printf("psi[%d][%f,%f] = %f\n", j, r, s,
				// 	lp2d.GetPolynomialEvaluation(r, s, j))
				A.Set(i, j, lp2d.GetPolynomialEvaluation(r, s, j))
			}
		}
		checkIfUnitMatrix(t, A)
	}
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
		v, _ = rt.basisEvaluation(r, s, j)
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
		// v = [2]float64{s2 * xi, s2 * eta} // Ervin
		v = [2]float64{0.5 * s2, 0.5 * s2} // Unit normal
		return
	}
	e2 := func(xi, eta float64) (v [2]float64) {
		v = [2]float64{-1, 0} // Unit normal
		// v = [2]float64{xi - 1, eta - 0.5} // Ervin
		// v = [2]float64{xi - 1, eta}
		return
	}
	e3 := func(xi, eta float64) (v [2]float64) {
		v = [2]float64{0, -1} // Unit normal
		// v = [2]float64{xi - 0.5, eta - 1} // Ervin
		// v = [2]float64{xi, eta - 1}
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
		assert.InDeltaf(t, BasisVector[0].At(i, 0), BV0E.AtVec(i), 0.0001, "")
		assert.InDeltaf(t, BasisVector[1].At(i, 0), BV1E.AtVec(i), 0.0001, "")
	}
}

func TestRTElementLagrangePolynomials(t *testing.T) {
	// Check that the edge lagrange polynomials are zero at i!=j
	for P := 1; P < 7; P++ {
		rt := NewRTElement(P)
		offset := 2 * rt.NpInt
		for j := offset; j < 4*rt.NpEdge; j++ {
			for i := offset; i < 4*rt.NpEdge; i++ {
				r, s := rt.R.AtVec(i), rt.S.AtVec(i)
				var checkVal float64
				if i == j {
					checkVal = 1.
				} else {
					checkVal = 0.
				}
				val := rt.basisPolynomialValue(r, s, j)
				assert.InDeltaf(t, val, checkVal, 0.00001, "")
			}
		}
	}
}
func TestRTElementConstruction1(t *testing.T) {
	P := 3
	rt := NewRTElement(P)
	j1_1 := rt.BasisVector[0].At(0, 0)
	j1_2 := rt.BasisVector[1].At(0, 0)
	j2_1 := rt.BasisVector[0].At(rt.NpInt, 0)
	j2_2 := rt.BasisVector[1].At(rt.NpInt, 0)
	// Both basis 1 and 2 are at the same location in RT, the single interior
	// point, so we don't need to change the location
	dot1 := j1_1*j1_1 + j1_2*j1_2                 // Basis1 dotted with itself
	dot2 := j1_1*j2_1 + j1_2*j2_2                 // Basis1 dotted with Basis2
	assert.False(t, math.Abs(dot1-dot2) < 0.0001) // Should not be equal
	BasisMatrix := [2]utils.Matrix{utils.NewMatrix(rt.Np, rt.Np),
		utils.NewMatrix(rt.Np, rt.Np)}
	for i := 0; i < rt.Np; i++ {
		r, s := rt.R.AtVec(i), rt.S.AtVec(i)
		v := [2]float64{rt.BasisVector[0].At(i, 0), rt.BasisVector[1].At(i, 0)}
		// fmt.Printf("Base[%d] = [%f,%f]\n", i, v[0], v[1])
		for j := 0; j < rt.Np; j++ {
			v2, _ := rt.basisEvaluation(r, s, j)
			// fmt.Printf("Psi[%d,%d] = [%f,%f]\n", i, j, v2[0], v2[1])
			BasisMatrix[0].Set(i, j, v[0]*v2[0])
			BasisMatrix[1].Set(i, j, v[1]*v2[1])
		}
	}
	BasisDot := BasisMatrix[0].Add(BasisMatrix[1])
	assert.InDeltaf(t, BasisDot.At(0, 0), dot1, 0.0001, "")
	assert.InDeltaf(t, BasisDot.At(0, rt.NpInt), dot2, 0.0001, "")
	// BasisMatrix[0].Print("BasisMatrix 0")
	// BasisMatrix[1].Print("BasisMatrix 1")
	// BasisDot.Print("BasisDot")
	// BasisDot.InverseWithCheck().Print("BasisDot Inverse")
}

func TestRTElementConstruction(t *testing.T) {
	// Define an RT element at order P
	P := 2
	rt := NewRTElement(P)

	// Verify that edges start at the bottom edge (edge 1) and proceed in a
	// counterclockwise fashion as we increase the index
	Print := func(label string, iip *int) (edge [2][]float64) {
		edge[0] = make([]float64, rt.NpEdge)
		edge[1] = make([]float64, rt.NpEdge)
		fmt.Printf("%s", label)
		for i := 0; i < rt.NpEdge; i++ {
			fmt.Printf("[%f,%f] ", rt.R.DataP[*iip], rt.S.DataP[*iip])
			edge[0][i] = rt.R.DataP[*iip]
			edge[1][i] = rt.S.DataP[*iip]
			*iip++
		}
		fmt.Printf("\n")
		return
	}
	// Bottom edge (edge 1) - S should be -1, edge runs left to right
	ii := 2 * rt.NpInt
	edge1a := Print("Edge 1:", &ii)
	assert.True(t, nearVec(edge1a[0], []float64{-.774597, 0, 0.774597},
		0.00001))
	assert.True(t, nearVec(edge1a[1], []float64{-1, -1, -1}, 0.00001))
	// Hypotenuse edge (edge 2) - R should go from right to left, opposite for s
	edge2a := Print("Edge 2:", &ii)
	assert.True(t, nearVec(edge2a[0], []float64{.774597, 0, -0.774597},
		0.00001))
	assert.True(t, nearVec(edge2a[1], []float64{-.774597, 0, 0.774597},
		0.00001))
	// Left edge (edge 3) - R should be -1, edge runs right to left
	edge3a := Print("Edge 3:", &ii)
	assert.True(t, nearVec(edge3a[0], []float64{-1, -1, -1}, 0.00001))
	assert.True(t, nearVec(edge3a[1], []float64{.774597, 0, -0.774597},
		0.00001))

	// Note that Ervin's edges are numbered differently:
	// Ervin's E1 => E2 (Hypotenuse)
	// Ervin's E2 => E3 (Left)
	// Ervin's E3 => E1 (Bottom)

	// 1D Edge Polynomials
	// Get the edge values for edge1,2,3
	assert.Panics(t, func() { rt.getEdgeXiParameter(0, 0, 0) })

	// Test polynomial bases

	// TEST1: Demonstrate each of the ways to calculate Polynomial Terms
	// 2D Interior Polynomials
	PInt := rt.P - 1
	RInt, SInt := NodesEpsilon(PInt)
	BasisTest := func(basis *JacobiBasis2D) (PSI utils.Vector, P_Alt []float64) {
		PSI = basis.GetAllPolynomials()
		CalcTerm := func(r, s float64, P int) (psi float64) {
			for i := 0; i <= P; i++ {
				for j := 0; j <= (P - i); j++ {
					pTerm := basis.PolynomialTerm(r, s, i, j)
					psi += pTerm
				}
			}
			return
		}
		r1, s1 := RInt.AtVec(0), SInt.AtVec(0)
		r2, s2 := RInt.AtVec(1), SInt.AtVec(1)
		r3, s3 := RInt.AtVec(2), SInt.AtVec(2)
		P_Alt = []float64{
			CalcTerm(r1, s1, PInt),
			CalcTerm(r2, s2, PInt),
			CalcTerm(r3, s3, PInt),
		}
		for i := 0; i < len(P_Alt); i++ {
			fmt.Printf("P_Alt[%d] = %f\n", i, P_Alt[i])
		}
		return
	}
	jb2d := NewJacobiBasis2D(rt.P-1, rt.RInt, rt.SInt, 0, 0)
	PSI, P_Alt := BasisTest(jb2d)
	assert.True(t, nearVec(P_Alt, PSI.DataP, 0.00001))

	j := 0
	i := 0
	r, s := rt.RInt.AtVec(i), rt.SInt.AtVec(i)
	fmt.Printf("Poly(%d)[%f,%f] = %f\n", j, r, s,
		rt.basisPolynomialValue(r, s, j))
	// Build a polynomial matrix for interior polynomials
	A := utils.NewMatrix(2*rt.NpInt, 2*rt.NpInt)
	for i = 0; i < 2*rt.NpInt; i++ {
		r, s = rt.R.AtVec(i), rt.S.AtVec(i)
		for j = 0; j < 2*rt.NpInt; j++ {
			// fmt.Printf("NpInt, j = %d, %d\n", rt.NpInt, j)
			A.Set(i, j, rt.basisPolynomialValue(r, s, j))
		}
	}
	A.Print("Interior Poly")
	// Let's call out some important features, the polynomial for the first
	// 0 <= j < NpInt points should be the same as the polynomial for the
	// NpInt <= j < 2*NpInt points, specifically the Alpha and Beta params
	// should make that the case. The [r,s] coordinates for 0 <= i < NpInt
	// are the same as the [r,s] coordinates for NpInt <= i < 2*NpInt, so the
	// The polynomial matrix should look like this:
	// a b c a b c
	// d e f d e f
	// g h i g h i
	// a b c a b c
	// d e f d e f
	// g h i g h i
	// The top NpInt terms repeat in the bottom, and the left/right
	// The terms multiplying the e4 and e5 vectors should be the same.
	// In the Ervin paper, the 2D polynomial bj[r,s] is distinct over all values
	// of j from 0 to NpInt (equivalent to 1, k(k+1)/2) which this
	// construction achieves:
	// P = 2 = k, NpInt=k(k+1)/2=2(3)/2=3, which is NpInt here (not 2*NpInt)
	assert.True(t, near(A.At(0, 0), A.At(0, rt.NpInt), 0.00001))
	assert.True(t, near(A.At(0, 0), A.At(rt.NpInt, 0), 0.00001))
	// Check the RT2 polynomial (special case for RT2)
	xi := make([]float64, 3)
	eta := make([]float64, 3)
	for i = 0; i < 2; i++ {
		xi[i] = 0.5 * (rt.RInt.AtVec(i) + 1)
		eta[i] = 0.5 * (rt.SInt.AtVec(i) + 1)
	}
	for i = 0; i < 2; i++ {
		pp := []float64{1 - xi[i] - eta[i], xi[i], eta[i]}
		assert.True(t, nearVec(A.Row(i).DataP[0:2], pp, 0.00001))
	}

	// Build polynomial derivative matrices
	dr := utils.NewMatrix(rt.Np, rt.Np)
	ds := utils.NewMatrix(rt.Np, rt.Np)
	for i = 0; i < rt.Np; i++ {
		r, s = rt.R.AtVec(i), rt.S.AtVec(i)
		fmt.Printf("[%f,%f]%d\n", r, s, i)
	}
	for i = 0; i < rt.Np; i++ {
		r, s = rt.R.AtVec(i), rt.S.AtVec(i)
		for j = 0; j < rt.Np; j++ {
			dr.Set(i, j, rt.basisPolynomialValue(r, s, j, Dr))
			ds.Set(i, j, rt.basisPolynomialValue(r, s, j, Ds))
		}
	}
	fmt.Printf("RT%d\n", P)
	dr.Print("Poly Dr")
	ds.Print("Poly Ds")
}

func TestErvinBasisFunctions2(t *testing.T) {
	R := []float64{1. / 3., 0.5, 2. / 3.}
	// assert.Equal(t, 1., DG1D.Lagrange1DPoly(1./3., R, 0))
	// assert.Equal(t, 0., DG1D.Lagrange1DPoly(1./3., R, 1))
	// assert.Equal(t, 0., DG1D.Lagrange1DPoly(1./3., R, 2))
	// assert.Panics(t, func() { DG1D.Lagrange1DPoly(1./3., R, 3) })
	// assert.InDeltaf(t, -9., DG1D.Lagrange1DPoly(1./3., R, 0, 1), 0.000001, "")
	// assert.InDeltaf(t, 12., DG1D.Lagrange1DPoly(1./3., R, 1, 1), 0.000001, "")
	// assert.InDeltaf(t, -3., DG1D.Lagrange1DPoly(1./3., R, 2, 1), 0.000001, "")

	// Generate Gauss Lobato points for P=5 to compare with the online article:
	// https://math.stackexchange.com/questions/1105160/evaluate-derivative-of-lagrange-polynomials-at-construction-points
	R = DG1D.JacobiGL(0, 0, 6).DataP
	// One row (i) is the evaluation of the j-th derivative at each i-th point
	validation_deriv := make([][]float64, len(R))
	for i := range validation_deriv {
		validation_deriv[i] = make([]float64, len(R))
	}
	validation_deriv[0] = []float64{-10.5, -2.4429, 0.6253, -0.3125, 0.2261,
		-0.2266, 0.5}
	validation_deriv[1] = []float64{14.2016, 0, -2.2158, 0.9075, -0.6164,
		0.6022, -1.3174}
	validation_deriv[2] = []float64{-5.669, 3.4558, 0, -2.007, 1.0664, -0.9613,
		2.05}
	validation_deriv[3] = []float64{3.2, -1.5986, 2.2667, 0, -2.2667, 1.5986,
		-3.2}
	validation_deriv[4] = []float64{-2.05, 0.9613, -1.0664, 2.007, 0, -3.4558,
		5.669}
	validation_deriv[5] = []float64{1.3174, -0.6022, 0.6164, -0.9075, 2.2158, 0,
		-14.2016}
	validation_deriv[6] = []float64{-0.5, 0.2266, -0.2261, 0.3125, -0.6253,
		2.4429, 10.5}
	// testVec := make([]float64, len(R))
	// for j, validate := range validation_deriv {
	// 	for i, r := range R {
	// 		testVec[i] = DG1D.Lagrange1DPoly(r, R, j, 1)
	// 	}
	// 	assert.True(t, nearVec(testVec, validate, 0.0001))
	// }
}
func TestErvinBasisFunctions1(t *testing.T) {
	// This tests the basic basis functions e1,e2,e3 for edges and e4,
	// e5 interior
	var (
		sr2   = math.Sqrt(2.)
		oosr2 = 1. / sr2
	)
	conv := func(r, s float64) (xiEta [2]float64) {
		xiEta = [2]float64{0.5 * (r + 1), 0.5 * (s + 1)}
		return
	}
	assert.Equal(t, conv(-1, -1), [2]float64{0., 0.})
	assert.Equal(t, conv(-1, 1), [2]float64{0., 1.})
	assert.Equal(t, conv(1, -1), [2]float64{1., 0.})
	rt := NewRTElement(2)

	// ---------------------------------------------------
	// the edge basis vector function [v] varies along the edge.
	// It is the product of a 1D edge function f(xi) and [v], so we have:
	// div(edgeFunction) = div(f(xi)*[v]) =
	//       [df(xi)/dr,df(xi)/ds] ⋅ [v] + f(xi) * ([div] ⋅ [v])
	//
	// div = df(xi)/dxi * (v1*dxi/dr + v2*dxi/ds) + f(xi) * (dv1/dr + dv2/ds)
	//
	// Conversion from Ervin coordinates:
	// xi  = 0.5 * (r + 1)
	// eta = 0.5 * (s + 1)
	// r = 2 * xi  - 1
	// s = 2 * eta - 1
	//
	// Left Edge (1) divergence:
	// 			[v] = [(r - 1)/2, s/2]
	// v1 = (r-1)/2, v2 = s/2, dv1/dr = 1/2, dv2/ds = 1/2
	// The left edge  in a counter-clockwise direction is parameterized:
	// r = -1 (constant), xi = -s => dxi/dr = 0, dxi/ds = -1,
	// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
	//     = df/dxi*(v1*(  0   )+v2*(  -1  )) + f(xi) * (  1/2  +   1/2 )
	//     = df/dxi*(          v2           ) + f(xi)
	//         div(edge1) = (df/dxi) * v2 + f(xi)
	//
	// Hypotenuse (2) divergence:
	// 			[v] = [Sqrt2/2 * (r+1), Sqrt2/2 * (s+1)]
	// v1 = Sqrt2/2 * (r+1), v2 = Sqrt2/2 * (s+1), dv1/dr = Sqrt2/2 = dv2/ds
	//
	// The hypotenuse in a counter-clockwise direction is parameterized:
	// xi = -r = s, => dxi/dr = -1, dxi/ds = 1
	//
	// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
	//     = df/dxi*(v1*( -1   )+v2*(   1  )) + f(xi) * (Sqrt2/2+Sqrt2/2)
	//     = df/dxi*(         v2-v1         ) + f(xi) * Sqrt2
	//         div(edge2) = (df/dxi) * (v2-v1) + Sqrt2 * f(xi)
	//
	// Bottom Edge (3) divergence:
	// 			[v] = [r/2, (s - 1)/2]
	// v1 = r/2, v2 = (s-1)/2, dv1/dr = 1/2, dv2/ds = 1/2
	// The bottom edge  in a counter-clockwise direction is parameterized:
	// xi = r, s = -1 (const) => dxi/dr = 1, dxi/ds = 0
	// div = df/dxi*(v1*(dxi/dr)+v2*(dxi/ds)) + f(xi) * (dv1/dr + dv2/ds)
	//     = df/dxi*(v1*(  1   )+v2*(   0  )) + f(xi) * (  1/2  +   1/2 )
	//     = df/dxi*(          v1           ) + f(xi)
	//        div(edge3) = (df/dxi) * v1 + f(xi)

	// Test the dot product of the e4 and e5 interior basis vectors against all
	// edge normals at the edge normal locations. All interior basis functions
	// should have zero dot product on all edges
	dot := func(v1, v2 [2]float64) (dp float64) {
		dp = v1[0]*v2[0] + v1[1]*v2[1]
		return
	}

	getMid := func(start, end [2]float64) (mid [2]float64) {
		mid[0] = (start[0] + end[0]) / 2
		mid[1] = (start[1] + end[1]) / 2
		return
	}

	offsetEdge1 := 2 * rt.NpInt             // Bottom edge
	offsetEdge2 := 2*rt.NpInt + rt.NpEdge   // Hypotenuse
	offsetEdge3 := 2*rt.NpInt + 2*rt.NpEdge // Left edge

	// Edge 1 endpoint 1 (start of path around triangle)
	e1Nd1 := [2]float64{-1, -1}
	e1Nd2 := [2]float64{1, -1}
	e2Nd1 := e1Nd2
	e2Nd2 := [2]float64{-1, 1}
	e3Nd1 := e2Nd2
	e3Nd2 := e1Nd1
	e1Mid := getMid(e1Nd1, e1Nd2)
	e2Mid := getMid(e2Nd1, e2Nd2)
	e3Mid := getMid(e3Nd1, e3Nd2)
	ef1 := rt.baseBasisVectors(e1Mid[0], e1Mid[1], offsetEdge1)
	ef2 := rt.baseBasisVectors(e2Mid[0], e2Mid[1], offsetEdge2)
	ef3 := rt.baseBasisVectors(e3Mid[0], e3Mid[1], offsetEdge3)

	ef4 := rt.baseBasisVectors(e1Nd1[0], e1Nd1[1], 0)
	ef5 := rt.baseBasisVectors(e1Nd1[0], e1Nd1[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef1))
	assert.Equal(t, 0., dot(ef5, ef1))
	// Edge 1 midpoint
	ef4 = rt.baseBasisVectors(e1Mid[0], e1Mid[1], 0)
	ef5 = rt.baseBasisVectors(e1Mid[0], e1Mid[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef1))
	// Edge 1 endpoint 2
	ef4 = rt.baseBasisVectors(e1Nd2[0], e1Nd2[1], 0)
	ef5 = rt.baseBasisVectors(e1Nd2[0], e1Nd2[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef1))
	assert.Equal(t, 0., dot(ef5, ef1))

	// Edge 2 midpoint
	ef4 = rt.baseBasisVectors(e2Mid[0], e2Mid[1], 0)
	ef5 = rt.baseBasisVectors(e2Mid[0], e2Mid[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef2))
	assert.Equal(t, 0., dot(ef5, ef2))
	// Edge 2 endpoint 1
	ef4 = rt.baseBasisVectors(e2Nd1[0], e2Nd1[1], 0)
	ef5 = rt.baseBasisVectors(e2Nd1[0], e2Nd1[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef2))
	assert.Equal(t, 0., dot(ef5, ef2))
	// Edge 2 endpoint 2
	ef4 = rt.baseBasisVectors(e2Nd2[0], e2Nd2[1], 0)
	ef5 = rt.baseBasisVectors(e2Nd2[0], e2Nd2[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef2))
	assert.Equal(t, 0., dot(ef5, ef2))

	// Edge 3 midpoint
	ef4 = rt.baseBasisVectors(e3Mid[0], e3Mid[1], 0)
	ef5 = rt.baseBasisVectors(e3Mid[0], e3Mid[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef3))
	assert.Equal(t, 0., dot(ef5, ef3))
	// Edge 3 endpoint 1
	ef4 = rt.baseBasisVectors(e3Nd1[0], e3Nd1[1], 0)
	ef5 = rt.baseBasisVectors(e3Nd1[0], e3Nd1[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef3))
	assert.Equal(t, 0., dot(ef5, ef3))
	// Edge 3 endpoint 2
	ef4 = rt.baseBasisVectors(e3Nd2[0], e3Nd2[1], 0)
	ef5 = rt.baseBasisVectors(e3Nd2[0], e3Nd2[1], rt.NpInt)
	assert.Equal(t, 0., dot(ef4, ef3))
	assert.Equal(t, 0., dot(ef5, ef3))

	// Test the midpoint of each edge for values of the base edge functions
	// Edge midpoints, they should be the unit normals
	assert.Equal(t, [2]float64{0, -1}, ef1)
	assert.InDeltaf(t, oosr2, ef2[0], 0.0000001, "")
	assert.InDeltaf(t, oosr2, ef2[1], 0.0000001, "")
	assert.Equal(t, [2]float64{-1, 0}, ef3)
}

func TestLagrangePolynomial(t *testing.T) {
	numSamples := 100
	rd := make([]float64, numSamples)
	xmin, xmax := -1., 1.
	fmin, fmax := -0.5, 1.25
	inc := (xmax - xmin) / float64(numSamples-1.)
	for i := 0; i < numSamples; i++ {
		rd[i] = xmin + float64(i)*inc
	}
	SamplesR := utils.NewVector(numSamples, rd)
	f := make([]float64, SamplesR.Len())
	var plot bool
	plot = false
	if plot {
		chart := utils.NewLineChart(1920, 1080, xmin, xmax, fmin, fmax)
		// TODO: Make a pluggable basis underneath the RT (and Lagrange) elements - Lagrange, Hesthaven, Spectral?
		var delay time.Duration
		Nmax := 4
		lineInc := 2. / float64(Nmax-2)
		lineColor := -1. // colormap goes from -1,1
		for n := 1; n < Nmax; n++ {
			Np := n + 1
			R := utils.NewVector(Np)
			inc = (xmax - xmin) / float64(Np-1)
			for i := 0; i < Np; i++ {
				R.DataP[i] = xmin + float64(i)*inc
			}
			lp := NewLagrangeBasis1D(R.DataP)
			_ = lp
			for j := 0; j < R.Len(); j++ {
				for i, r := range SamplesR.DataP {
					// f[i] = DG1D.JacobiP(utils.NewVector(1, []float64{r}), 0, 0, j)[0]
					f[i] = lp.BasisPolynomial([]float64{r}, j)[0]
				}
				if n == Nmax-1 && j == R.Len()-1 {
					delay = 120 * time.Second
				}
				name := "JacobiP[" + strconv.Itoa(n) + "," + strconv.Itoa(j) + "]"
				fmt.Printf("Chart Name: [%s], lineColor = %5.3f\n", name, lineColor)
				chart.Plot(delay, SamplesR.DataP, f, lineColor, name)
			}
			lineColor += lineInc
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
	if true { // Check Divergence for polynomial vector fields of order < N
		// against analytical solution
		// Nend := 8
		// for N := 1; N < Nend; N++ {
		N := 1
		rt := NewRTElement(N)
		for cOrder := 0; cOrder <= N; cOrder++ {
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
