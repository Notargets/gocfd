package DG2D

import (
	"math"
)

type ScalarTestField interface {
	P(r, s float64, P int) (val float64)
	Gradient(r, s float64, P int) (grad [2]float64)
}

type VectorTestField interface {
	F(r, s float64, P int) (f1, f2 float64)
	Divergence(r, s float64, P int) (div float64)
	String() string
}

type SinCosVectorField struct{}

func (scf SinCosVectorField) String() string { return "Sin / Cos Field" }

func (scf SinCosVectorField) F(r, s float64, P int) (f1, f2 float64) {
	var (
		Pi = math.Pi
	)
	f1, f2 = math.Sin(Pi*(r+1)), math.Cos(Pi*(s+1))
	return
}

func (scf SinCosVectorField) Divergence(r, s float64, P int) (div float64) {
	var (
		Pi = math.Pi
	)
	div = Pi * (math.Cos(Pi*(r+1)) - math.Sin(Pi*(s+1)))
	return
}

type PolyVectorField struct{}

func (lpf PolyVectorField) String() string { return "R + S Permutation Field" }

func (lpf PolyVectorField) F(r, s float64, P int) (f1, f2 float64) {
	var (
		p = float64(P)
	)
	f1, f2 = math.Pow(r+s+10, p), math.Pow(10*(r+s), p)
	return
}

func (lpf PolyVectorField) Divergence(r, s float64, P int) (div float64) {
	var (
		p = float64(P)
	)
	if P > 0 {
		div = p * (math.Pow(r+s+10, p-1) + 10*math.Pow(10*(r+s), p-1))
	}
	return
}

type PolyVectorField3 struct{}

func (lpf PolyVectorField3) String() string { return "[s+10,10r]^p Zero Divergence Field" }

func (lpf PolyVectorField3) F(r, s float64, P int) (f1, f2 float64) {
	var (
		p = float64(P)
	)
	f1, f2 = math.Pow(s+10, p), math.Pow(10*r, p)
	return
}

func (lpf PolyVectorField3) Divergence(r, s float64, P int) (div float64) {
	return
}

type PolyVectorField2 struct{}

func (lpf PolyVectorField2) String() string { return "[r,s]^p Simple Field" }

func (lpf PolyVectorField2) F(r, s float64, P int) (f1, f2 float64) {
	var (
		p = float64(P)
	)
	f1, f2 = math.Pow(r, p), math.Pow(s, p)
	return
}

func (lpf PolyVectorField2) Divergence(r, s float64, P int) (div float64) {
	var (
		p = float64(P)
	)
	if P > 0 {
		div = p * (math.Pow(r, p-1) + math.Pow(s, p-1))
	}
	return
}

type PolyScalarField struct{}

func (lpf PolyScalarField) P(r, s float64, P int) (val float64) {
	var (
		p = float64(P)
	)
	val = math.Pow(10*r+s+10, p)
	return
}

func (lpf PolyScalarField) Gradient(r, s float64, P int) (grad [2]float64) {
	var (
		p = float64(P)
	)
	if P > 0 {
		grad = [2]float64{
			10. * p * math.Pow(10*r+s+10, p-1.),
			p * math.Pow(10*r+s+10, p-1.),
		}
	}
	return
}
