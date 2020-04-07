package utils

import (
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"
)

type Vector struct {
	V *mat.VecDense
}

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (v Vector) Dims() (r, c int)         { return v.V.Dims() }
func (v Vector) At(i, j int) float64      { return v.V.At(i, j) }
func (v Vector) T() mat.Matrix            { return v.V.T() }
func (v Vector) AtVec(i int) float64      { return v.V.AtVec(i) }
func (v Vector) RawVector() blas64.Vector { return v.V.RawVector() }
func (v Vector) SubVec(a, b Vector)       { v.V.SubVec(a.V, b.V) }
func (v Vector) Len() int                 { return v.V.Len() }

// Chainable (extended) methods
func (v Vector) Sub(a Vector) Vector { v.V.SubVec(v.V, a.V); return v }
func (v Vector) Add(a float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i := range data {
		data[i] += a
	}
	return v
}

func (v Vector) Apply(f func(float64) float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i, val := range data {
		data[i] = f(val)
	}
	return v
}

func (v Vector) POW(p int) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i, val := range data {
		data[i] = POW(val, p)
	}
	return v
}

func (v Vector) Min() (min float64) {
	var (
		data = v.V.RawVector().Data
	)
	min = data[0]
	for _, val := range data {
		if val < min {
			min = val
		}
	}
	return
}

func (v Vector) Max() (max float64) {
	var (
		data = v.V.RawVector().Data
	)
	max = data[0]
	for _, val := range data {
		if val > max {
			max = val
		}
	}
	return
}

// Row is a user-defined Row vector.
type Row []float64

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (v Row) Dims() (r, c int)    { return 1, len(v) }
func (v Row) At(_, j int) float64 { return v[j] }
func (v Row) T() mat.Matrix       { return Column(v) }

// RawVector allows fast path computation with the vector.
func (v Row) RawVector() blas64.Vector {
	return blas64.Vector{N: len(v), Data: v, Inc: 1}
}

// Column is a user-defined Column vector.
type Column []float64

// Dims, At and T minimally satisfy the mat.Matrix interface.
func (v Column) Dims() (r, c int)    { return len(v), 1 }
func (v Column) At(i, _ int) float64 { return v[i] }
func (v Column) T() mat.Matrix       { return Row(v) }

// RawVector allows fast path computation with the vector.
func (v Column) RawVector() blas64.Vector {
	return blas64.Vector{N: len(v), Data: v, Inc: 1}
}
