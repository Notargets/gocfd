package utils

import (
	"gonum.org/v1/gonum/blas/blas64"
	"gonum.org/v1/gonum/mat"
)

type Vector struct {
	V *mat.VecDense
}

func NewVector(n int, dataO ...[]float64) Vector {
	var (
		data []float64
	)
	if len(dataO) != 0 {
		data = dataO[0]
	}
	return Vector{
		mat.NewVecDense(n, data),
	}
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
func (v Vector) Transpose() Matrix {
	var (
		nr, nc = v.V.Dims()
		m      *mat.Dense
	)
	m = mat.NewDense(nc, nr, v.V.RawVector().Data)
	return Matrix{m}
}

func (v Vector) ToIndex() (I Index) {
	var (
		data = v.V.RawVector().Data
	)
	I = make(Index, v.Len())
	for i, val := range data {
		I[i] = int(val)
	}
	return
}

func (v Vector) ToMatrix() Matrix {
	var (
		m = mat.NewDense(v.V.Len(), 1, v.V.RawVector().Data)
	)
	return Matrix{m}
}

func (v Vector) Subtract(a Vector) Vector {
	var (
		data  = v.V.RawVector().Data
		dataA = a.V.RawVector().Data
	)
	for i := range data {
		data[i] -= dataA[i]
	}
	return v
}

func (v Vector) Subset(I Index) Vector {
	var (
		data  = v.V.RawVector().Data
		dataR = make([]float64, len(I))
		r     *mat.VecDense
	)
	for i, ind := range I {
		dataR[i] = data[ind]
	}
	r = mat.NewVecDense(len(dataR), dataR)
	return Vector{r}
}

func (v Vector) Scale(a float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i := range data {
		data[i] *= a
	}
	return v
}

func (v Vector) AddScalar(a float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i := range data {
		data[i] += a
	}
	return v
}

func (v Vector) Set(a float64) Vector {
	var (
		data = v.V.RawVector().Data
	)
	for i := range data {
		data[i] = a
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
